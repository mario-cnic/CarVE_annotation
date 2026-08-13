# 02_cohort_enrich.R
# Cohort variant table with per-variant OR enrichment.
# Filters master to HiC + both variants and computes ORs vs three comparators:
#   int   — internal cohort controls
#   nfe   — gnomAD v4 NFE (AN ÷ 2, alleles → individuals)
#   total — gnomAD v4 joint total (AN ÷ 2, alleles → individuals)
#
# Rules:
#   - ac_case = 0  → OR = NA  (Haldane correction suppressed; OR uninformative)
#   - NA AN         → positional imputation via impute_an_neighbors()
#   - gnomAD AN     → divided by 2 before use (alleles → individuals)
#
# Usage:   Rscript src/R/02_cohort_enrich.R <GENE>
# Example: Rscript src/R/02_cohort_enrich.R KLHL24
#
# Outputs (data_curated/):
#   {gene}_cohort.rds
#   {gene}_cohort.csv
#
# Per-phenotype OR block (16 cols × 7 phenotypes = 112 OR cols):
#   {PHENO}_ac_case, {PHENO}_an_case,
#   {PHENO}_ac_ctrl_int, {PHENO}_an_ctrl_int,
#   {PHENO}_or_int,   {PHENO}_ci_low_int,   {PHENO}_ci_high_int,   {PHENO}_p_int,
#   {PHENO}_or_nfe,   {PHENO}_ci_low_nfe,   {PHENO}_ci_high_nfe,   {PHENO}_p_nfe,
#   {PHENO}_or_total, {PHENO}_ci_low_total, {PHENO}_ci_high_total, {PHENO}_p_total
#
# gnomAD AC/AN controls appear ONCE (block_gnomad_counts) — not repeated per phenotype.

suppressPackageStartupMessages({
  library(tidyverse)
})
# statsJPO: source the development version (installed package is stale)
# Provides: compute_or, impute_an_neighbors
source("src/R/statsJPO.R")

# ── Parameters ────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript src/R/02_cohort_enrich.R <GENE> <RAW_ROOT>")
gene <- toupper(args[1])
raw_root <- args[2]

out_dir_data <- file.path(raw_root, "data_curated", gene)
dir.create(out_dir_data, showWarnings = FALSE, recursive = TRUE)

master_path <- file.path(out_dir_data, paste0(tolower(gene), "_master.rds"))
out_rds     <- file.path(out_dir_data, paste0(tolower(gene), "_cohort.rds"))
out_csv     <- file.path(out_dir_data, paste0(tolower(gene), "_cohort.csv"))

if (!file.exists(master_path)) stop("Master not found: ", master_path)

# ── Load & filter to cohort variants ─────────────────────────────────────────
message("\n── Loading master ───────────────────────────────────────────")
master <- readRDS(master_path)
message("Master: ", nrow(master), " rows x ", ncol(master), " cols")

# Extract phenotypes dynamically based on available column format
pheno_cols_legacy <- grep("_AC_CA$", names(master), value = TRUE)
pheno_cols_current <- grep("_P_CA$", names(master), value = TRUE)

if (length(pheno_cols_legacy) > 0) {
  schema <- "legacy"
  PHENOTYPES <- gsub("_AC_CA$", "", pheno_cols_legacy)
} else if (length(pheno_cols_current) > 0) {
  schema <- "current"
  PHENOTYPES <- gsub("_P_CA$", "", pheno_cols_current)
} else {
  stop("No phenotype columns found (expected format: {PHENO}_AC_CA or {PHENO}_P_CA)")
}

message("Schema detected: ", schema)
message("Detected phenotypes: ", paste(PHENOTYPES, collapse = ", "))

cohort <- filter(master, SOURCE_ORIGIN %in% c("HiC", "both"))
message("Cohort variants: ", nrow(cohort),
        "  (HiC: ", sum(cohort$SOURCE_ORIGIN == "HiC"),
        " | both: ", sum(cohort$SOURCE_ORIGIN == "both"), ")")

# ── Drop pre-computed OR/P/IC columns from source XLSX ───────────────────────
stale_or_cols <- as.vector(outer(PHENOTYPES, c("_OR", "_P_VALUE", "_IC"), paste0))
n_stale       <- length(intersect(stale_or_cols, names(cohort)))
cohort        <- select(cohort, -any_of(stale_or_cols))
message("Dropped ", n_stale, " pre-computed OR/P/IC columns from XLSX")

# ── Prepare gnomAD controls (alleles → individuals, with AN imputation) ───────
message("\n── Preparing gnomAD controls ────────────────────────────────")

# Positional AN imputation: variants absent from gnomAD have NA AN;
# their locus was still sequenced — estimate from positional neighbors.
an_nfe_raw   <- cohort$gnomADv4_AN_joint_nfe
an_total_raw <- cohort$gnomADv4_AN_joint

an_nfe_ind   <- impute_an_neighbors(cohort$Locus, an_nfe_raw)   / 2
an_total_ind <- impute_an_neighbors(cohort$Locus, an_total_raw) / 2

gnomad_ac_nfe   <- replace_na(cohort$gnomADv4_AC_joint_nfe, 0)
gnomad_ac_total <- replace_na(cohort$gnomADv4_AC_joint,     0)

message("gnomAD NFE   AN imputed: ", sum(is.na(an_nfe_raw)),   " variants")
message("gnomAD total AN imputed: ", sum(is.na(an_total_raw)), " variants")

# ── safe_or helper ────────────────────────────────────────────────────────────
# Returns NA when ac_case = 0 (OR uninformative; Haldane correction suppressed).
# Also returns NA for missing/invalid inputs or zero-denominator AN.
na_or <- tibble(OR = NA_real_, OR_lower = NA_real_,
                OR_upper = NA_real_, p_value = NA_real_)

safe_or <- function(ac_case, an_case, ac_ctrl, an_ctrl) {
  if (any(is.na(c(ac_case, an_case, ac_ctrl, an_ctrl)))) return(na_or)
  if (an_case == 0 || an_ctrl == 0)                       return(na_or)
  if (ac_case == 0)                                        return(na_or)
  tryCatch(
    compute_or(ac_case, an_case, ac_ctrl, an_ctrl) |>
      select(OR, OR_lower, OR_upper, p_value),
    error = function(e) na_or
  )
}

# ── Compute per-variant OR enrichment ─────────────────────────────────────────
# Per-phenotype block:
#   ac_case, an_case, ac_ctrl_int, an_ctrl_int  (counts)
#   or_int / nfe / total + CI + p               (ORs vs 3 comparators)
message("\n── Computing ORs ────────────────────────────────────────────")

or_blocks <- map(PHENOTYPES, function(ph) {
  if (schema == "legacy") {
    ac_case_raw  <- cohort[[paste0(ph, "_AC_CA")]]
    an_case_raw  <- cohort[[paste0(ph, "_AN_CA")]]
    ac_ctrl_raw  <- cohort[[paste0(ph, "_AC_CON")]]
    an_ctrl_raw  <- cohort[[paste0(ph, "_AN_CON")]]
  } else {
    ac_case_raw  <- cohort[[paste0(ph, "_P_CA")]]
    nac_case_raw <- cohort[[paste0(ph, "_NP_CA")]]
    ac_ctrl_raw  <- cohort[[paste0(ph, "_P_CON")]]
    nac_ctrl_raw <- cohort[[paste0(ph, "_NP_CON")]]
    an_case_raw  <- ac_case_raw + nac_case_raw
    an_ctrl_raw  <- ac_ctrl_raw + nac_ctrl_raw
  }

  # NA AC → 0 for OR calculation (variant absent = 0 carriers)
  ac_case <- replace(ac_case_raw, is.na(ac_case_raw), 0)
  ac_ctrl <- replace(ac_ctrl_raw, is.na(ac_ctrl_raw), 0)

  # impute when NA (patient not tested for this phenotype)
  an_case     <- impute_an_neighbors(cohort$Locus, an_case_raw)
  an_ctrl     <- impute_an_neighbors(cohort$Locus, an_ctrl_raw)

  # Use res_* names to avoid tibble()-internal name collision:
  # once or_int=res_int$OR is defined as a column, a bare 'or_int' in subsequent
  # arguments would resolve to that new numeric column, not to res_int.
  res_int   <- pmap(list(ac_case, an_case, ac_ctrl,         an_ctrl),    safe_or) |> bind_rows()
  res_nfe   <- pmap(list(ac_case, an_case, gnomad_ac_nfe,   an_nfe_ind), safe_or) |> bind_rows()
  res_total <- pmap(list(ac_case, an_case, gnomad_ac_total, an_total_ind), safe_or) |> bind_rows()

  n_or <- sum(!is.na(res_int$OR))
  message(sprintf("  %-6s  OR_int: %d / %d computed  (ac_case=0: %d)",
                  ph, n_or, nrow(cohort), sum(ac_case == 0)))

  block <- tibble(
    ac_case       = ac_case,
    an_case       = an_case,
    ac_ctrl_int   = ac_ctrl,
    an_ctrl_int   = an_ctrl,
    or_int        = res_int$OR,
    ci_low_int    = res_int$OR_lower,
    ci_high_int   = res_int$OR_upper,
    p_int         = res_int$p_value,
    or_nfe        = res_nfe$OR,
    ci_low_nfe    = res_nfe$OR_lower,
    ci_high_nfe   = res_nfe$OR_upper,
    p_nfe         = res_nfe$p_value,
    or_total      = res_total$OR,
    ci_low_total  = res_total$OR_lower,
    ci_high_total = res_total$OR_upper,
    p_total       = res_total$p_value
  )
  names(block) <- paste0(ph, "_", names(block))
  block

}) |> bind_cols()

# ── Column ordering ───────────────────────────────────────────────────────────
# Per filtering_order.txt:
#   identity → gnomAD freq → OR blocks (per phenotype) → gnomAD counts →
#   functional scores → splicing → everything else
#
# gnomAD AC/AN controls (NFE and total) appear ONCE here, not per phenotype.

or_col_names <- names(or_blocks)   # 16 × 7 = 112

block_identity <- c(
  "gene",
  "SYMBOL", "Locus", "SOURCE_ORIGIN",
  "Consequence", "VARIANT_CLASS",
  "HGVSc", "HGVSp", "CDNA_NAME",
  "IMPACT", "NEW_IMPACT", "lof_truncating",
  "FILTER",
  "EXON", "INTRON",
  "Protein_position", "Amino_acids", "Existing_variation",
  "ID_MUTACIÓN", "PATOGENICIDAD"
)

block_gnomad_freq <- c(
  "gnomADv4_AF_joint",
  "gnomADv4_AF_grpmax_joint",
  "gnomADv4_grpmax_joint",
  "gnomADv4_faf95_joint",
  "gnomADv4_fafmax_faf95_max_joint"
)

# gnomAD AC/AN once (raw values; ÷2 applied internally for OR calculation)
block_gnomad_counts <- c(
  "gnomADv4_AC_joint",
  "gnomADv4_AN_joint",
  "gnomADv4_AC_joint_nfe",
  "gnomADv4_AN_joint_nfe"
)

block_scores <- c(
  "REVEL_rankscore",         "REVEL_score",
  "AlphaMissense_rankscore", "AlphaMissense_score",
  "BayesDel_addAF_rankscore","BayesDel_addAF_score",
  "MetaRNN_rankscore",       "MetaRNN_score",
  "PHACTboost_rankscore",    "PHACTboost_score"
)

block_splicing <- c(
  "spliceAI_MAX",
  "SPiP_interpretation", "SPiP_prediction",
  "5UTR_annotation",     "5UTR_consequence"
)

# Assemble: ordered blocks first, then everything else
cohort <- bind_cols(cohort, or_blocks)

ordered_cols    <- c(block_identity, block_gnomad_freq, or_col_names,
                     block_gnomad_counts, block_scores, block_splicing)
ordered_present <- ordered_cols[ordered_cols %in% names(cohort)]
rest_cols       <- setdiff(names(cohort), ordered_present)
cohort_out      <- select(cohort, all_of(c(ordered_present, rest_cols)))

# ── Summary ───────────────────────────────────────────────────────────────────
message("\n── Output summary ───────────────────────────────────────────")
message("Dimensions: ", nrow(cohort_out), " rows x ", ncol(cohort_out), " cols")
message("OR cols   : ", length(or_col_names),
        "  (", length(PHENOTYPES), " phenotypes × 16)")
message("Ordered cols present: ", length(ordered_present),
        " / ", length(ordered_cols))

# ── Sanity checks ─────────────────────────────────────────────────────────────
uniquely_stale <- as.vector(outer(PHENOTYPES, c("_P_VALUE", "_IC"), paste0))

stopifnot(
  "Row count"               = nrow(cohort_out) == nrow(cohort),
  "Locus complete"          = all(!is.na(cohort_out$Locus)),
  "OR cols present"         = all(or_col_names %in% names(cohort_out)),
  "Stale P_VALUE/IC absent" = !any(uniquely_stale %in% names(cohort_out)),
  "ac_case=0 → OR_int NA"   = all(is.na(
    cohort_out[[paste0(PHENOTYPES[1], "_or_int")]][cohort_out[[paste0(PHENOTYPES[1], "_ac_case")]] == 0]
  ))
)
message("Sanity checks: PASS")

# ── Export ────────────────────────────────────────────────────────────────────
saveRDS(cohort_out, out_rds)
write_csv(cohort_out, out_csv)
message("\nSaved RDS: ", out_rds)
message("Saved CSV: ", out_csv)
