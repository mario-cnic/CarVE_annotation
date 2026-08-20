# 03_burden.R
# Burden test: aggregate AC by functional class × phenotype, compute OR.
#
# Usage:   Rscript src/R/03_burden.R <GENE>
# Example: Rscript src/R/03_burden.R KLHL24
#
# Outputs (results/):
#   burden_{gene}.rds
#   burden_{gene}.csv
#
# ── Pre-burden filters ────────────────────────────────────────────────────────
#   MAF < 0.00005  (gnomADv4_AF_joint; NA = absent from gnomAD → kept)
#   gnomAD PASS    — all HiC/both variants pass by definition (set in 01)
#
# ── Variant classes (independent, overlapping) ────────────────────────────────
#   truncating    — lof_truncating == TRUE
#   missense      — Consequence contains "missense_variant" (any REVEL)
#   missense_high — missense_variant AND REVEL_score > 0.6
#   synonymous    — Consequence contains "synonymous_variant"
#
# ── gnomAD control: ALL variants in class (cohort + gnomAD-only) ──────────────
#   gnomAD-only variants contribute 0 to case AC but their full gnomAD AC to
#   the control — required for an unbiased gene-level burden comparison.
#
# ── OR logic ──────────────────────────────────────────────────────────────────
#   ac_case = 0 (aggregate) → OR = NA  (same rule as 02)
#   Three comparators per class × phenotype: internal, gnomAD NFE, gnomAD total
#   gnomAD AN: positional imputation + ÷2 (alleles → individuals); rounded
#   gnomAD AN denominator: median across ALL class variants after MAF filter
#   Internal AN: max(AC + NAC) across cohort class variants per phenotype

suppressPackageStartupMessages(library(tidyverse))
# statsJPO: source the development version (installed package is stale)
source("src/R/statsJPO.R")

# ── Parameters ────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript src/R/03_burden.R <GENE> <RAW_ROOT>")
gene <- toupper(args[1])
raw_root <- args[2]
MAF_THRESHOLD <- 0.00005

out_dir_data <- file.path(raw_root, "data_curated", gene)
out_dir_res  <- file.path(raw_root, "results",      gene)
dir.create(out_dir_res, showWarnings = FALSE, recursive = TRUE)

master_path <- file.path(out_dir_data, paste0(tolower(gene), "_master.rds"))
out_rds     <- file.path(out_dir_res,  paste0("burden_", tolower(gene), ".rds"))
out_csv     <- file.path(out_dir_res,  paste0("burden_", tolower(gene), ".csv"))

if (!file.exists(master_path)) stop("Master not found: ", master_path)

# ── Load & filter ─────────────────────────────────────────────────────────────
message("\n── Loading master ───────────────────────────────────────────")
master <- readRDS(master_path)
message("Master: ", nrow(master), " rows x ", ncol(master), " cols")

# Extract phenotypes dynamically based on available column format for 03 burden test
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

# all_filt: ALL variants passing MAF (cohort + gnomAD-only) — used for gnomAD ctrl
all_filt <- master |>
  filter(is.na(gnomADv4_AF_joint) | gnomADv4_AF_joint < MAF_THRESHOLD)

# cohort: HiC + both only — used for case AC and internal ctrl
cohort <- all_filt |>
  filter(SOURCE_ORIGIN %in% c("HiC", "both"))

message("After MAF < ", MAF_THRESHOLD, " filter:")
message("  All variants : ", nrow(all_filt),
        "  (removed: ", nrow(master) - nrow(all_filt), ")")
message("  Cohort       : ", nrow(cohort),
        "  (HiC: ", sum(cohort$SOURCE_ORIGIN == "HiC"),
        " | both: ", sum(cohort$SOURCE_ORIGIN == "both"), ")")
message("  gnomAD-only  : ", nrow(all_filt) - nrow(cohort))

# ── gnomAD AN: positional imputation across ALL filtered variants ─────────────
message("\n── Imputing gnomAD ANs ──────────────────────────────────────")
an_total_raw <- all_filt$gnomADv4_AN_joint
an_nfe_raw   <- all_filt$gnomADv4_AN_joint_nfe
an_total_imp <- impute_an_neighbors(all_filt$Locus, an_total_raw)
an_nfe_imp   <- impute_an_neighbors(all_filt$Locus, an_nfe_raw)
message("gnomAD total AN imputed: ", sum(is.na(an_total_raw)), " variants")
message("gnomAD NFE   AN imputed: ", sum(is.na(an_nfe_raw)),   " variants")

# ── Variant class indices ─────────────────────────────────────────────────────
# idx_all: in all_filt (gnomAD AC denominator/numerator)
# idx_coh: in cohort   (case AC + internal ctrl)
define_classes <- function(df) list(
  truncating    = which(df$lof_truncating == TRUE),
  missense      = which(str_detect(df$Consequence, "missense_variant")),
  missense_high = which(str_detect(df$Consequence, "missense_variant") &
                          !is.na(df$REVEL_score) & df$REVEL_score > 0.6),
  synonymous    = which(str_detect(df$Consequence, "synonymous_variant"))
)

CLASSES_ALL <- define_classes(all_filt)
CLASSES_COH <- define_classes(cohort)

message("\n── Variant class sizes ──────────────────────────────────────")
walk(names(CLASSES_ALL), function(cl)
  message(sprintf("  %-14s  all: %3d  (cohort: %d  gnomAD-only: %d)",
                  cl,
                  length(CLASSES_ALL[[cl]]),
                  length(CLASSES_COH[[cl]]),
                  length(CLASSES_ALL[[cl]]) - length(CLASSES_COH[[cl]]))))

# ── safe_burden_or ────────────────────────────────────────────────────────────
na_or <- tibble(
  OR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_,
  p_value = NA_real_, haldane_applied = NA
)

safe_burden_or <- function(ac_case, an_case, ac_ctrl, an_ctrl) {
  if (any(is.na(c(ac_case, an_case, ac_ctrl, an_ctrl)))) return(na_or)
  if (an_case == 0 || an_ctrl == 0)                       return(na_or)
  if (ac_case == 0)                                        return(na_or)
  tryCatch(
    compute_or(ac_case, an_case, ac_ctrl, an_ctrl) |>
      select(OR, OR_lower, OR_upper, p_value, haldane_applied),
    error = function(e) na_or
  )
}

# ── Burden calculation ────────────────────────────────────────────────────────
message("\n── Calculating burden ───────────────────────────────────────")

results <- map_dfr(names(CLASSES_ALL), function(class_name) {
  idx_all <- CLASSES_ALL[[class_name]]
  idx_coh <- CLASSES_COH[[class_name]]
  n_all   <- length(idx_all)
  n_coh   <- length(idx_coh)

  # gnomAD control: ALL class variants (cohort + gnomAD-only)
  # AN denominator: median imputed AN / 2, rounded to integer
  an_gnomad_total <- round(median(an_total_imp[idx_all], na.rm = TRUE) / 2)
  an_gnomad_nfe   <- round(median(an_nfe_imp[idx_all],   na.rm = TRUE) / 2)
  ac_gnomad_total <- sum(replace_na(all_filt$gnomADv4_AC_joint[idx_all],     0))
  ac_gnomad_nfe   <- sum(replace_na(all_filt$gnomADv4_AC_joint_nfe[idx_all], 0))

  map_dfr(PHENOTYPES, function(ph) {
    if (schema == "legacy") {
      ac_case_vec  <- replace_na(cohort[[paste0(ph, "_AC_CA")]][idx_coh], 0)
      an_case_vec  <- cohort[[paste0(ph, "_AN_CA")]][idx_coh]
      ac_ctrl_vec  <- replace_na(cohort[[paste0(ph, "_AC_CON")]][idx_coh], 0)
      an_ctrl_vec  <- cohort[[paste0(ph, "_AN_CON")]][idx_coh]

      ac_case <- sum(ac_case_vec)
      ac_ctrl <- sum(ac_ctrl_vec)

      # AN limit: using max across vector values
      an_case <- max(an_case_vec, na.rm = TRUE)
      an_ctrl <- max(an_ctrl_vec, na.rm = TRUE)

    } else {
      # Case and internal ctrl: cohort variants only (current schema)
      ac_case_vec  <- replace_na(cohort[[paste0(ph, "_P_CA")]][idx_coh],  0)
      nac_case_vec <- cohort[[paste0(ph, "_NP_CA")]][idx_coh]
      ac_ctrl_vec  <- replace_na(cohort[[paste0(ph, "_P_CON")]][idx_coh], 0)
      nac_ctrl_vec <- cohort[[paste0(ph, "_NP_CON")]][idx_coh]

      ac_case <- sum(ac_case_vec)
      ac_ctrl <- sum(ac_ctrl_vec)

      # AN: cohort size per phenotype; constant → max() is safe
      an_case <- max(ac_case_vec + replace_na(nac_case_vec, 0), na.rm = TRUE)
      an_ctrl <- max(ac_ctrl_vec + replace_na(nac_ctrl_vec, 0), na.rm = TRUE)
    }

    if (is.infinite(an_case)) an_case <- NA_real_
    if (is.infinite(an_ctrl)) an_ctrl <- NA_real_

    res_int   <- safe_burden_or(ac_case, an_case, ac_ctrl,         an_ctrl)
    res_nfe   <- safe_burden_or(ac_case, an_case, ac_gnomad_nfe,   an_gnomad_nfe)
    res_total <- safe_burden_or(ac_case, an_case, ac_gnomad_total, an_gnomad_total)

    tibble(
      gene              = gene,
      variant_class     = class_name,
      phenotype         = ph,
      n_variants_cohort = n_coh,
      n_variants_all    = n_all,   # cohort + gnomAD-only

      ac_case           = ac_case,
      an_case           = an_case,
      ac_exceeds_an     = ac_case > an_case,

      ac_ctrl_int       = ac_ctrl,
      an_ctrl_int       = an_ctrl,
      or_int            = res_int$OR,
      ci_low_int        = res_int$OR_lower,
      ci_high_int       = res_int$OR_upper,
      p_int             = res_int$p_value,
      haldane_int       = res_int$haldane_applied,

      ac_ctrl_nfe       = ac_gnomad_nfe,
      an_ctrl_nfe       = an_gnomad_nfe,
      or_nfe            = res_nfe$OR,
      ci_low_nfe        = res_nfe$OR_lower,
      ci_high_nfe       = res_nfe$OR_upper,
      p_nfe             = res_nfe$p_value,
      haldane_nfe       = res_nfe$haldane_applied,

      ac_ctrl_total     = ac_gnomad_total,
      an_ctrl_total     = an_gnomad_total,
      or_total          = res_total$OR,
      ci_low_total      = res_total$OR_lower,
      ci_high_total     = res_total$OR_upper,
      p_total           = res_total$p_value,
      haldane_total     = res_total$haldane_applied
    )
  })
})

# ── Summary ───────────────────────────────────────────────────────────────────
message("\n── Output: ", nrow(results), " rows  (",
        n_distinct(results$variant_class), " classes × ",
        n_distinct(results$phenotype), " phenotypes)")

or_summary <- results |>
  group_by(variant_class) |>
  summarise(
    n_all    = first(n_variants_all),
    n_cohort = first(n_variants_cohort),
    ac_gnomad_nfe = first(ac_ctrl_nfe),
    n_or_int  = sum(!is.na(or_int)),
    or_range  = if (any(!is.na(or_int)))
      paste0(round(min(or_int, na.rm=TRUE), 2), "–", round(max(or_int, na.rm=TRUE), 2))
    else "all NA",
    .groups = "drop"
  )

message("\n  class          n_all cohort  ac_nfe  ORs  range_int")
walk(seq_len(nrow(or_summary)), function(i)
  message(sprintf("  %-14s  %3d    %3d     %3d   %d/7  %s",
                  or_summary$variant_class[i],
                  or_summary$n_all[i],
                  or_summary$n_cohort[i],
                  or_summary$ac_gnomad_nfe[i],
                  or_summary$n_or_int[i],
                  or_summary$or_range[i])))

# ── Export ────────────────────────────────────────────────────────────────────
saveRDS(results, out_rds)
write_csv(results, out_csv)
message("\nSaved RDS: ", out_rds)
message("Saved CSV: ", out_csv)
