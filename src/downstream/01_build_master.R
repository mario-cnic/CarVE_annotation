# 01_build_master.R
# Build gene-level master table from parsed XLSX.
# Retains ALL source columns (original names) plus two derived fields:
#   gene           — gene identifier (first column)
#   lof_truncating — TRUE for canonical LOF variants passing QC
#
# Usage:   Rscript scripts/01_build_master.R <GENE>
# Example: Rscript scripts/01_build_master.R KLHL24
#
# Outputs (data_curated/):
#   {gene}_master.rds
#   {gene}_master.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(arrow)
})
# statsJPO: source the development version (installed package is stale)
# Provides: compute_or, compute_enrichment, gnomad_pass_unified, impute_an_neighbors
source("src/R/statsJPO.R")

# ── Parameters ────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript src/R/01_build_master.R <GENE>")
gene <- toupper(args[1])
# ── Paths ─────────────────────────────────────────────────────────────────────
raw_root <- args[2]

# DSP and TCAP share a folder
gene_folder <- if (gene %in% c("DSP", "TCAP")) "TCAP_DSP" else gene

input_csv <- file.path(raw_root, paste0(gene, ".parsed.clean.csv"))
input_xlsx <- file.path(raw_root, paste0(gene, ".parsed.clean.xlsx"))
input_pq <- file.path(raw_root, paste0(gene, ".parsed.clean.pq"))

if (file.exists(input_pq)) {
  input_file_path <- input_pq
} else if (file.exists(input_csv)) {
  input_file_path <- input_csv
} else if (file.exists(input_xlsx)) {
  input_file_path <- input_xlsx
} else {
  stop("File not found. Neither .csv, .xlsx nor .pq exists for ", gene, " in ", raw_root)
}

out_dir_data <- file.path(raw_root, "data_curated", gene)
out_dir_res  <- file.path(raw_root, "results",      gene)
dir.create(out_dir_data, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_res,  showWarnings = FALSE, recursive = TRUE)

out_rds  <- file.path(out_dir_data, paste0(tolower(gene), "_master.rds"))
out_csv  <- file.path(out_dir_data, paste0(tolower(gene), "_master.csv"))
qc_path  <- file.path(out_dir_res,  paste0("audit_fields_", gene, ".txt"))

# ── Load ──────────────────────────────────────────────────────────────────────
message("\n── Reading source ───────────────────────────────────────────")
message("Gene:   ", gene)
message("Source: ", input_file_path)

if (grepl("\\.pq$", input_file_path)) {
  raw <- read_parquet(input_file_path)
} else if (grepl("\\.csv$", input_file_path)) {
  raw <- read_csv(input_file_path, show_col_types = FALSE)
} else {
  raw <- read_excel(input_file_path, guess_max = 7000)
}
message("Source: ", nrow(raw), " rows x ", ncol(raw), " cols")

# ── Warn about duplicate column names ─────────────────────────────────────────
# readxl appends a suffix (e.g. .3) to columns that appear more than once.
dup_cols <- names(raw)[grepl("\\.[0-9]+$", names(raw))]
if (length(dup_cols) > 0) {
  message("\nWARNING — ", length(dup_cols),
          " duplicate-flagged column(s) from source (readxl suffix):")
  walk(dup_cols, ~ message("  ", .x))
}

# ── Fix HiC FILTER → "PASS" ───────────────────────────────────────────────────
# HiC variants pass internal QC by definition; the FILTER field is absent/NA
# for variants not present in gnomAD.
master <- raw |>
  mutate(FILTER = if_else(SOURCE_ORIGIN == "HiC", "PASS", FILTER))

# ── Add lof_truncating ────────────────────────────────────────────────────────
# Definition: HIGH VEP impact AND canonical LOF consequence AND passes QC filter.
#   - gnomAD-only variants: must have FILTER == "PASS"
#   - HiC / both: always pass (FILTER already fixed above)
LOF_CONSEQUENCES <- c(
  "frameshift_variant",
  "splice_acceptor_variant",
  "splice_donor_variant",
  "stop_gained"
)

master <- master |>
  mutate(
    lof_truncating = (
      IMPACT == "HIGH" &
      str_detect(Consequence, paste(LOF_CONSEQUENCES, collapse = "|")) &
      (SOURCE_ORIGIN %in% c("HiC", "both") | FILTER == "PASS")
    ) |> replace_na(FALSE)  # NA Consequence → FALSE
  )

n_lof <- sum(master$lof_truncating, na.rm = TRUE)
message("\nlof_truncating: ", n_lof, " TRUE / ", nrow(master), " total")

# ── Add gene identifier ───────────────────────────────────────────────────────
master <- master |> mutate(gene = gene, .before = 1)

# ── Type coercion ─────────────────────────────────────────────────────────────
# Excel occasionally imports numeric fields as character.
NUMERIC_COLS <- c(
  "spliceAI_MAX",
  "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL",
  "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
  "SPiP_prediction",
  "gnomADv4_AF_joint",    "gnomADv4_AF_grpmax_joint",
  "gnomADv4_grpmax_joint",
  "gnomADv4_faf95_joint", "gnomADv4_fafmax_faf95_max_joint",
  "gnomADv4_faf95_joint_nfe",
  "gnomADv4_AC_joint",    "gnomADv4_AN_joint",
  "gnomADv4_AC_joint_nfe","gnomADv4_AN_joint_nfe",
  "gnomADv4_AF_joint_nfe",
  "REVEL_score",       "REVEL_rankscore",
  "AlphaMissense_score","AlphaMissense_rankscore",
  "BayesDel_addAF_score","BayesDel_addAF_rankscore",
  "BayesDel_noAF_score",
  "MetaLR_score",
  "MetaRNN_score",     "MetaRNN_rankscore",
  "MPC_score",
  "VEST4_score",
  "ClinPred_score",
  "PHACTboost_score",  "PHACTboost_rankscore",
  "CADD_PHRED"
)

cols_to_coerce <- intersect(NUMERIC_COLS, names(master))
master <- master |>
  mutate(across(all_of(cols_to_coerce), ~ suppressWarnings(as.numeric(.x))))

message("Coerced ", length(cols_to_coerce), " columns to numeric")

# ── Sanity checks ─────────────────────────────────────────────────────────────
stopifnot(
  "Locus must be complete"         = all(!is.na(master$Locus)),
  "SOURCE_ORIGIN must be complete" = all(!is.na(master$SOURCE_ORIGIN)),
  "SOURCE_ORIGIN values are valid" = all(
    master$SOURCE_ORIGIN %in% c("gnomAD", "HiC", "both")
  )
)

n_na_csq <- sum(is.na(master$Consequence))
if (n_na_csq > 0)
  message("NOTE — ", n_na_csq, " rows have NA Consequence (lof_truncating set to FALSE)")

n_cohort <- sum(master$SOURCE_ORIGIN %in% c("both", "HiC"), na.rm = TRUE)
n_gnomad <- sum(master$SOURCE_ORIGIN == "gnomAD", na.rm = TRUE)

message("\n── Master summary ───────────────────────────────────────────")
message("Dimensions : ", nrow(master), " rows x ", ncol(master), " cols")
message("source_origin — cohort (both/HiC): ", n_cohort,
        " | gnomAD only: ", n_gnomad)

# ── Export ────────────────────────────────────────────────────────────────────
saveRDS(master, out_rds)
write_csv(master, out_csv)
message("\nSaved RDS: ", out_rds)
message("Saved CSV: ", out_csv)

# ── QC report ─────────────────────────────────────────────────────────────────
na_pct <- master |>
  summarise(across(everything(), ~ round(mean(is.na(.x)) * 100, 1))) |>
  pivot_longer(everything(), names_to = "column", values_to = "pct_na") |>
  arrange(desc(pct_na))

qc_lines <- c(
  paste0("QC report — ", gene, "  (", Sys.time(), ")"),
  paste0("Dimensions : ", nrow(master), " rows x ", ncol(master), " cols"),
  paste0("Source file: ", input_file_path),
  "",
  "source_origin counts:",
  paste0("  gnomAD : ", sum(master$SOURCE_ORIGIN == "gnomAD")),
  paste0("  HiC    : ", sum(master$SOURCE_ORIGIN == "HiC")),
  paste0("  both   : ", sum(master$SOURCE_ORIGIN == "both")),
  "",
  "lof_truncating:",
  paste0("  TRUE  : ", sum(master$lof_truncating,  na.rm = TRUE)),
  paste0("  FALSE : ", sum(!master$lof_truncating, na.rm = TRUE)),
  "",
  if (length(dup_cols) > 0)
    c(paste0("Duplicate-flagged cols (", length(dup_cols), "):"),
      paste0("  ", dup_cols))
  else
    "Duplicate cols: none",
  "",
  "% NA per column (top 40 most incomplete):",
  paste0(sprintf("  %-40s %5.1f%%", na_pct$column[1:min(40, nrow(na_pct))],
                 na_pct$pct_na[1:min(40, nrow(na_pct))]))
)

writeLines(qc_lines, qc_path)
message("QC report: ", qc_path)
