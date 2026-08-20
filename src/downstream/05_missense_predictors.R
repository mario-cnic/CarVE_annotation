# 05_missense_predictors.R
# Unified script for In-Silico Predictor Profiling (Radiography & Stats)
# Generalizes the KLHL24-specific analysis (91, 92 scripts) to any gene.
#
# Usage:   Rscript src/R/05_missense_predictors.R <GENE> <RAW_ROOT>
# Example: Rscript src/R/05_missense_predictors.R KLHL24 .

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript src/R/05_missense_predictors.R <GENE> <RAW_ROOT>")
gene <- toupper(args[1])
raw_root <- args[2]

MAF_THRESHOLD <- 0.00005

out_dir_data <- file.path(raw_root, "data_curated", gene, "missense_predictors")
out_dir_fig  <- file.path(raw_root, "figures", gene, "missense_predictors")
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

master_path <- file.path(raw_root, "data_curated", gene, paste0(tolower(gene), "_master.rds"))
if (!file.exists(master_path)) stop("Master data not found: ", master_path)

message("\n── Loading master ───────────────────────────────────────────")
master <- readRDS(master_path)

# Extract phenotypes dynamically
pheno_cols_legacy <- grep("_AC_CA$", names(master), value = TRUE)
pheno_cols_current <- grep("_P_CA$", names(master), value = TRUE)

if (length(pheno_cols_legacy) > 0) {
  schema <- "legacy"
  PHENOTYPES <- gsub("_AC_CA$", "", pheno_cols_legacy)
} else if (length(pheno_cols_current) > 0) {
  schema <- "current"
  PHENOTYPES <- gsub("_P_CA$", "", pheno_cols_current)
} else {
  stop("No phenotype columns found.")
}

message("Target Gene: ", gene)
message("Detected phenotypes: ", paste(PHENOTYPES, collapse = ", "))

# 1. Clean & Filter to base missense variants
message("\n── Filtering Missense variants ──────────────────────────────")
df_base <- master %>%
  dplyr::filter(
    str_detect(Consequence, "missense_variant"),
    is.na(FILTER) | FILTER == "PASS" | FILTER == ".",
    is.na(gnomADv4_AF_joint) | gnomADv4_AF_joint <= MAF_THRESHOLD
  ) %>%
  dplyr::mutate(
    # Extract Amino Acid position from HGVSp (e.g., "p.Arg502Gln" -> 502)
    # Handling cases where HGVSp might contain multiple components separated by %3D or others
    HGVSp_clean = sub(".*p\\.", "p.", HGVSp),
    AA_Pos = as.numeric(str_match(HGVSp_clean, "p\\.\\(?([a-zA-Z]+)?([0-9]+)")[,3]),
    
    # Cast scores to numeric
    REVEL = suppressWarnings(as.numeric(REVEL_score)),
    AlphaMissense = suppressWarnings(as.numeric(AlphaMissense_score)),
    BayesDel = suppressWarnings(as.numeric(BayesDel_addAF_score)),
    
    gnomAD_AC_NFE = tidyr::replace_na(as.numeric(gnomADv4_AC_joint_nfe), 0)
  ) %>%
  dplyr::filter(!is.na(AA_Pos))

message("Found ", nrow(df_base), " rare, PASS missense variants with AA positional mapping.")

# 2. Plot: Predictor Radiography
message("\n── Generating Predictor Radiography Plot ────────────────────")
df_plot_multi <- df_base %>%
  dplyr::select(Locus, AA_Pos, REVEL, AlphaMissense, BayesDel) %>%
  tidyr::pivot_longer(
    cols = c(REVEL, AlphaMissense, BayesDel),
    names_to = "Predictor",
    values_to = "Score"
  ) %>%
  dplyr::filter(!is.na(Score))

if (nrow(df_plot_multi) > 0) {
  p_multi <- ggplot(df_plot_multi, aes(x = AA_Pos, y = Score, color = Predictor)) +
    geom_point(alpha = 0.5, size = 2.5) +
    geom_smooth(method = "loess", span = 0.3, se = FALSE, color = "black", linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~ Predictor, scales = "free_y", ncol = 1) + 
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Missense Predictors Radiography -", gene),
      subtitle = "Pathogenicity scores across amino acid sequence",
      x = "Amino Acid Position (N-term -> C-term)",
      y = "Pathogenicity Score"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold", color = "white"),
      strip.background = element_rect(fill = "#2C3E50"),
      panel.grid.minor = element_blank()
    )
    
  out_rad <- file.path(out_dir_fig, paste0("radiography_", tolower(gene), "_predictors.png"))
  ggsave(out_rad, p_multi, width = 9, height = 8, dpi = 150)
  message("  Saved radiography plot: ", out_rad)
} else {
  message("  Not enough data for radiography plot.")
}

# 3. Stats: Casos vs Controles (gnomAD NFE) for each phenotype
message("\n── Comparing Scores (Cases vs gnomAD NFE) ───────────────────")

for (ph in PHENOTYPES) {
  # Get case metric conditionally
  if (schema == "legacy") {
    ac_case_col <- paste0(ph, "_AC_CA")
  } else {
    ac_case_col <- paste0(ph, "_P_CA")
  }
  
  if (!ac_case_col %in% names(df_base)) next
  
  df_casos <- df_base %>% 
    dplyr::mutate(AC_CASE = tidyr::replace_na(as.numeric(!!sym(ac_case_col)), 0)) %>%
    dplyr::filter(AC_CASE > 0) %>% 
    tidyr::uncount(AC_CASE) %>% 
    dplyr::mutate(Grupo = "Cases")
    
  df_controles <- df_base %>% 
    dplyr::filter(gnomAD_AC_NFE > 0) %>% 
    tidyr::uncount(gnomAD_AC_NFE) %>% 
    dplyr::mutate(Grupo = "gnomAD NFE")
    
  df_plot <- dplyr::bind_rows(df_casos, df_controles)
  
  if(nrow(df_plot) == 0) {
    message("  Skipping ", ph, " (no case/control variants).")
    next
  }
  
  df_long <- df_plot %>%
    dplyr::select(Locus, Grupo, REVEL, AlphaMissense, BayesDel) %>%
    tidyr::pivot_longer(cols = c(REVEL, AlphaMissense, BayesDel), names_to = "Predictor", values_to = "Score") %>%
    dplyr::filter(!is.na(Score))
    
  if(nrow(df_long) == 0) next
  
  # Wilcoxon stats
  stats_res <- df_long %>%
    dplyr::group_by(Predictor) %>%
    dplyr::summarise(
      Alleles_Cases = sum(Grupo == "Cases"),
      Alleles_Controls = sum(Grupo == "gnomAD NFE"),
      Median_Cases = median(Score[Grupo == "Cases"], na.rm = TRUE),
      Median_Controls = median(Score[Grupo == "gnomAD NFE"], na.rm = TRUE),
      P_Value = if (Alleles_Cases > 0 & Alleles_Controls > 0) {
                  wilcox.test(Score ~ Grupo, exact = FALSE)$p.value
                } else { NA },
      .groups = "drop"
    )
    
  out_csv <- file.path(out_dir_data, paste0("stats_predictors_", tolower(gene), "_", ph, ".csv"))
  write.csv(stats_res, out_csv, row.names = FALSE)
  message("  Saved stats for ", ph, ": ", out_csv)
  
  # Boxplot
  p_box <- ggplot(df_long, aes(x = Grupo, y = Score, fill = Grupo)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    facet_wrap(~ Predictor, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste("Predictor Comparison -", gene, "-", ph),
      x = "Group", y = "Score"
    ) +
    theme(legend.position = "none")
    
  out_box <- file.path(out_dir_fig, paste0("boxplot_", tolower(gene), "_", ph, "_predictors.png"))
  ggsave(out_box, p_box, width = 7, height = 5, dpi = 150)
  message("  Saved boxplot: ", out_box)
}

message("\nDone.")
