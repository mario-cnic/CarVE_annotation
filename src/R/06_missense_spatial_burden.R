# 06_missense_spatial_burden.R
# Unified script for Sliding Window Spatial Burden
# Generalizes KLHL24 specific spatial analysis (93, 94, 95 scripts)
# 
# Usage: Rscript src/R/06_missense_spatial_burden.R <GENE> <RAW_ROOT> [WINDOW_SIZE] [STEP_SIZE] [MAF_THRESHOLD] [REVEL_THRESHOLD]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})
source("src/R/statsJPO.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript src/R/06_missense_spatial_burden.R <GENE> <RAW_ROOT> [WINDOW_SIZE] [STEP_SIZE] [MAF_THRESHOLD] [REVEL_THRESHOLD]")
gene <- toupper(args[1])
raw_root <- args[2]

# Optional Parameters with defaults
window_size <- if (length(args) >= 3) as.numeric(args[3]) else 50
step_size <- if (length(args) >= 4) as.numeric(args[4]) else ceiling(window_size / 2)
MAF_THRESHOLD <- if (length(args) >= 5) as.numeric(args[5]) else 0.00005
revel_threshold <- if (length(args) >= 6) as.numeric(args[6]) else 0.6
thresholds <- list("Todas_Missense"=0, "REVEL_Alto"=revel_threshold)

out_dir_data <- file.path(raw_root, "data_curated", gene, "spatial_burden")
out_dir_fig  <- file.path(raw_root, "figures", gene, "spatial_burden")
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

master_path <- file.path(raw_root, "data_curated", gene, paste0(tolower(gene), "_master.rds"))
if (!file.exists(master_path)) stop("Master data not found: ", master_path)

message("\n── Loading master ───────────────────────────────────────────")
master <- readRDS(master_path)

# Phenotype detection
pheno_cols_legacy <- grep("_AC_CA$", names(master), value = TRUE)
pheno_cols_current <- grep("_P_CA$", names(master), value = TRUE)
if (length(pheno_cols_legacy) > 0) {
  schema <- "legacy"
  PHENOTYPES <- gsub("_AC_CA$", "", pheno_cols_legacy)
} else if (length(pheno_cols_current) > 0) {
  schema <- "current"
  PHENOTYPES <- gsub("_P_CA$", "", pheno_cols_current)
} else stop("No phenotype columns found")

# Base Filtering and Mapping
df_base <- master %>%
  dplyr::filter(
    str_detect(Consequence, "missense_variant"),
    is.na(FILTER) | FILTER == "PASS" | FILTER == ".",
    is.na(gnomADv4_AF_joint) | gnomADv4_AF_joint <= MAF_THRESHOLD
  ) %>%
  dplyr::mutate(
    HGVSp_clean = sub(".*p\\.", "p.", HGVSp),
    AA_Pos = as.numeric(str_match(HGVSp_clean, "p\\.\\(?([a-zA-Z]+)?([0-9]+)")[,3]),
    REVEL = suppressWarnings(as.numeric(REVEL_score)),
    gnomAD_AC_NFE = tidyr::replace_na(as.numeric(gnomADv4_AC_joint_nfe), 0),
    gnomAD_AN_NFE = as.numeric(gnomADv4_AN_joint_nfe)
  ) %>%
  dplyr::filter(!is.na(AA_Pos))

max_aa_len <- suppressWarnings(max(df_base$AA_Pos, na.rm = TRUE))
if(is.infinite(max_aa_len)) stop("No valid AA positions found.")

# Safe Stats Calculation Helper (with Pseudocounts/Haldane)
safe_burden_or <- function(ac_case, an_case, ac_ctrl, an_ctrl) {
  if (is.na(an_ctrl) || an_ctrl <= 0 || is.na(an_case) || an_case <= 0) {
    return(data.frame(OR=NA, CI_Lo=NA, CI_Hi=NA, P_val=NA))
  }
  
  ref_case <- max(0, an_case - ac_case)
  ref_ctrl <- max(0, an_ctrl - ac_ctrl)
  
  if (ac_case == 0 || ac_ctrl == 0) {
    ratio_size <- an_ctrl / an_case
    a <- ac_case + 0.5; b <- ref_case + 0.5
    c_val <- ac_ctrl + (0.5 * ratio_size); d_val <- ref_ctrl + (0.5 * ratio_size)
  } else {
    a <- ac_case; b <- ref_case; c_val <- ac_ctrl; d_val <- ref_ctrl
  }
  
  denom <- b * c_val
  if (denom == 0) return(data.frame(OR=NA, CI_Lo=NA, CI_Hi=NA, P_val=NA))
  
  or_val <- (a * d_val) / denom
  se_log_or <- sqrt(1/a + 1/b + 1/c_val + 1/d_val)
  ci_low <- exp(log(or_val) - 1.96 * se_log_or)
  ci_up  <- exp(log(or_val) + 1.96 * se_log_or)
  
  m <- matrix(c(round(a), round(c_val), round(b), round(d_val)), nrow=2)
  p_val <- tryCatch(fisher.test(m)$p.value, error=function(e) NA)
  
  return(data.frame(OR=or_val, CI_Lo=ci_low, CI_Hi=ci_up, P_val=p_val))
}

# Parameters
window_size <- 50
step_size <- ceiling(window_size / 2)
thresholds <- list("Todas_Missense"=0, "REVEL_Alto"=0.6)

message("\n── Generating spatial burdens (window=", window_size, " step=", step_size, ") ──")

for(ph in PHENOTYPES) {
  if (schema == "legacy") {
    ac_col <- paste0(ph, "_AC_CA")
    an_col <- paste0(ph, "_AN_CA")
  } else {
    ac_col <- paste0(ph, "_P_CA")
    an_col <- "MCH_AN_CA" # Fallback if _NP_CA sum is not straightforward here
  }
  
  if(!ac_col %in% names(df_base)) next
  
  # Impute global AN for cases/controls for the spatial test dynamically
  df_sub <- df_base %>%
    dplyr::mutate(AC_CASE = tidyr::replace_na(as.numeric(!!sym(ac_col)), 0))
    
  if (schema == "current" && paste0(ph, "_NP_CA") %in% names(df_sub)) {
     df_sub <- df_sub %>% dplyr::mutate(AN_CASE = AC_CASE + as.numeric(!!sym(paste0(ph, "_NP_CA"))))
  } else if (an_col %in% names(df_sub)) {
     df_sub <- df_sub %>% dplyr::mutate(AN_CASE = as.numeric(!!sym(an_col)))
  } else {
     df_sub <- df_sub %>% dplyr::mutate(AN_CASE = max(AC_CASE, na.rm=TRUE)) # Rough fallback
  }

  AN_casos <- max(df_sub$AN_CASE, na.rm=TRUE)
  valid_nfe_an <- df_sub$gnomAD_AN_NFE[df_sub$gnomAD_AN_NFE > 0 & !is.na(df_sub$gnomAD_AN_NFE)]
  AN_gnomad <- if(length(valid_nfe_an) > 0) mean(valid_nfe_an, na.rm=TRUE)/2 else 0 # Divided by 2 for ind
  
  if (is.infinite(AN_casos) || AN_gnomad == 0) {
    message("  Skipping ", ph, " (Invalid AN vectors).")
    next
  }

  results_espacial <- list()
  
  for(t_name in names(thresholds)) {
    min_revel <- thresholds[[t_name]]
    df_t <- df_sub %>% dplyr::filter(REVEL >= min_revel | is.na(REVEL) & min_revel == 0)
    
    for(start_pos in seq(1, max_aa_len, by=step_size)) {
      end_pos <- min(start_pos + window_size - 1, max_aa_len)
      
      df_w <- df_t %>% dplyr::filter(AA_Pos >= start_pos & AA_Pos <= end_pos)
      
      ac_case_w <- sum(df_w$AC_CASE, na.rm=TRUE)
      ac_ctrl_w <- sum(df_w$gnomAD_AC_NFE, na.rm=TRUE)
      
      stats <- safe_burden_or(ac_case_w, AN_casos, ac_ctrl_w, AN_gnomad)
      
      results_espacial[[length(results_espacial)+1]] <- data.frame(
        Phenotype = ph,
        Class = t_name,
        Start = start_pos,
        End = end_pos,
        AC_Case = ac_case_w,
        AC_Ctrl = ac_ctrl_w,
        OR = stats$OR,
        CI_Lo = stats$CI_Lo,
        CI_Hi = stats$CI_Hi,
        P_val = stats$P_val
      )
    }
  }
  
  res_df <- bind_rows(results_espacial) %>%
    mutate(LogP = -log10(P_val))
  
  out_csv <- file.path(out_dir_data, paste0("spatial_burden_", tolower(gene), "_", ph, ".csv"))
  write.csv(res_df, out_csv, row.names=FALSE)
  
  # Spatial Plot (Generic Bar Plot)
  p_spatial <- ggplot(res_df, aes(x=(Start+End)/2, y=LogP, fill=OR)) +
    geom_bar(stat="identity", width=step_size) +
    facet_wrap(~Class, ncol=1) +
    geom_hline(yintercept=-log10(0.05), color="red", linetype="dashed") +
    theme_bw() +
    scale_fill_gradient2(low="lightblue", high="darkred", mid="white", midpoint=1, na.value="grey") +
    labs(
      title = paste("Spatial Burden Analysis -", gene, "-", ph),
      x = "Amino Acid Position",
      y = "-log10(P-Value)"
    )
    
  out_plt <- file.path(out_dir_fig, paste0("spatial_burden_", tolower(gene), "_", ph, ".png"))
  ggsave(out_plt, p_spatial, width=10, height=6, dpi=150)
  message("  Saved spatial burden for ", ph, ": ", out_plt)
}

message("\nDone.")
