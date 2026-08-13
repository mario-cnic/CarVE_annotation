# 07_domain_burden.R
# Script for Fixed Intervals (Domain) Spatial Burden
# Calculates burden using predefined biological domains/regions natively (like DSP)
# 
# Usage: Rscript src/R/07_domain_burden.R <GENE> <RAW_ROOT> [MAF_THRESHOLD] [REVEL_THRESHOLD]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

# Load external stats functions if available, else define inline
if(file.exists("src/R/statsJPO.R")) {
  source("src/R/statsJPO.R")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript src/R/07_domain_burden.R <GENE> <RAW_ROOT> [MAF_THRESHOLD] [REVEL_THRESHOLD]")
gene <- toupper(args[1])
raw_root <- args[2]

# Optional Parameters with defaults
MAF_THRESHOLD <- if (length(args) >= 3) as.numeric(args[3]) else 0.00005
revel_threshold <- if (length(args) >= 4) as.numeric(args[4]) else 0.6

# Define predictor thresholds
predictors <- list(
  "All_missense" = list(col = "ALL", val = 0),
  "REVEL"        = list(col = "REVEL", val = revel_threshold),
  "CADD"         = list(col = "CADD", val = 20),
  "AlphaMissense"= list(col = "AlphaMissense", val = 0.5)
)
names(predictors) <- c("All missense", paste0("REVEL \u2265 ", revel_threshold), "CADD \u2265 20", "AlphaMissense \u2265 0.5")

out_dir_data <- file.path(raw_root, "data_curated", gene, "domain_burden")
out_dir_fig  <- file.path(raw_root, "figures", gene, "domain_burden")
dir.create(out_dir_data, showWarnings = FALSE, recursive = TRUE)
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
    CADD = suppressWarnings(as.numeric(CADD_phred)),
    AlphaMissense = suppressWarnings(as.numeric(AlphaMissense_score)),
    gnomAD_AC_NFE = tidyr::replace_na(as.numeric(gnomADv4_AC_joint_nfe), 0),
    gnomAD_AN_NFE = as.numeric(gnomADv4_AN_joint_nfe)
  ) %>%
  dplyr::filter(!is.na(AA_Pos))

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

# ── Define Domains (Hardcoded from Image provided for DSP or generic fallback) ──
if (gene == "DSP") {
  domain_df <- data.frame(
    Level = c("Major region", "Subregion", "Subregion", "Subregion", "Subregion", "Subregion", "Subregion", "Subregion", "Subregion", "Subregion", "Subregion", "Subregion", "Major region", "Major region", "Subregion", "Subregion", "Subregion", "Subregion", "Subregion", "Subregion", "Subregion"),
    RegionName = c("Head domain", "N-terminal", "SR3", "SR4", "SR5 (N-terminal)", "SH3", "SR5 (C-terminal)", "SR6", "Short arm of plakin domain", "SR7", "SR8", "CT domain", "Rod domain", "Tail domain", "PRD-A", "PRD-B", "LM", "Spacer between linker and PRD-C", "PRD-C", "GSR motif", "Pergola hotspot core within GSR"),
    StartPos = c(1, 1, 179, 272, 376, 461, 517, 546, 655, 655, 771, 883, 1028, 1960, 1960, 2208, 2456, 2566, 2609, 2822, 2828),
    EndPos = c(1028, 178, 271, 375, 460, 516, 545, 654, 882, 770, 882, 1028, 1945, 2871, 2208, 2456, 2565, 2608, 2822, 2871, 2844)
  )
} else {
  message("WARNING: No predefined intervals for gene ", gene, " available natively. Using dummy ranges.")
  max_aa <- suppressWarnings(max(df_base$AA_Pos, na.rm = TRUE))
  domain_df <- data.frame(
    Level = c("Major region", "Major region"),
    RegionName = c("N-term Half", "C-term Half"),
    StartPos = c(1, floor(max_aa/2)),
    EndPos = c(floor(max_aa/2 - 1), max_aa)
  )
}
domain_df$Length <- domain_df$EndPos - domain_df$StartPos + 1

message("\n── Generating domain burdens ──")

for(ph in PHENOTYPES) {
  if (schema == "legacy") {
    ac_col <- paste0(ph, "_AC_CA")
    an_col <- paste0(ph, "_AN_CA")
  } else {
    ac_col <- paste0(ph, "_P_CA")
    an_col <- "MCH_AN_CA" 
  }
  
  if(!ac_col %in% names(df_base)) next
  
  df_sub <- df_base %>%
    dplyr::mutate(AC_CASE = tidyr::replace_na(as.numeric(!!sym(ac_col)), 0))
    
  if (schema == "current" && paste0(ph, "_NP_CA") %in% names(df_sub)) {
     df_sub <- df_sub %>% dplyr::mutate(AN_CASE = AC_CASE + as.numeric(!!sym(paste0(ph, "_NP_CA"))))
  } else if (an_col %in% names(df_sub)) {
     df_sub <- df_sub %>% dplyr::mutate(AN_CASE = as.numeric(!!sym(an_col)))
  } else {
     df_sub <- df_sub %>% dplyr::mutate(AN_CASE = max(AC_CASE, na.rm=TRUE))
  }

  AN_casos <- max(df_sub$AN_CASE, na.rm=TRUE)
  valid_nfe_an <- df_sub$gnomAD_AN_NFE[df_sub$gnomAD_AN_NFE > 0 & !is.na(df_sub$gnomAD_AN_NFE)]
  AN_gnomad <- if(length(valid_nfe_an) > 0) mean(valid_nfe_an, na.rm=TRUE)/2 else 0
  
  if (is.infinite(AN_casos) || AN_gnomad == 0) {
    message("  Skipping ", ph, " (Invalid AN vectors).")
    next
  }

  domain_results <- list()
  
  for(t_name in names(predictors)) {
    pred_col <- predictors[[t_name]]$col
    pred_val <- predictors[[t_name]]$val
    
    if(pred_col == "ALL") {
      df_t <- df_sub
    } else {
      df_t <- df_sub %>% dplyr::filter(!is.na(!!sym(pred_col)) & !!sym(pred_col) >= pred_val)
    }
    
    class_tot_case <- sum(df_t$AC_CASE, na.rm=TRUE)
    class_tot_ctrl <- sum(df_t$gnomAD_AC_NFE, na.rm=TRUE)
    class_label <- paste0(t_name, "\n(Var. Cases = ", class_tot_case, ", Var. Ctrls = ", class_tot_ctrl, ")")
    
    for(i in 1:nrow(domain_df)) {
      dom <- domain_df[i, ]
      
      df_w <- df_t %>% dplyr::filter(AA_Pos >= dom$StartPos & AA_Pos <= dom$EndPos)
      
      ac_case_w <- sum(df_w$AC_CASE, na.rm=TRUE)
      ac_ctrl_w <- sum(df_w$gnomAD_AC_NFE, na.rm=TRUE)
      
      stats <- safe_burden_or(ac_case_w, AN_casos, ac_ctrl_w, AN_gnomad)
      
      domain_results[[length(domain_results)+1]] <- data.frame(
        Phenotype = ph,
        Class = class_label,
        Level = dom$Level,
        Region = dom$RegionName,
        Start = dom$StartPos,
        End = dom$EndPos,
        Length = dom$Length,
        AC_Case = ac_case_w,
        AC_Ctrl = ac_ctrl_w,
        OR = stats$OR,
        CI_Lo = stats$CI_Lo,
        CI_Hi = stats$CI_Hi,
        P_val = stats$P_val
      )
    }
  }
  
  res_df <- bind_rows(domain_results) %>%
    group_by(Class) %>%
    mutate(
      FDR = p.adjust(P_val, method="fdr"),
      LogP = -log10(P_val),
      Significance = case_when(
        FDR < 0.05 ~ "FDR < 0.05",
        P_val < 0.05 ~ "P < 0.05",
        TRUE ~ "Not Sig"
      )
    ) %>%
    ungroup()
  
  out_csv <- file.path(out_dir_data, paste0("domain_burden_", tolower(gene), "_", ph, ".csv"))
  write.csv(res_df, out_csv, row.names=FALSE)
  
  # Fancy Plot 1: Forest Plot of Extracted Domains
  res_subregions <- res_df %>% filter(Level == "Subregion")
  
  p_forest <- ggplot(res_subregions, aes(x = OR, y = reorder(Region, Start))) +
    geom_vline(xintercept = 1, linetype="dashed", color="grey50") +
    geom_errorbar(aes(xmin = CI_Lo, xmax = CI_Hi), width=0.2, color="darkblue") +
    geom_point(aes(color=Significance, size=-log10(P_val))) +
    facet_wrap(~Class, scales="free_x") +
    scale_x_log10() +
    theme_bw() +
    scale_color_manual(values=c("FDR < 0.05"="red", "P < 0.05"="orange", "Not Sig"="grey40")) +
    labs(
      title = paste("Subregion Burden Forest Plot -", gene, "-", ph),
      x = "Odds Ratio (log scale)",
      y = "Domain"
    )
    
  out_forest <- file.path(out_dir_fig, paste0("domain_forest_", tolower(gene), "_", ph, ".png"))
  ggsave(out_forest, p_forest, width=12, height=10, dpi=300, bg="white")
  
  # Fancy Plot 2: Linear Gene Map
  plot_df <- res_df %>%
    mutate(
      ymin = ifelse(Level == "Major region", 1.2, 0),
      ymax = ifelse(Level == "Major region", 2.0, 1.0)
    )

  p_map <- ggplot(plot_df) +
    geom_rect(aes(xmin = Start, xmax = End, ymin = ymin, ymax = ymax, fill = LogP), color = "black", alpha = 0.9) +
    facet_wrap(~ Class, ncol = 1) +
    scale_fill_gradient2(low = "lightblue", high = "darkred", mid = "white", midpoint = 1.3, na.value = "grey") +
    geom_text_repel(
      aes(x = (Start + End) / 2, 
          y = ifelse(Level == "Major region", ymax, ymin), 
          label = Region),
      direction = "y",
      nudge_y = ifelse(plot_df$Level == "Major region", 1.5, -1.5),
      hjust = 0.5,
      vjust = 0.5,
      size = 3.5,
      segment.size = 0.3,
      min.segment.length = 0,
      max.overlaps = Inf
    ) +
    scale_y_continuous(limits = c(-3, 4)) +
    theme_bw() +
    theme(
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white"),
      axis.text.y = element_blank(), 
      axis.title.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey95")
    ) +
    labs(
      title = paste("Domain Gene Map Burden -", gene, "-", ph),
      x = "Amino Acid Position",
      fill = "-log10(P)"
    )
    
  out_map <- file.path(out_dir_fig, paste0("domain_genemap_", tolower(gene), "_", ph, ".png"))
  ggsave(out_map, p_map, width=12, height=8, dpi=300, bg="white")
  
  message("  Saved domain burden and plots for ", ph)
}

message("\nDone.")
