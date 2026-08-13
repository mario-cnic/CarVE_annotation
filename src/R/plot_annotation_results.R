#!/usr/bin/env Rscript

# Plot Annotation Results Script
# Generates publication-grade PDF and PNG plots from parsed variant dataset

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- NULL
output_dir <- "plots"

i <- 1
while (i <= length(args)) {
  if (args[i] == "-i" || args[i] == "--input") {
    input_file <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "-o" || args[i] == "--output-dir") {
    output_dir <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(input_file)) {
  stop("Input file must be supplied via --input <file>", call. = FALSE)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading input data from ", input_file, "...")
if (grepl("\\.pq$|\\.parquet$", input_file)) {
  if (requireNamespace("arrow", quietly = TRUE)) {
    df <- arrow::read_parquet(input_file)
  } else {
    temp_tsv <- tempfile(fileext = ".tsv")
    system(paste("python3 -c \"import pandas as pd; pd.read_parquet('", input_file, "').to_csv('", temp_tsv, "', sep='\\t', index=False)\"", sep = ""))
    df <- read_tsv(temp_tsv, show_col_types = FALSE)
    unlink(temp_tsv)
  }
} else {
  df <- read_tsv(input_file, show_col_types = FALSE)
}

# Color Palette Definitions
theme_publication <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 13),
    legend.title = element_text(face = "bold")
  )

# 1. Variant Consequence Distribution Plot
if ("Consequence" %in% colnames(df)) {
  message("Generating Consequence Distribution Plot...")
  cons_df <- df %>%
    count(Consequence) %>%
    arrange(desc(n)) %>%
    head(15)

  p1 <- ggplot(cons_df, aes(x = reorder(Consequence, n), y = n, fill = n)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    scale_fill_viridis_c(option = "magma") +
    labs(title = "Variant Consequence Distribution", x = "VEP Consequence", y = "Variant Count") +
    theme_publication

  ggsave(file.path(output_dir, "consequence_distribution.png"), p1, width = 9, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "consequence_distribution.pdf"), p1, width = 9, height = 6)
}

# 2. Pathogenicity Predictor Correlation (REVEL vs AlphaMissense)
revel_col <- intersect(c("REVEL_score", "REVEL", "revel"), colnames(df))[1]
am_col <- intersect(c("am_pathogenicity", "AlphaMissense", "alphamissense"), colnames(df))[1]

if (!is.na(revel_col) && !is.na(am_col)) {
  message("Generating Pathogenicity Predictor Correlation Plot...")
  patho_df <- df %>%
    mutate(
      REVEL = as.numeric(.data[[revel_col]]),
      AlphaMissense = as.numeric(.data[[am_col]])
    ) %>%
    filter(!is.na(REVEL) & !is.na(AlphaMissense))

  if (nrow(patho_df) > 0) {
    p2 <- ggplot(patho_df, aes(x = REVEL, y = AlphaMissense)) +
      geom_point(alpha = 0.6, color = "#2c3e50") +
      geom_smooth(method = "lm", color = "#e74c3c", se = TRUE) +
      geom_vline(xintercept = 0.75, linetype = "dashed", color = "#7f8c8d") +
      geom_hline(yintercept = 0.80, linetype = "dashed", color = "#7f8c8d") +
      labs(title = "REVEL vs AlphaMissense Pathogenicity Scores", x = "REVEL Score", y = "AlphaMissense Score") +
      theme_publication

    ggsave(file.path(output_dir, "pathogenicity_correlation.png"), p2, width = 8, height = 6, dpi = 300)
    ggsave(file.path(output_dir, "pathogenicity_correlation.pdf"), p2, width = 8, height = 6)
  }
}

# 3. Allele Frequency Spectrum
af_col <- intersect(c("AF", "gnomAD_AF", "gnomAD_AF_joint", "AF_joint"), colnames(df))[1]
if (!is.na(af_col)) {
  message("Generating Allele Frequency Spectrum Plot...")
  af_df <- df %>%
    mutate(AF_num = as.numeric(.data[[af_col]])) %>%
    filter(!is.na(AF_num) & AF_num > 0)

  if (nrow(af_df) > 0) {
    p3 <- ggplot(af_df, aes(x = AF_num)) +
      geom_histogram(bins = 40, fill = "#3498db", color = "#2980b9", alpha = 0.8) +
      scale_x_log10() +
      labs(title = "gnomAD Allele Frequency Spectrum (Log Scale)", x = "gnomAD Allele Frequency (AF)", y = "Variant Count") +
      theme_publication

    ggsave(file.path(output_dir, "allele_frequency_spectrum.png"), p3, width = 8, height = 5, dpi = 300)
    ggsave(file.path(output_dir, "allele_frequency_spectrum.pdf"), p3, width = 8, height = 5)
  }
}

message("Static publication plots successfully generated in: ", output_dir)
