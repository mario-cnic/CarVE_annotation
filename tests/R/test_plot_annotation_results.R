# Unit test for src/R/plot_annotation_results.R

cat("Running R unit test for plot_annotation_results.R...\n")

test_dir <- tempdir()
tsv_path <- file.path(test_dir, "sample.tsv")
plots_dir <- file.path(test_dir, "test_plots")

# Create mini input dataset
df <- data.frame(
  Consequence = c("missense_variant", "splice_donor_variant", "stop_gained", "missense_variant"),
  REVEL = c(0.85, 0.50, NA, 0.92),
  AlphaMissense = c(0.90, 0.30, 0.98, 0.95),
  gnomAD_AF_joint = c(0.0001, 0.00005, 0.001, 0.00001)
)

write.table(df, file = tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)

cmd <- paste("Rscript src/R/plot_annotation_results.R --input", tsv_path, "--output-dir", plots_dir)
status <- system(cmd)

if (status != 0) {
  cat("  [FAIL] Rscript execution failed with status:", status, "\n")
  quit(status = 1)
}

png1 <- file.path(plots_dir, "consequence_distribution.png")
png2 <- file.path(plots_dir, "pathogenicity_correlation.png")
png3 <- file.path(plots_dir, "allele_frequency_spectrum.png")

if (file.exists(png1) && file.exists(png2) && file.exists(png3)) {
  cat("  [PASS] All static PNG/PDF figures successfully generated.\n")
} else {
  cat("  [FAIL] Missing generated plot files in output directory.\n")
  quit(status = 1)
}
