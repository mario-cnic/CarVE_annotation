# Unit tests for src/R/statsJPO.R helper functions

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

source("src/R/statsJPO.R")

cat("Running R unit tests for statsJPO.R...\n")

# Test 1: Data frame check in contingency_analysis
test_df <- data.frame(
  group = factor(c("A", "A", "B", "B", "A", "B")),
  status = factor(c("Case", "Control", "Case", "Control", "Case", "Control"))
)

res_contingency <- tryCatch({
  contingency_analysis(
    dataframe = test_df,
    variable = "status",
    group_var = "group",
    category_filter = "Case",
    palette = c("red", "blue"),
    plot_title = "Test Contingency",
    pairwise_comparisons = FALSE
  )
}, error = function(e) {
  cat("ERROR in contingency_analysis:", e$message, "\n")
  return(NULL)
})

if (!is.null(res_contingency) && !is.null(res_contingency$chi2_test)) {
  cat("  [PASS] contingency_analysis completed successfully.\n")
} else {
  cat("  [FAIL] contingency_analysis failed.\n")
  quit(status = 1)
}

# Test 2: means_analysis
num_df <- data.frame(
  group = factor(c("A", "A", "A", "B", "B", "B")),
  val = c(10.5, 12.1, 11.8, 25.4, 23.9, 26.1)
)

res_means <- tryCatch({
  means_analysis(
    dataframe = num_df,
    variable = "val",
    group_var = "group",
    palette = c("A" = "red", "B" = "blue"),
    plot_title = "Test Means"
  )
}, error = function(e) {
  cat("ERROR in means_analysis:", e$message, "\n")
  return(NULL)
})

if (!is.null(res_means) && !is.null(res_means$summary_by_group)) {
  cat("  [PASS] means_analysis completed successfully.\n")
} else {
  cat("  [FAIL] means_analysis failed.\n")
  quit(status = 1)
}

cat("All statsJPO.R unit tests passed cleanly!\n")
