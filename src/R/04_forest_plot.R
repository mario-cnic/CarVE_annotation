# 04_forest_plot.R
# Forest plots: one per gene × phenotype.
# Y axis: 4 variant classes (facets) × 3 comparators (int / gnomAD NFE / gnomAD total).
# X axis: OR on log scale with 95% CI.
#
# Usage:   Rscript scripts/04_forest_plot.R <GENE>
# Example: Rscript scripts/04_forest_plot.R KLHL24
#
# Outputs (figures/):
#   forest_{gene}_{phenotype}.png   — one per phenotype (7 per gene)
#   forest_{gene}_all.pdf           — combined, all phenotypes

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript scripts/04_forest_plot.R <GENE>")
gene <- toupper(args[1])
raw_root <- args[2]
burden_path <- file.path(raw_root, "results", gene, paste0("burden_", tolower(gene), ".rds"))
if (!file.exists(burden_path)) stop("Burden not found: ", burden_path)

out_dir_fig <- file.path(raw_root, "figures", gene)
dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

burden <- readRDS(burden_path)

# ── Constants ─────────────────────────────────────────────────────────────────
# Dynamically extract phenotypes from the burden data
PHENOTYPES   <- sort(unique(burden$phenotype))

# Identify non-empty phenotypes for the PDF (at least one case AC > 0)
NON_EMPTY_PHENOTYPES <- burden |>
  group_by(phenotype) |>
  summarise(total_ac = sum(ac_case, na.rm = TRUE), .groups = "drop") |>
  filter(total_ac > 0) |>
  pull(phenotype)

CLASS_ORDER  <- c("truncating", "missense_high", "missense", "synonymous")
CLASS_LABELS <- c(
  truncating    = "Truncating",
  missense_high = "Missense\n(REVEL > 0.6)",
  missense      = "Missense\n(all)",
  synonymous    = "Synonymous"
)

COMP_ORDER  <- c("int", "nfe", "total")
COMP_LABELS <- c(int = "Internal", nfe = "gnomAD NFE", total = "gnomAD total")
COMP_COLORS <- c(int = "#E41A1C", nfe = "#377EB8", total = "#4DAF4A")

# ── Reshape burden to long format ─────────────────────────────────────────────
long <- burden |>
  select(
    gene, variant_class, phenotype,
    n_variants_cohort, n_variants_all, ac_case,
    or_int,    ci_low_int,    ci_high_int,    p_int,
    or_nfe,    ci_low_nfe,    ci_high_nfe,    p_nfe,
    or_total,  ci_low_total,  ci_high_total,  p_total
  ) |>
  pivot_longer(
    cols          = matches("^(or|ci_low|ci_high|p)_(int|nfe|total)$"),
    names_to      = c(".value", "comparator"),
    names_pattern = "^(or|ci_low|ci_high|p)_(int|nfe|total)$"
  ) |>
  mutate(
    variant_class = factor(variant_class, levels = CLASS_ORDER),
    comparator    = factor(comparator,    levels = COMP_ORDER,
                           labels = unname(COMP_LABELS))
  )

# Significance flag and p-value label (shown only when p < 0.05)
long <- long |>
  mutate(
    sig     = coalesce(!is.na(p) & p < 0.05, FALSE),
    p_label = case_when(
      is.na(p)  ~ "",
      p < 0.001 ~ formatC(p, format = "e", digits = 1),
      p < 0.05  ~ paste0("p=", sprintf("%.3f", p)),
      TRUE      ~ ""
    )
  )

# ── Plot factory ──────────────────────────────────────────────────────────────
make_forest_plot <- function(d, phenotype_name) {

  # Subtitle: case AC per variant class
  subtitle <- d |>
    distinct(variant_class, ac_case) |>
    arrange(variant_class) |>
    mutate(lbl = paste0(CLASS_LABELS[as.character(variant_class)], ": AC=", ac_case)) |>
    pull(lbl) |>
    paste(collapse = "   ·   ")

  # x-axis limits derived from CI range; always include OR=1
  ci_vals <- c(d$ci_low[!is.na(d$ci_low)], d$ci_high[!is.na(d$ci_high)])
  if (length(ci_vals) == 0) {
    x_min <- 0.1; x_max <- 10
  } else {
    x_min <- max(0.05, min(ci_vals) * 0.6)
    x_max <- min(500,  max(ci_vals) * 2.0)
  }
  # Always bracket OR=1
  x_min <- min(x_min, 0.8)
  x_max <- max(x_max, 1.2)

  ggplot(d, aes(y = comparator, x = or, color = comparator)) +
    facet_grid(
      variant_class ~ .,
      labeller = labeller(variant_class = CLASS_LABELS),
      scales   = "free_y",
      space    = "free_y"
    ) +
    # Reference line at OR = 1
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    # CI error bars (non-NA only)
    geom_errorbar(
      data      = filter(d, !is.na(or)),
      aes(xmin  = ci_low, xmax = ci_high),
      width     = 0.25, linewidth = 0.7,
      orientation = "y"
    ) +
    # Points: filled = significant, open = not
    geom_point(
      data = filter(d, !is.na(or)),
      aes(shape = sig),
      size = 3.5
    ) +
    # Label where ac_case = 0 (OR is NA by design)
    geom_text(
      data   = filter(d, is.na(or)),
      aes(label = "AC=0", x = 1),
      color  = "grey60", size = 2.8, hjust = 0.5
    ) +
    # p-value annotation for significant results
    geom_text(
      data  = filter(d, sig),
      aes(x = ci_high, label = paste0("  ", p_label)),
      hjust = 0, size = 2.6, show.legend = FALSE
    ) +
    scale_x_log10(
      limits = c(x_min, x_max),
      breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100),
      labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20", "50", "100"),
      oob    = scales::squish
    ) +
    scale_color_manual(
      values = setNames(COMP_COLORS, unname(COMP_LABELS)),
      guide  = guide_legend(override.aes = list(shape = 16, size = 3))
    ) +
    scale_shape_manual(
      values = c("FALSE" = 1, "TRUE" = 16),
      guide  = "none"
    ) +
    labs(
      title    = paste0(gene, "  -  ", phenotype_name),
      subtitle = paste0("Case AC:  ", subtitle),
      x        = "Odds Ratio  (log scale, 95% CI)",
      y        = NULL,
      color    = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background   = element_rect(fill = "grey92", color = NA),
      strip.text         = element_text(face = "bold", size = 9.5),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position    = "bottom",
      legend.key.size    = unit(0.9, "lines"),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 8, color = "grey40"),
      plot.margin        = margin(8, 14, 8, 8)
    )
}

# ── Generate plots ────────────────────────────────────────────────────────────
message("\n── Generating forest plots ──────────────────────────────────")

plots <- map(PHENOTYPES, function(ph) {
  d <- filter(long, phenotype == ph)
  make_forest_plot(d, ph)
})
names(plots) <- PHENOTYPES

# Individual PNGs
walk(PHENOTYPES, function(ph) {
  out <- file.path(out_dir_fig, paste0("forest_", tolower(gene), "_", ph, ".png"))
  ggsave(out, plots[[ph]], width = 7, height = 7.5, dpi = 150)
  message("  Saved: ", out)
})

# Combined PDF (all non-empty phenotypes, one per page)
pdf_path <- file.path(out_dir_fig, paste0("forest_", tolower(gene), "_all.pdf"))
pdf(pdf_path, width = 7, height = 7.5)
walk(NON_EMPTY_PHENOTYPES, ~ print(plots[[.x]]))
invisible(dev.off())
message("  Saved PDF: ", pdf_path, " (included non-empty phenotypes: ", length(NON_EMPTY_PHENOTYPES), ")")

message("\nDone: ", length(PHENOTYPES), " total individual plots for gene ", gene)
