# 08_describe_ttn_results.R
# Basic descriptive summary for TTN parsed.clean.pq tables.
#
# Usage:
#   Rscript src/R/08_describe_ttn_results.R [TTN_RESULTS_DIR] [OUT_DIR]
#
# Defaults:
#   TTN_RESULTS_DIR = results/TTN
#   OUT_DIR         = TTN_RESULTS_DIR

args <- commandArgs(trailingOnly = TRUE)

ttn_dir <- if (length(args) >= 1) args[1] else file.path("results", "TTN")
out_dir <- if (length(args) >= 2) args[2] else ttn_dir

if (!dir.exists(ttn_dir)) {
  stop("TTN directory not found: ", ttn_dir)
}

if (!requireNamespace("arrow", quietly = TRUE)) {
  stop("Package 'arrow' is required for parquet input. Install with install.packages('arrow').")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pq_files <- list.files(
  ttn_dir,
  pattern = "\\.parsed\\.clean\\.pq$",
  full.names = TRUE
)

if (length(pq_files) == 0) {
  stop("No .parsed.clean.pq files found in: ", ttn_dir)
}

read_header <- function(path) {
  tab <- arrow::read_parquet(path, as_data_frame = FALSE)
  tab$schema$names
}

read_selected_columns <- function(path, selected_cols, header_cols) {
  if (length(selected_cols) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  df <- arrow::read_parquet(
    path,
    col_select = tidyselect::all_of(selected_cols),
    as_data_frame = TRUE
  )

  # Keep character-like handling consistent with the original TSV script.
  for (nm in names(df)) {
    if (!is.character(df[[nm]])) {
      df[[nm]] <- as.character(df[[nm]])
    }
  }
  df
}

clean_values <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | x == "" | x == "." | x == "-"] <- NA_character_
  x
}

variant_key <- function(df) {
  has_id <- "ID_MUTACION" %in% names(df)
  if (has_id) {
    ids <- clean_values(df$ID_MUTACION)
    if (sum(!is.na(ids)) > 0) {
      return(ids)
    }
  }

  needed <- c("CHROM", "POS", "REF", "ALT")
  if (!all(needed %in% names(df))) {
    return(rep(NA_character_, nrow(df)))
  }

  chrom <- clean_values(df$CHROM)
  pos <- clean_values(df$POS)
  ref <- clean_values(df$REF)
  alt <- clean_values(df$ALT)
  ifelse(
    is.na(chrom) | is.na(pos) | is.na(ref) | is.na(alt),
    NA_character_,
    paste(chrom, pos, ref, alt, sep = ":")
  )
}

summarize_file <- function(path) {
  header_cols <- read_header(path)

  wanted <- c(
    "CHROM", "POS", "REF", "ALT", "ID_MUTACION",
    "Feature", "Consequence", "IMPACT", "CANONICAL",
    "HGVSp", "Protein_position", "REVEL", "gnomADv4"
  )
  selected <- intersect(wanted, header_cols)

  df <- read_selected_columns(path, selected, header_cols)

  file_label <- sub("\\.parsed\\.clean\\.pq$", "", basename(path))
  vkey <- variant_key(df)

  consequence <- if ("Consequence" %in% names(df)) clean_values(df$Consequence) else rep(NA_character_, nrow(df))
  impact <- if ("IMPACT" %in% names(df)) clean_values(df$IMPACT) else rep(NA_character_, nrow(df))
  feature <- if ("Feature" %in% names(df)) clean_values(df$Feature) else rep(NA_character_, nrow(df))
  canonical <- if ("CANONICAL" %in% names(df)) clean_values(df$CANONICAL) else rep(NA_character_, nrow(df))
  hgvsp <- if ("HGVSp" %in% names(df)) clean_values(df$HGVSp) else rep(NA_character_, nrow(df))
  protein_pos <- if ("Protein_position" %in% names(df)) clean_values(df$Protein_position) else rep(NA_character_, nrow(df))

  revel_num <- rep(NA_real_, nrow(df))
  if ("REVEL" %in% names(df)) {
    revel_num <- suppressWarnings(as.numeric(clean_values(df$REVEL)))
  }

  gnomad_num <- rep(NA_real_, nrow(df))
  if ("gnomADv4" %in% names(df)) {
    gnomad_num <- suppressWarnings(as.numeric(clean_values(df$gnomADv4)))
  }

  data.frame(
    file = file_label,
    n_rows = nrow(df),
    n_unique_variants = length(unique(vkey[!is.na(vkey)])),
    n_transcripts = length(unique(feature[!is.na(feature)])),
    n_unique_hgvsp = length(unique(hgvsp[!is.na(hgvsp)])),
    n_unique_protein_positions = length(unique(protein_pos[!is.na(protein_pos)])),
    n_consequences = length(unique(consequence[!is.na(consequence)])),
    n_impact_levels = length(unique(impact[!is.na(impact)])),
    n_canonical_yes = sum(canonical == "YES", na.rm = TRUE),
    n_missense = sum(grepl("missense_variant", consequence, fixed = TRUE), na.rm = TRUE),
    n_impact_high = sum(impact == "HIGH", na.rm = TRUE),
    n_impact_moderate = sum(impact == "MODERATE", na.rm = TRUE),
    n_revel_available = sum(!is.na(revel_num)),
    n_revel_ge_0_50 = sum(revel_num >= 0.50, na.rm = TRUE),
    n_revel_ge_0_70 = sum(revel_num >= 0.70, na.rm = TRUE),
    n_gnomad_available = sum(!is.na(gnomad_num)),
    stringsAsFactors = FALSE
  )
}

extract_transcript_counts <- function(path) {
  header_cols <- read_header(path)
  file_label <- sub("\\.parsed\\.clean\\.pq$", "", basename(path))

  if (!("Feature" %in% header_cols)) {
    return(data.frame(
      file = character(0),
      transcript = character(0),
      n_variants = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  df <- read_selected_columns(path, "Feature", header_cols)
  transcript <- clean_values(df$Feature)
  transcript <- transcript[!is.na(transcript)]

  if (length(transcript) == 0) {
    return(data.frame(
      file = character(0),
      transcript = character(0),
      n_variants = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  counts <- sort(table(transcript), decreasing = TRUE)
  data.frame(
    file = file_label,
    transcript = names(counts),
    n_variants = as.integer(counts),
    stringsAsFactors = FALSE
  )
}

message("Summarizing TTN files in: ", ttn_dir)

summary_list <- lapply(pq_files, summarize_file)
summary_df <- do.call(rbind, summary_list)

transcript_list <- lapply(pq_files, extract_transcript_counts)
transcript_df <- do.call(rbind, transcript_list)

summary_out <- file.path(out_dir, "TTN_descriptive_summary_parquet.csv")
transcript_out <- file.path(out_dir, "TTN_transcript_counts_parquet.csv")

utils::write.csv(summary_df, summary_out, row.names = FALSE)
utils::write.csv(transcript_df, transcript_out, row.names = FALSE)

message("Saved summary: ", summary_out)
message("Saved transcript counts: ", transcript_out)

print(summary_df)