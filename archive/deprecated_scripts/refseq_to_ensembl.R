# Install biomaRt if you haven't already:
# BiocManager::install("biomaRt")

library(biomaRt)

# 1. Read the data efficiently into a data frame
mapping_file <- "_RAW/data_HIC/data_hic_transcripts.txt"
# Assuming tab-separated with no header based on your original code
input_data <- read.table(mapping_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(input_data) <- c("gene_name", "original_refseq")

# 2. Strip the version numbers for the BioMart query
# This regex removes the dot and any numbers at the very end of the string
input_data$query_refseq <- gsub("\\.[0-9]+$", "", input_data$original_refseq)

# 3. Connect to the human Ensembl dataset
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 4. Query BioMart ONCE for all unique IDs
# This prevents the MySQL server from dropping your connection
message("Querying BioMart... Please wait.")
mapping <- getBM(
  attributes = c('refseq_mrna', 'ensembl_transcript_id'),
  filters = 'refseq_mrna',
  values = unique(input_data$query_refseq),
  mart = ensembl
)

# 5. Merge the results back with your original data
# all.x = TRUE ensures that even if a mapping isn't found, the original row is kept (with NA)
final_mapping <- merge(
  x = input_data,
  y = mapping,
  by.x = "query_refseq",
  by.y = "refseq_mrna",
  all.x = TRUE 
)

# 6. Clean up and reorder columns (Gene Name, Original RefSeq, Ensembl ID)
final_mapping <- final_mapping[, c("gene_name", "original_refseq", "ensembl_transcript_id")]

# Optional: Print a warning for any IDs that STILL didn't map
missing_mappings <- final_mapping$original_refseq[is.na(final_mapping$ensembl_transcript_id)]
if(length(missing_mappings) > 0) {
  warning(paste(length(missing_mappings), "RefSeq IDs still could not be mapped. They will have 'NA' in the output."))
}

# 7. Save the mapping to a file
output_file <- "_RAW/data_HIC/data_hic_transcript_mapping.txt"
write.table(final_mapping,
            file = output_file,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

message("Mapping complete! Saved to: ", output_file)