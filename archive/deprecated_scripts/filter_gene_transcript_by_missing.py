import pandas as pd

# Input files
mapping_file = 'gene_transcript_mapping_datahic.txt'
missing_file = 'missing_transcript_datahic.txt'
output_file = 'gene_transcript_mapping_datahic_missing.txt'

# If not found, try local path
import os
if not os.path.exists(missing_file):
    missing_file = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'missing_transcript.txt'))
    print(f"Trying missing_file at {missing_file}")

# Read missing genes
with open(missing_file) as f:
    missing_genes = set(line.strip() for line in f if line.strip())

# Read mapping file
df = pd.read_csv(mapping_file)

# Filter rows where GEN is in missing_genes
filtered = df[df['GEN'].isin(missing_genes)]

# Save result
filtered.to_csv(output_file, index=False)
print(f'Saved filtered mapping to {output_file} ({len(filtered)} rows)')