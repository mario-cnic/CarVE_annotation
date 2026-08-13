import pandas as pd
import mygene
import warnings

def map_refseq_to_ensembl(input_file, output_file):
    # 1. Read the input file (assuming tab-separated, no headers)
    print(f"Reading data from {input_file}...")
    try:
        df = pd.read_csv(input_file, sep="\t", header=None, names=["gene_name", "original_refseq", "count"])
    except FileNotFoundError:
        print(f"Error: Could not find the file {input_file}")
        return

    # 2. Strip the version numbers (e.g., changing NM_005159.4 to NM_005159)
    # The regex \.\d+$ looks for a dot followed by numbers at the very end of the string
    df['query_refseq'] = df['original_refseq'].str.replace(r'\.\d+$', '', regex=True)

    # Extract unique IDs to minimize API calls
    unique_refseqs = df['query_refseq'].dropna().unique().tolist()

    # 3. Query the MyGene.info API
    print(f"Querying MyGene.info for {len(unique_refseqs)} unique RefSeq IDs...")
    mg = mygene.MyGeneInfo()
    
    # querymany handles batching automatically and is very fast
    results = mg.querymany(
        unique_refseqs,
        scopes='refseq',
        fields='ensembl.transcript',
        species='human',
        as_dataframe=True,
        returnall=False
    )

    # 4. Clean up the results
    # Reset index so 'query' (our stripped RefSeq ID) becomes a standard column
    results = results.reset_index()
    
    # MyGene's dataframe output can sometimes nest results if a gene has multiple Ensembl matches.
    # We need a helper function to safely extract the first transcript string.
    # We need a helper function to safely extract the first transcript string.
    def extract_transcript(val):
        # 1. If it's a list, grab the first item immediately (bypassing pd.isna)
        if isinstance(val, list):
            return str(val[0]) if len(val) > 0 else None
            
        # 2. If it's a float (which is what pandas uses for NaN) or None, it's missing
        if isinstance(val, float) or val is None:
            return None
            
        # 3. Otherwise, return the value (should be a string)
        return str(val)
    # Depending on the query results, the column might be named 'ensembl.transcript'
    if 'ensembl.transcript' in results.columns:
        results['ensembl_transcript_id'] = results['ensembl.transcript'].apply(extract_transcript)
    else:
        results['ensembl_transcript_id'] = None

    # 5. Merge the results back with the original dataframe
    # We do a left merge to ensure we keep every row from the original input file
    final_df = pd.merge(
        df, 
        results[['query', 'ensembl_transcript_id']], 
        left_on='query_refseq', 
        right_on='query', 
        how='left'
    )

    # Keep only the requested columns
    final_df = final_df[['gene_name', 'original_refseq', 'ensembl_transcript_id', 'count']]
    # Check for missing mappings and issue a warning
    missing_count = final_df['ensembl_transcript_id'].isna().sum()
    final_df.rename(columns={'gene_name': 'GEN', 'original_refseq': 'NM', 'ensembl_transcript_id': 'ENST'}, inplace=True)
    final_df.sort_values(['GEN','count','NM'], ascending=[True, False, True], inplace=True)
    if missing_count > 0:
        warnings.warn(f"\nWarning: {missing_count} RefSeq IDs could not be mapped. They will appear as empty/NaN in the output.")

    # 6. Save to output file
    final_df.to_csv(output_file, sep=",", index=False)
    print(f"Success! Mapping saved to {output_file}")
    
    # Print a preview
    print("\nPreview of mapped data:")
    print(final_df.head())

# --- Run the script ---
if __name__ == "__main__":
    INPUT_FILE = "_RAW/data_HIC/data_hic_transcripts.txt"
    OUTPUT_FILE = "gene_transcript_mapping_datahic.txt"
    
    map_refseq_to_ensembl(INPUT_FILE, OUTPUT_FILE)