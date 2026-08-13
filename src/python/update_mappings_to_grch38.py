#!/usr/bin/env python3
import os
import time
import urllib.request
import json
import pandas as pd

# Load mapping file
mapping_path = "resources/gene_transcript_mapping.txt"
df_map = pd.read_csv(mapping_path)

# Load local database
local_db_path = "/home/mruizp/data_lab_PGP/shared/utils/data/ensembl_to_refseq.tsv.gz"
if os.path.exists(local_db_path):
    print(f"Loading local reference database from: {local_db_path}")
    df_local = pd.read_csv(local_db_path, sep='\t')
else:
    df_local = None
    print("No local reference database found. Proceeding via Ensembl REST API.")

def fetch_canonical_enst_api(gene_symbol):
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=5) as r:
            data = json.loads(r.read().decode('utf-8'))
        for t in data.get('Transcript', []):
            if t.get('is_canonical') == 1:
                return t.get('id')
        if data.get('Transcript'):
            return data['Transcript'][0].get('id')
    except Exception as e:
        print(f"API Error for {gene_symbol}: {e}")
    return None

updated_count = 0
for idx, row in df_map.iterrows():
    gene = str(row['Gen']).strip()
    curr_enst = str(row['ENST']).strip()
    
    # Try resolving canonical GRCh38 ENST ID from local db first
    resolved_enst = None
    if df_local is not None:
        matches = df_local[(df_local['Gene name'].str.upper() == gene.upper()) & (df_local['Ensembl Canonical'] == 1.0)]
        if not matches.empty:
            # Prioritize rows that have a RefSeq MANE Select mapping
            matches = matches.sort_values(by='RefSeq match transcript (MANE Select)', na_position='last')
            row = matches.iloc[0]
            resolved_enst = row.get('Transcript stable ID')
            
    # Fallback to API if not resolved locally
    if not resolved_enst:
        resolved_enst = fetch_canonical_enst_api(gene)
        time.sleep(0.1) # Rate limit safety
        
    if resolved_enst and resolved_enst != curr_enst:
        print(f"Updating {gene}: {curr_enst} -> {resolved_enst}")
        df_map.at[idx, 'ENST'] = resolved_enst
        updated_count += 1

if updated_count > 0:
    # Save the updated mapping back to the file
    df_map.to_csv(mapping_path, index=False)
    print(f"\nSuccessfully updated {updated_count} transcript mappings to GRCh38 in {mapping_path}")
else:
    print("\nAll transcripts are already up-to-date with GRCh38 canonical transcripts.")
