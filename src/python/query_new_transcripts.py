#!/usr/bin/env python3
"""
query_new_transcripts.py

Finds and resolves RefSeq (NM) and Ensembl (ENST) transcript mappings for new genes,
supporting both raw Parquet folder audits and single gene symbols. Uses local reference
databases first for speed and cluster safety, falling back to Ensembl REST API.

Usage:
    python3 src/python/query_new_transcripts.py <INPUT> [--mapping-file PATH] [--auto-append]
    
    # Examples:
    python3 src/python/query_new_transcripts.py _RAW/data_Valencia_plus_hic/oracle
    python3 src/python/query_new_transcripts.py ACADVL --auto-append
"""

import os
import glob
import re
import sys
import time
import argparse
import urllib.request
import json
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Query transcript mappings for a single gene symbol or a folder/file of raw cases files."
    )
    parser.add_argument(
        "input",
        type=str,
        help="Input gene symbol (e.g. ACADVL) or path to raw cases folder/file"
    )
    parser.add_argument(
        "--mapping-file",
        type=str,
        default="resources/gene_transcript_mapping.txt",
        help="Path to gene-to-transcript mapping file"
    )
    parser.add_argument(
        "--auto-append",
        action="store_true",
        help="Automatically append resolved missing transcripts to the mapping file"
    )
    return parser.parse_args()

def get_local_db_path():
    """Returns local ensembl_to_refseq.tsv.gz path if found in resources, local workspace, or shared mounts."""
    paths = [
        "resources/ensembl_to_refseq.tsv.gz",
        "ensembl_to_refseq.tsv.gz",
        "shared/utils/data/ensembl_to_refseq.tsv.gz",
        "/home/mruizp/data_lab_PGP/shared/utils/data/ensembl_to_refseq.tsv.gz",
        "/data_lab_PGP/shared/utils/data/ensembl_to_refseq.tsv.gz"
    ]
    for p in paths:
        if os.path.exists(p):
            return p
    return None

def extract_refseq_from_raw(file_path):
    """Loads raw dataset (parquet, excel, or csv) and extracts the first non-null RefSeq NM accession found in any column."""
    if not os.path.exists(file_path):
        return None
    try:
        ext = os.path.splitext(file_path)[1].lower()
        if ext in ['.parquet', '.pq']:
            df = pd.read_parquet(file_path)
        elif ext == '.xlsx':
            df = pd.read_excel(file_path)
        elif ext == '.csv':
            df = pd.read_csv(file_path)
        else:
            return None
        
        # Search every column for NM_ accession
        for col in df.columns:
            col_series = df[col].dropna().astype(str)
            for val in col_series:
                match = re.search(r'(NM_\d+(?:\.\d+)?)', val)
                if match:
                    return match.group(1)
    except Exception as e:
        print(f"Warning: Failed reading {file_path}: {e}")
    return None

def fetch_canonical_enst_api(gene_symbol):
    """Queries Ensembl REST API for canonical transcript Ensembl ID."""
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=5) as r:
                data = json.loads(r.read().decode('utf-8'))
                
            for t in data.get('Transcript', []):
                if t.get('is_canonical') == 1:
                    return t.get('id')
            if data.get('Transcript'):
                return data['Transcript'][0].get('id')
        except Exception as e:
            time.sleep(1)
    return None

def fetch_refseq_xref_api(enst_id):
    """Queries Ensembl REST API for RefSeq mRNA xrefs of an ENST ID."""
    if not enst_id or enst_id == "N/A":
        return None
    for db in ['RefSeq_mRNA', 'RefSeq_mRNA_predicted']:
        xref_url = f"https://rest.ensembl.org/xrefs/id/{enst_id}?external_db={db}"
        xref_req = urllib.request.Request(xref_url, headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(xref_req, timeout=5) as xr:
                xref_data = json.loads(xr.read().decode('utf-8'))
            if xref_data:
                refseq_nm = xref_data[0].get('display_id') or xref_data[0].get('primary_id')
                if refseq_nm:
                    return refseq_nm
        except Exception:
            continue
    return None

def resolve_gene(gene, raw_file_path=None, df_local=None):
    """Resolves RefSeq NM and Ensembl ENST canonical transcripts for a gene."""
    refseq = None
    enst = None
    
    # 1. If we have a raw file, extract the exact RefSeq NM used in the dataset
    if raw_file_path:
        refseq = extract_refseq_from_raw(raw_file_path)
        
    # 2. Try to resolve from local database
    if df_local is not None:
        matches = df_local[(df_local['Gene name'].str.upper() == gene.upper()) & (df_local['Ensembl Canonical'] == 1.0)]
        if not matches.empty:
            row = matches.iloc[0]
            enst = row.get('Transcript stable ID')
            if not refseq:
                refseq = row.get('RefSeq match transcript (MANE Select)')
            print(f"  {gene} resolved locally (MANE Select)")
            
    # 3. Fallback to API queries if missing information
    if not enst or enst == "N/A" or not refseq or refseq == "N/A":
        print(f"  {gene} info incomplete. Querying Ensembl REST API...")
        if not enst:
            enst = fetch_canonical_enst_api(gene)
        if not refseq:
            refseq = fetch_refseq_xref_api(enst)
            
    # Normalize defaults
    enst = enst if enst else "N/A"
    refseq = refseq if refseq else "N/A"
    
    return refseq, enst

def main():
    args = parse_args()
    
    # Check mapping file
    if not os.path.exists(args.mapping_file):
        print(f"Error: Mapping file {args.mapping_file} does not exist.")
        sys.exit(1)
        
    df_map = pd.read_csv(args.mapping_file)
    mapped_genes = set(df_map['Gen'].dropna().str.strip().str.upper())
    
    # Determine mode: RAW folder/file vs Gene Symbol
    is_path = os.path.exists(args.input)
    
    targets = [] # list of tuples: (gene_symbol, file_path_or_None)
    
    if is_path:
        # Path mode: scan files
        if os.path.isdir(args.input):
            extensions = ['*.parquet', '*.pq', '*.xlsx', '*.csv', '*.PARQUET', '*.PQ', '*.XLSX', '*.CSV']
            files = []
            for ext in extensions:
                files.extend(glob.glob(os.path.join(args.input, ext)))
            files = sorted(list(set(files)))
            for f in files:
                base_name = os.path.splitext(os.path.basename(f))[0]
                gene = base_name.split('_')[0].strip().upper()
                if gene not in mapped_genes:
                    targets.append((gene, f))
        else:
            base_name = os.path.splitext(os.path.basename(args.input))[0]
            gene = base_name.split('_')[0].strip().upper()
            if gene not in mapped_genes:
                targets.append((gene, args.input))
                
        if not targets:
            print(f"All genes in '{args.input}' are already mapped in {args.mapping_file}. Nothing to resolve.")
            sys.exit(0)
    else:
        # Single gene symbol mode
        gene = args.input.strip().upper()
        if gene in mapped_genes:
            print(f"Warning: Gene {gene} is already mapped in {args.mapping_file}.")
        targets.append((gene, None))
        
    # Load local database if available
    df_local = None
    local_db_path = get_local_db_path()
    if local_db_path:
        print(f"Loading local reference database from: {local_db_path} ...")
        try:
            df_local = pd.read_csv(local_db_path, sep='\t')
        except Exception as e:
            print(f"Warning: Failed to load local database: {e}. Falling back to REST API.")
    else:
        print("No local reference database found. Proceeding via Ensembl REST API.")
        
    resolved_lines = []
    print(f"\nProcessing {len(targets)} target(s)...")
    for gene, path in targets:
        refseq, enst = resolve_gene(gene, path, df_local)
        resolved_lines.append(f"{gene},{refseq},{enst}")
        print(f"  -> {gene}: RefSeq={refseq} | Ensembl={enst}")
        time.sleep(0.1)
        
    # Output formatting
    print("\n==================================================")
    print("Resolved Transcript Lines:")
    print("==================================================")
    for line in resolved_lines:
        print(line)
    print("==================================================\n")
    
    # Auto-append if requested
    if args.auto_append:
        lines_to_append = []
        for line in resolved_lines:
            parts = line.split(',')
            g = parts[0].strip().upper()
            refseq = parts[1].strip()
            enst = parts[2].strip()
            # Only append if the gene is not already mapped and we successfully resolved it (not N/A)
            if g not in mapped_genes and refseq != "N/A" and enst != "N/A":
                lines_to_append.append(line)
        
        if lines_to_append:
            try:
                # Read the file contents to see if it ends with a newline
                with open(args.mapping_file, "r") as f_in:
                    content = f_in.read()
                needs_newline = len(content) > 0 and not content.endswith("\n")
                
                with open(args.mapping_file, "a") as f_out:
                    if needs_newline:
                        f_out.write("\n")
                    for line in lines_to_append:
                        f_out.write(line + "\n")
                print(f"Successfully appended {len(lines_to_append)} transcript mappings to {args.mapping_file}")
            except Exception as e:
                print(f"Error appending to mapping file: {e}")
        else:
            print("No new mappings to append (all targets already mapped).")

if __name__ == "__main__":
    main()
