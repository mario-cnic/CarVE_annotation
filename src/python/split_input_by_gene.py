#!/usr/bin/env python3
"""
Multi-Gene Input Splitter for Annotation Pipeline
Splits a multi-gene input file (.parquet, .vcf, .tsv, .csv, .xlsx) into gene-specific files.
"""

import os
import sys
import argparse
import logging
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("split_input_by_gene")

def parse_args():
    parser = argparse.ArgumentParser(description="Split multi-gene input file by gene")
    parser.add_argument("--input", required=True, help="Input multi-gene file (.pq, .parquet, .tsv, .csv, .xlsx, .vcf)")
    parser.add_argument("--output-dir", required=True, help="Directory to save gene-specific input files")
    parser.add_argument("--gene-col", default=None, help="Column name containing gene symbols (auto-detected if None)")
    return parser.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    logger.info(f"Reading multi-gene input file: {args.input}...")
    ext = args.input.lower().split('.')[-1]

    if ext in ('pq', 'parquet'):
        df = pd.read_parquet(args.input)
    elif ext == 'tsv':
        df = pd.read_csv(args.input, sep='\t')
    elif ext in ('csv', 'txt'):
        df = pd.read_csv(args.input)
    elif ext == 'xlsx':
        df = pd.read_excel(args.input)
    elif args.input.endswith('.vcf') or args.input.endswith('.vcf.gz'):
        # Parse VCF into DataFrame with INFO/GENE parsing
        records = []
        import gzip
        open_fn = gzip.open if args.input.endswith('.gz') else open
        with open_fn(args.input, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    chrom, pos, id_val, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                    info = parts[7] if len(parts) > 7 else ""
                    gene = "UNKNOWN"
                    for item in info.split(';'):
                        if item.startswith('GENE=') or item.startswith('Gene='):
                            gene = item.split('=')[1]
                            break
                    records.append({'CHROM': chrom, 'POS': pos, 'ID': id_val, 'REF': ref, 'ALT': alt, 'GENE': gene, 'INFO': info})
        df = pd.DataFrame(records)
    else:
        raise ValueError(f"Unsupported file format: {args.input}")

    # Auto-detect gene column
    gene_col = args.gene_col
    if not gene_col:
        for c in df.columns:
            if c.lower() in ('gene', 'gen', 'symbol', 'hgnc_symbol', 'gene_name'):
                gene_col = c
                break

    if not gene_col or gene_col not in df.columns:
        logger.error("Could not auto-detect gene column. Please specify with --gene-col.")
        sys.exit(1)

    unique_genes = df[gene_col].dropna().astype(str).unique()
    logger.info(f"Found {len(unique_genes)} unique genes: {list(unique_genes)}")

    for gene, gdf in df.groupby(gene_col):
        gene_clean = str(gene).strip().upper().replace('/', '_')
        out_file = os.path.join(args.output_dir, f"{gene_clean}.pq")
        gdf.to_parquet(out_file, index=False)
        logger.info(f"Saved {len(gdf)} variants for gene {gene_clean} to: {out_file}")

if __name__ == "__main__":
    main()
