#!/usr/bin/env python3
"""
noonan_string_audit.py

Exploratory script to scan final parsed clean Parquet files for string corruptions
introduced by the global substring replacement bug (e.g., 'Noo.', 'domi.t', 'malig.t').
Generates a markdown report summarizing findings by gene, column, and values.

Usage:
    python3 src/python/noonan_string_audit.py --input results/combined --output results/combined/noonan_analysis_report.md
"""

import os
import glob
import argparse
import pandas as pd
from datetime import datetime

def parse_arguments():
    parser = argparse.ArgumentParser(description="Audit parquet files for text corruption.")
    parser.add_argument(
        "--input",
        type=str,
        default="results/combined",
        help="Path to a parquet file or a directory containing parquet files to audit"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results/combined/noonan_analysis_report.md",
        help="Path to write the markdown audit report"
    )
    return parser.parse_args()

def scan_file_for_corruption(fpath):
    """
    Scans a single parquet file for typical corrupted substrings:
    - 'Noo.' (corrupted 'Noonan')
    - 'domi.t' (corrupted 'dominant')
    - 'malig.t' (corrupted 'malignant')
    """
    results = []
    try:
        df = pd.read_parquet(fpath)
    except Exception as e:
        print(f"Error reading {fpath}: {e}")
        return results

    corruptions = {
        "Noo.": "Noonan",
        "domi.t": "dominant",
        "malig.t": "malignant"
    }

    for col in df.columns:
        # We only check columns that are of object/string dtype
        if df[col].dtype == object or pd.api.types.is_string_dtype(df[col]):
            # Get series representation
            series = df[col].astype(str)
            for pattern, original in corruptions.items():
                mask = series.str.contains(pattern, case=False, na=False)
                if mask.any():
                    matched_rows = df[mask]
                    count = len(matched_rows)
                    unique_matches = matched_rows[col].dropna().unique()
                    
                    results.append({
                        "column": col,
                        "pattern": pattern,
                        "inferred_original": original,
                        "count": count,
                        "examples": [str(x) for x in unique_matches[:3]]
                    })
    return results

def main():
    args = parse_arguments()

    print(f"Auditing parquet data from: {args.input}")
    
    # Resolve target files
    if os.path.isdir(args.input):
        pq_files = sorted(glob.glob(os.path.join(args.input, "*.parquet")) + glob.glob(os.path.join(args.input, "*.pq")))
    else:
        pq_files = [args.input] if os.path.exists(args.input) else []

    if not pq_files:
        print(f"No parquet files found at: {args.input}")
        return

    print(f"Found {len(pq_files)} files to audit.")

    audit_summary = {}
    total_corruptions_found = 0

    for fpath in pq_files:
        gene = os.path.basename(fpath).split(".")[0]
        file_results = scan_file_for_corruption(fpath)
        if file_results:
            audit_summary[gene] = file_results
            total_corruptions_found += len(file_results)

    # Make output directory if it doesn't exist
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Write Markdown Report
    with open(args.output, "w") as f:
        f.write("# String Corruption Audit Report (Noonan Syndrome & 'nan' Bug)\n\n")
        f.write(f"**Date/Time of Audit**: {datetime.now().strftime('%Y-%m-%dT%H:%M:%S')}\n")
        f.write(f"**Target Directory**: `{args.input}`\n")
        f.write(f"**Files Scanned**: {len(pq_files)}\n\n")

        f.write("## 1. Executive Summary\n\n")
        if total_corruptions_found == 0:
            f.write("✅ **No string corruptions found.** All scanned files are clean and do not contain corrupted substrings like `'Noo.'`, `'domi.t'`, or `'malig.t'`.\n\n")
        else:
            f.write(f"⚠️ **Corrupted strings detected in {len(audit_summary)} gene files.**\n")
            f.write("The global substring replacement `.replace(\"nan\", \".\")` inside `vcf_parser_pysam.py` converted occurrences of `'nan'` inside actual words into `'.'`. ")
            f.write("This resulted in words like *Noonan*, *dominant*, and *malignant* being truncated/corrupted in the final parquet outputs.\n\n")

            f.write("## 2. Detailed Findings per Gene\n\n")
            for gene, results in sorted(audit_summary.items()):
                f.write(f"### Gene: **{gene}**\n\n")
                f.write("| Column Name | Detected Pattern | Inferred Original | Match Count | Example Values |\n")
                f.write("| :--- | :--- | :--- | :--- | :--- |\n")
                for res in results:
                    examples_str = "; ".join(f"`{ex}`" for ex in res["examples"])
                    f.write(f"| `{res['column']}` | `{res['pattern']}` | {res['inferred_original']} | {res['count']} | {examples_str} |\n")
                f.write("\n")

        f.write("## 3. Resolution Details\n\n")
        f.write("The bug has been fixed in the shared utility file `vcf_parser_pysam.py` by replacing the global substring search with a safe exact value comparison (`val_str.lower() == 'nan'`). ")
        f.write("This prevents any future pipeline runs from corrupting words containing `'nan'` while preserving the correct representation of missing values.\n")

    print(f"Audit completed. Markdown report written to: {args.output}")

if __name__ == "__main__":
    main()
