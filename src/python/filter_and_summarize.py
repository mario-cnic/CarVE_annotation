#!/usr/bin/env python3
"""
Customizable Variant Filtering and Summarization Engine

Filters parsed variant tables based on user-defined criteria:
  - gnomAD Allele Frequency (AF) thresholds
  - Pathogenicity scores (REVEL, AlphaMissense, CADD)
  - Splicing scores (SPiP, SpliceAI)
  - Consequence categories (High, Moderate, Splicing, etc.)

Outputs filtered dataset (.pq, .tsv, .xlsx) and detailed filtering metrics summary.
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["ARROW_IO_THREADS"] = "1"

import sys
import argparse
import logging
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("filter_and_summarize")

def parse_args():
    parser = argparse.ArgumentParser(description="Customizable Variant Filtering Engine")
    parser.add_argument("--input", required=True, help="Input .pq, .tsv file, or directory containing parsed files")
    parser.add_argument("--output-dir", required=True, help="Directory to save filtered results and summaries")
    parser.add_argument("--max-af", type=float, default=None, help="Maximum gnomAD allele frequency threshold (e.g. 0.001)")
    parser.add_argument("--min-revel", type=float, default=None, help="Minimum REVEL pathogenicity score threshold (e.g. 0.75)")
    parser.add_argument("--min-alphamissense", type=float, default=None, help="Minimum AlphaMissense score threshold (e.g. 0.80)")
    parser.add_argument("--min-spip", type=float, default=None, help="Minimum SPiP splicing score threshold (e.g. 0.50)")
    parser.add_argument("--min-cadd", type=float, default=None, help="Minimum CADD phred score threshold (e.g. 20.0)")
    parser.add_argument("--min-qual", type=float, default=None, help="Minimum VCF QUAL phred score threshold (e.g. 30.0)")
    parser.add_argument("--pass-qc-only", action="store_true", help="Filter out variants with LOW_QUAL status in QC_STATUS")
    parser.add_argument("--consequences", type=str, default=None, help="Comma-separated list of target VEP consequences")
    parser.add_argument("--output-format", choices=["pq", "parquet", "tsv", "xlsx"], default="pq", help="Output file format")
    parser.add_argument("--output-prefix", type=str, default=None, help="Prefix (e.g. gene name) for output filenames, to avoid collisions when multiple runs share --output-dir")
    return parser.parse_args()

def load_data(input_path):
    if os.path.isdir(input_path):
        files = [os.path.join(input_path, f) for f in os.listdir(input_path) if f.endswith('.pq') or f.endswith('.parquet') or f.endswith('.tsv')]
        if not files:
            raise FileNotFoundError(f"No .pq or .tsv files found in directory {input_path}")
        logger.info(f"Loading {len(files)} parsed files from directory {input_path}...")
        dfs = []
        for f in files:
            if f.endswith('.pq') or f.endswith('.parquet'):
                dfs.append(pd.read_parquet(f))
            else:
                dfs.append(pd.read_csv(f, sep='\t'))
        return pd.concat(dfs, ignore_index=True)
    elif input_path.endswith('.pq') or input_path.endswith('.parquet'):
        return pd.read_parquet(input_path)
    elif input_path.endswith('.tsv') or input_path.endswith('.csv'):
        sep = '\t' if input_path.endswith('.tsv') else ','
        return pd.read_csv(input_path, sep=sep)
    else:
        raise ValueError(f"Unsupported input file format: {input_path}")

def filter_dataframe(df, max_af=None, min_revel=None, min_am=None, min_spip=None, min_cadd=None, min_qual=None, pass_qc_only=False, consequences=None):
    total_start = len(df)
    logger.info(f"Starting filtering on dataset with {total_start} variants...")
    filtered_df = df.copy()

    # Track filter steps
    stats = {"Total_Input_Variants": total_start}

    # 0. Base Caller Quality Controls
    # Missing QC/QUAL data is never coerced to a value (e.g. 0) that would make it look
    # like a known-bad or known-good result. A variant with no recorded QC_STATUS/QUAL
    # always passes these opt-in filters (we cannot exclude on evidence we don't have),
    # but is explicitly flagged as QC-not-evaluated so it's never silently indistinguishable
    # from a variant that was actually checked and passed.
    if pass_qc_only:
        if "QC_STATUS" in filtered_df.columns:
            qc_missing_mask = filtered_df["QC_STATUS"].isna()
            filtered_df["QC_STATUS_evaluated"] = ~qc_missing_mask
            if qc_missing_mask.any():
                logger.warning(
                    f"{int(qc_missing_mask.sum())} variant(s) have no QC_STATUS value; "
                    "retained by --pass-qc-only (not excluded for lack of evidence) but "
                    "flagged QC_STATUS_evaluated=False for downstream review."
                )
            mask = (filtered_df["QC_STATUS"] != "LOW_QUAL") | qc_missing_mask
            filtered_df = filtered_df[mask]
            logger.info(f"Filter pass_qc_only: Retained {len(filtered_df)} / {total_start} variants.")
            stats["PASS_QC_Only"] = len(filtered_df)
            stats["QC_STATUS_missing_count"] = int(qc_missing_mask.sum())
        else:
            logger.warning(
                "QC_STATUS column not present in this dataset at all — --pass-qc-only "
                "cannot be evaluated for any variant; no variants excluded on this basis."
            )
            stats["QC_STATUS_column_present"] = False

    if min_qual is not None:
        if "QUAL" in filtered_df.columns:
            qual_vals = pd.to_numeric(filtered_df["QUAL"], errors="coerce")
            qual_missing_mask = qual_vals.isna()
            filtered_df["QUAL_evaluated"] = ~qual_missing_mask
            if qual_missing_mask.any():
                logger.warning(
                    f"{int(qual_missing_mask.sum())} variant(s) have no QUAL value; "
                    "retained by --min-qual (not excluded for lack of evidence, never "
                    "coerced to 0) but flagged QUAL_evaluated=False for downstream review."
                )
            mask = (qual_vals >= min_qual) | qual_missing_mask
            filtered_df = filtered_df[mask]
            logger.info(f"Filter min_qual >= {min_qual}: Retained {len(filtered_df)} variants.")
            stats[f"QUAL_>=_{min_qual}"] = len(filtered_df)
            stats["QUAL_missing_count"] = int(qual_missing_mask.sum())
        else:
            logger.warning(
                "QUAL column not present in this dataset at all — --min-qual cannot be "
                "evaluated for any variant; no variants excluded on this basis."
            )
            stats["QUAL_column_present"] = False

    # 1. gnomAD Allele Frequency Filter
    if max_af is not None:
        af_col = None
        for col in ['AF', 'gnomAD_AF', 'gnomAD_AF_joint', 'gnomad_af', 'AF_joint']:
            if col in filtered_df.columns:
                af_col = col
                break
        if af_col:
            af_vals = pd.to_numeric(filtered_df[af_col], errors='coerce').fillna(0.0)
            mask = af_vals <= max_af
            filtered_df = filtered_df[mask]
            logger.info(f"Filter max_af <= {max_af}: Retained {len(filtered_df)} / {total_start} variants.")
            stats[f"AF_<=_{max_af}"] = len(filtered_df)
        else:
            logger.warning("Allele frequency column not found in dataset. Skipping AF filter.")

    # 2. REVEL Score Filter
    if min_revel is not None:
        revel_col = None
        for col in ['REVEL_score', 'REVEL', 'revel']:
            if col in filtered_df.columns:
                revel_col = col
                break
        if revel_col:
            revel_vals = pd.to_numeric(filtered_df[revel_col], errors='coerce')
            mask = (revel_vals >= min_revel) | (revel_vals.isna())
            filtered_df = filtered_df[mask]
            logger.info(f"Filter min_revel >= {min_revel}: Retained {len(filtered_df)} variants.")
            stats[f"REVEL_>=_{min_revel}"] = len(filtered_df)

    # 3. AlphaMissense Score Filter
    if min_am is not None:
        am_col = None
        for col in ['am_pathogenicity', 'AlphaMissense', 'alphamissense']:
            if col in filtered_df.columns:
                am_col = col
                break
        if am_col:
            am_vals = pd.to_numeric(filtered_df[am_col], errors='coerce')
            mask = (am_vals >= min_am) | (am_vals.isna())
            filtered_df = filtered_df[mask]
            logger.info(f"Filter min_alphamissense >= {min_am}: Retained {len(filtered_df)} variants.")
            stats[f"AlphaMissense_>=_{min_am}"] = len(filtered_df)

    # 4. SPiP Score Filter
    if min_spip is not None:
        spip_col = None
        for col in ['SPiP_score', 'SPiP', 'spip_score']:
            if col in filtered_df.columns:
                spip_col = col
                break
        if spip_col:
            spip_vals = pd.to_numeric(filtered_df[spip_col], errors='coerce')
            mask = (spip_vals >= min_spip) | (spip_vals.isna())
            filtered_df = filtered_df[mask]
            logger.info(f"Filter min_spip >= {min_spip}: Retained {len(filtered_df)} variants.")
            stats[f"SPiP_>=_{min_spip}"] = len(filtered_df)

    # 5. CADD Score Filter
    if min_cadd is not None:
        cadd_col = None
        for col in ['CADD_PHRED', 'CADD', 'cadd_phred']:
            if col in filtered_df.columns:
                cadd_col = col
                break
        if cadd_col:
            cadd_vals = pd.to_numeric(filtered_df[cadd_col], errors='coerce')
            mask = (cadd_vals >= min_cadd) | (cadd_vals.isna())
            filtered_df = filtered_df[mask]
            logger.info(f"Filter min_cadd >= {min_cadd}: Retained {len(filtered_df)} variants.")
            stats[f"CADD_>=_{min_cadd}"] = len(filtered_df)

    # 6. Consequence Filter
    if consequences:
        c_list = [c.strip().lower() for c in consequences.split(',')]
        cons_col = None
        for col in ['Consequence', 'consequence', 'VEP_Consequence']:
            if col in filtered_df.columns:
                cons_col = col
                break
        if cons_col:
            mask = filtered_df[cons_col].astype(str).str.lower().apply(lambda val: any(c in val for c in c_list))
            filtered_df = filtered_df[mask]
            logger.info(f"Filter consequences ({consequences}): Retained {len(filtered_df)} variants.")
            stats["Consequence_Filtered"] = len(filtered_df)

    stats["Final_Retained_Variants"] = len(filtered_df)
    return filtered_df, pd.DataFrame([stats])

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    df = load_data(args.input)
    filtered_df, stats_df = filter_dataframe(
        df,
        max_af=args.max_af,
        min_revel=args.min_revel,
        min_am=args.min_alphamissense,
        min_spip=args.min_spip,
        min_cadd=args.min_cadd,
        min_qual=args.min_qual,
        pass_qc_only=args.pass_qc_only,
        consequences=args.consequences
    )

    out_ext = "pq" if args.output_format in ("pq", "parquet") else args.output_format
    name_prefix = f"{args.output_prefix}_" if args.output_prefix else ""
    out_file = os.path.join(args.output_dir, f"{name_prefix}filtered_variants.{out_ext}")
    stats_file = os.path.join(args.output_dir, f"{name_prefix}filtering_summary_metrics.tsv")

    if out_ext == "pq":
        filtered_df.to_parquet(out_file, index=False)
    elif out_ext == "xlsx":
        filtered_df.to_excel(out_file, index=False)
    else:
        filtered_df.to_csv(out_file, sep='\t', index=False)

    stats_df.to_csv(stats_file, sep='\t', index=False)
    logger.info(f"Saved filtered variants ({len(filtered_df)}) to: {out_file}")
    logger.info(f"Saved summary metrics to: {stats_file}")

if __name__ == "__main__":
    main()
