#!/usr/bin/env python3
"""
Backfill / Transform already-generated Parquet files to include unpacked `spliceai_custom_*` columns.
"""

import os
import sys
import glob
import logging
import pandas as pd
import numpy as np
from datetime import datetime

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("backfill_spliceai_custom")

def extract_fields(val):
    if pd.isna(val) or val is None or str(val).strip() in ("", ".", "nan", "None", "-"):
        return [np.nan] * 10
    first_entry = str(val).split(",")[0].strip()
    parts = first_entry.split("|")
    if len(parts) >= 10:
        try:
            symbol = parts[1] if parts[1] != "." else np.nan
            ds_ag = float(parts[2]) if parts[2] != "." else np.nan
            ds_al = float(parts[3]) if parts[3] != "." else np.nan
            ds_dg = float(parts[4]) if parts[4] != "." else np.nan
            ds_dl = float(parts[5]) if parts[5] != "." else np.nan
            dp_ag = float(parts[6]) if parts[6] != "." else np.nan
            dp_al = float(parts[7]) if parts[7] != "." else np.nan
            dp_dg = float(parts[8]) if parts[8] != "." else np.nan
            dp_dl = float(parts[9]) if parts[9] != "." else np.nan
            scores = [x for x in [ds_ag, ds_al, ds_dg, ds_dl] if not np.isnan(x)]
            max_score = max(scores) if scores else np.nan
            return [symbol, ds_ag, ds_al, ds_dg, ds_dl, dp_ag, dp_al, dp_dg, dp_dl, max_score]
        except Exception:
            return [np.nan] * 10
    return [np.nan] * 10

def update_parquet_file(pq_path):
    try:
        df = pd.read_parquet(pq_path)
    except Exception as e:
        logger.error(f"Error reading {pq_path}: {e}")
        return False, 0

    if "SpliceAI" not in df.columns:
        logger.warning(f"No 'SpliceAI' column found in {pq_path}, skipping.")
        return False, 0

    if "spliceai_custom_DS_AG" in df.columns and "spliceai_custom_MAX" in df.columns:
        logger.info(f"Already transformed: {pq_path}")
        return True, len(df)

    custom_cols = [
        "spliceai_custom_SYMBOL",
        "spliceai_custom_DS_AG",
        "spliceai_custom_DS_AL",
        "spliceai_custom_DS_DG",
        "spliceai_custom_DS_DL",
        "spliceai_custom_DP_AG",
        "spliceai_custom_DP_AL",
        "spliceai_custom_DP_DG",
        "spliceai_custom_DP_DL",
        "spliceai_custom_MAX",
    ]

    extracted = pd.DataFrame(
        df["SpliceAI"].apply(extract_fields).tolist(),
        index=df.index,
        columns=custom_cols
    )

    for c in custom_cols:
        df[c] = extracted[c]

    # Reorder columns to put custom spliceai columns near spliceAI_MAX
    cols = list(df.columns)
    # Remove newly added cols from end
    base_cols = [c for c in cols if c not in custom_cols]
    
    # Locate spliceAI_MAX or SpliceAI_status position
    insert_idx = None
    for target in ["SpliceAI_status", "spliceAI_MAX", "IMPACT"]:
        if target in base_cols:
            insert_idx = base_cols.index(target) + 1
            break
            
    if insert_idx is not None:
        new_col_order = base_cols[:insert_idx] + custom_cols + base_cols[insert_idx:]
    else:
        new_col_order = base_cols + custom_cols
        
    df = df[new_col_order]

    # Atomic write
    tmp_out = pq_path + f".tmp_{os.getpid()}"
    try:
        df.to_parquet(tmp_out, index=False, engine="pyarrow", compression="snappy")
        os.replace(tmp_out, pq_path)
        logger.info(f"Successfully transformed and updated {pq_path} ({len(df):,} variants, {len(df.columns)} columns)")
        return True, len(df)
    except Exception as e:
        if os.path.exists(tmp_out):
            os.remove(tmp_out)
        logger.error(f"Error saving updated parquet {pq_path}: {e}")
        return False, 0

def main():
    target_dirs = [
        "RUNS/run_20260813_1028/results",
        "RUNS/run_more_genes_20260817_1143/results"
    ]
    
    all_files = []
    for d in target_dirs:
        if os.path.exists(d):
            all_files.extend(glob.glob(os.path.join(d, "*.parsed.clean.pq")))

    logger.info(f"Found {len(all_files)} Parquet files to process across runs.")
    
    success_count = 0
    total_variants = 0
    for f in sorted(all_files):
        ok, n_vars = update_parquet_file(f)
        if ok:
            success_count += 1
            total_variants += n_vars
            
    logger.info(f"Transformation complete: {success_count}/{len(all_files)} files updated successfully ({total_variants:,} total variants processed).")

if __name__ == "__main__":
    main()
