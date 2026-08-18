#!/usr/bin/env python3
"""
Backfill / Transform already-generated Parquet files to include:
  1. Upgraded multi-omics NEW_IMPACT
  2. Clinical & Discovery PRIORITY_TIER (Tier 1 - 4)
  3. Continuous VARIANT_PRIORITY_SCORE (0.0 - 100.0)
"""

import os
import sys
import glob
import logging
import pandas as pd
import numpy as np
from datetime import datetime

# Import prioritization builder from filter_variants
sys.path.insert(0, "/home/mruizp/data_lab_PGP/shared/utils/src")
from filter_variants import build_newImpact, build_priority_tier

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("backfill_prioritization")

def update_parquet_prioritization(pq_path):
    try:
        df = pd.read_parquet(pq_path)
    except Exception as e:
        logger.error(f"Error reading {pq_path}: {e}")
        return False, 0

    # Apply upgraded NEW_IMPACT and PRIORITY_TIER
    df = build_newImpact(df)
    df = build_priority_tier(df)

    # Reorder columns to ensure priority metrics are prominent
    priority_cols = ["PRIORITY_TIER", "VARIANT_PRIORITY_SCORE", "NEW_IMPACT"]
    base_cols = [c for c in df.columns if c not in priority_cols]

    insert_idx = None
    for target in ["CDNA_NAME", "intron_offset_signed", "HGVSp", "HGVSc", "SYMBOL"]:
        if target in base_cols:
            insert_idx = base_cols.index(target) + 1
            break

    if insert_idx is not None:
        new_col_order = base_cols[:insert_idx] + priority_cols + base_cols[insert_idx:]
    else:
        new_col_order = priority_cols + base_cols

    df = df[new_col_order]

    # Atomic write
    tmp_out = pq_path + f".tmp_{os.getpid()}"
    try:
        df.to_parquet(tmp_out, index=False, engine="pyarrow", compression="snappy")
        os.replace(tmp_out, pq_path)
        logger.info(f"Successfully prioritized & updated {pq_path} ({len(df):,} variants, {len(df.columns)} columns)")
        return True, len(df)
    except Exception as e:
        if os.path.exists(tmp_out):
            os.remove(tmp_out)
        logger.error(f"Error saving prioritized parquet {pq_path}: {e}")
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

    logger.info(f"Found {len(all_files)} Parquet files to prioritize across runs.")
    
    success_count = 0
    total_variants = 0
    for f in sorted(all_files):
        ok, n_vars = update_parquet_prioritization(f)
        if ok:
            success_count += 1
            total_variants += n_vars
            
    logger.info(f"Prioritization backfill complete: {success_count}/{len(all_files)} files updated successfully ({total_variants:,} total variants prioritized).")

if __name__ == "__main__":
    main()
