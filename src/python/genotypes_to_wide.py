#!/usr/bin/env python3
"""Adds wide per-sample genotype columns (GT_<sample>, GQ_<sample>, DP_<sample>, AD_<sample>) to the
main variant table from the long genotype table, joined on `Locus`.

GT_<sample> uses the values the downstream filters read: HET, HOMALT, HOMREF, MISSING.
het, alt_plus_other_allele (one copy of this ALT) -> HET; hom_alt, haploid_alt -> HOMALT; hom_ref, haploid_ref -> HOMREF; every other status
(no_call, partial_no_call, multiallelic_other) -> MISSING. The long table keeps the exact status.

Writes a NEW Parquet table; the input table is never modified. Fails (exit 1) when the sample count
exceeds `--max-samples`, when the join changes the row count, or when any table row has no genotype
row (the output is then incomplete).
"""

import argparse
import sys

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

GT_CATEGORY = {
    "het": "HET", "alt_plus_other_allele": "HET", "hom_alt": "HOMALT", "haploid_alt": "HOMALT",
    "hom_ref": "HOMREF", "haploid_ref": "HOMREF",
}


def wide_genotypes(long_df):
    """One row per Locus with GT_/GQ_/DP_/AD_<sample> columns."""
    df = long_df[["Locus", "sample_id", "gt_status", "gq", "dp", "ad"]].copy()
    df["GT"] = df["gt_status"].map(GT_CATEGORY).fillna("MISSING")
    df["GQ"], df["DP"], df["AD"] = df["gq"], df["dp"], df["ad"]
    wide = df.pivot(index="Locus", columns="sample_id", values=["GT", "GQ", "DP", "AD"])
    wide.columns = [f"{field}_{sample}" for field, sample in wide.columns]
    return wide


def run(table_path, genotypes_path, output_path, max_samples, batch_rows=250000):
    if not table_path.endswith((".pq", ".parquet")) or not output_path.endswith((".pq", ".parquet")):
        print("ERROR: the wide view supports Parquet tables only (--output_format pq)", file=sys.stderr)
        return 1
    long_df = pd.read_parquet(genotypes_path, columns=["Locus", "sample_id", "gt_status", "gq", "dp", "ad"])
    n_samples = long_df["sample_id"].nunique()
    if n_samples > max_samples:
        print(f"ERROR: {n_samples} samples exceed --max-samples {max_samples}; use the long genotype table",
              file=sys.stderr)
        return 1
    if long_df.duplicated(["Locus", "sample_id"]).any():
        print("ERROR: duplicate Locus/sample_id pairs in the genotype table", file=sys.stderr)
        return 1
    wide = wide_genotypes(long_df)
    del long_df
    gt_cols = [c for c in wide.columns if c.startswith("GT_")]

    # The main table is joined batch by batch: a WGS table has tens of millions of rows x ~600 columns.
    source = pq.ParquetFile(table_path)
    writer, n_rows, unmatched = None, 0, 0
    try:
        for batch in source.iter_batches(batch_size=batch_rows):
            part = batch.to_pandas()
            joined = part.join(wide, on="Locus", how="left")
            unmatched += int(joined[gt_cols[0]].isna().sum()) if gt_cols else 0
            joined[gt_cols] = joined[gt_cols].fillna("MISSING")
            table = pa.Table.from_pandas(joined, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(output_path, table.schema)
            writer.write_table(table.cast(writer.schema))
            n_rows += len(joined)
    finally:
        if writer is not None:
            writer.close()
    if n_rows != source.metadata.num_rows:
        print(f"ERROR: joined {n_rows} rows, the table has {source.metadata.num_rows}", file=sys.stderr)
        return 1
    if unmatched:
        print(f"ERROR: {unmatched} row(s) of the table have no genotype row (Locus mismatch); "
              f"{output_path} is incomplete and must not be used", file=sys.stderr)
        return 1
    print(f"Wide genotype view: {n_rows} rows, {n_samples} sample(s) -> {output_path}")
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--table", required=True, help="Main variant table (.pq)")
    p.add_argument("--genotypes", required=True, help="Long genotype table (.pq) from extract_genotypes.py")
    p.add_argument("--output", required=True, help="Output table path (must differ from --table)")
    p.add_argument("--max-samples", type=int, default=10, help="Fail above this many samples (default 10)")
    a = p.parse_args()
    if a.output == a.table:
        sys.exit("error: --output must differ from --table")
    sys.exit(run(a.table, a.genotypes, a.output, a.max_samples))


if __name__ == "__main__":
    main()
