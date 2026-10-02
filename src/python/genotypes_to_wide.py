#!/usr/bin/env python3
"""Adds wide per-sample genotype columns (GT_<sample>, GQ_<sample>, DP_<sample>, AD_<sample>) to the
main variant table from the long genotype table, joined on `Locus`.

GT_<sample> uses the values the downstream filters read: HET, HOMALT, HOMREF, MISSING.
het -> HET; hom_alt, haploid_alt -> HOMALT; hom_ref, haploid_ref -> HOMREF; every other status
(no_call, partial_no_call, multiallelic_other) -> MISSING. The long table keeps the exact status.

Writes a NEW table; the input table is never modified. Fails when the sample count exceeds
`--max-samples` or when the join changes the row count.
"""

import argparse
import sys

import pandas as pd

GT_CATEGORY = {
    "het": "HET", "hom_alt": "HOMALT", "haploid_alt": "HOMALT",
    "hom_ref": "HOMREF", "haploid_ref": "HOMREF",
}


def read_table(path):
    if path.endswith((".pq", ".parquet")):
        return pd.read_parquet(path)
    if path.endswith((".tsv", ".tsv.gz")):
        return pd.read_csv(path, sep="\t", low_memory=False)
    raise ValueError(f"unsupported table format: {path} (use .pq or .tsv)")


def write_table(df, path):
    if path.endswith((".pq", ".parquet")):
        df.to_parquet(path, index=False)
    elif path.endswith((".tsv", ".tsv.gz")):
        df.to_csv(path, sep="\t", index=False)
    else:
        raise ValueError(f"unsupported table format: {path}")


def wide_genotypes(long_df):
    """One row per Locus with GT_/GQ_/DP_/AD_<sample> columns."""
    df = long_df.copy()
    df["GT"] = df["gt_status"].map(GT_CATEGORY).fillna("MISSING")
    df["GQ"], df["DP"], df["AD"] = df["gq"], df["dp"], df["ad"]
    wide = df.pivot(index="Locus", columns="sample_id", values=["GT", "GQ", "DP", "AD"])
    wide.columns = [f"{field}_{sample}" for field, sample in wide.columns]
    return wide.reset_index()


def run(table_path, genotypes_path, output_path, max_samples):
    table = read_table(table_path)
    long_df = pd.read_parquet(genotypes_path)
    n_samples = long_df["sample_id"].nunique()
    if n_samples > max_samples:
        print(f"ERROR: {n_samples} samples exceed --max-samples {max_samples}; use the long genotype table",
              file=sys.stderr)
        return 1
    if long_df.duplicated(["Locus", "sample_id"]).any():
        print("ERROR: duplicate Locus/sample_id pairs in the genotype table", file=sys.stderr)
        return 1
    wide = wide_genotypes(long_df)
    out = table.merge(wide, on="Locus", how="left", validate="many_to_one")
    if len(out) != len(table):
        print(f"ERROR: join changed the row count ({len(table)} -> {len(out)})", file=sys.stderr)
        return 1
    gt_cols = [c for c in out.columns if c.startswith("GT_")]
    out[gt_cols] = out[gt_cols].fillna("MISSING")
    unmatched = int(table["Locus"].isin(wide["Locus"]).eq(False).sum())
    write_table(out, output_path)
    print(f"Wide genotype view: {len(out)} rows, {n_samples} sample(s), {unmatched} row(s) without genotypes -> {output_path}")
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--table", required=True, help="Main variant table (.pq or .tsv)")
    p.add_argument("--genotypes", required=True, help="Long genotype table (.pq) from extract_genotypes.py")
    p.add_argument("--output", required=True, help="Output table path (must differ from --table)")
    p.add_argument("--max-samples", type=int, default=10, help="Fail above this many samples (default 10)")
    a = p.parse_args()
    if a.output == a.table:
        sys.exit("error: --output must differ from --table")
    sys.exit(run(a.table, a.genotypes, a.output, a.max_samples))


if __name__ == "__main__":
    main()
