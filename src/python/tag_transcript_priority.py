#!/usr/bin/env python3
"""
Adds TRANSCRIPT_PRIORITY_TIER to a VCF_TO_TABLE output table (whole-VCF entry point only —
the legacy per-gene path already filters to one transcript via --gene_set, so has no use for
this tagging step).

Tier 1: row's Feature (transcript) matches its gene's curated transcript in
        resources/gene_transcript_mapping.txt (the 217-panel's clinically-vetted choice).
Tier 2: row's own MANE_SELECT column (a standard VEP/CSQ field, already present when
        vcf_parser_pysam.py is run with --add_vep) is non-empty.
Tier 3: everything else — kept, not dropped, mirroring VEP's own --flag_pick semantics.

Transcript comparison ignores version suffixes (ENST00000545968.6 vs ENST00000545968), same
convention as src/hpc/gene_coords.sh's "drop version number from transcript".
"""

import argparse
import sys

import pandas as pd


def load_curated_mapping(path: str) -> dict[str, str]:
    """Returns {gene_symbol: transcript_id} from the 217-panel's Gen,NM,ENST file."""
    mapping = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            fields = line.strip().split(",")
            if len(fields) >= 3 and fields[0].strip():
                gene, enst = fields[0].strip().upper(), fields[2].strip()
                mapping[gene] = enst
    return mapping


def read_table(path: str) -> pd.DataFrame:
    if path.endswith(".pq") or path.endswith(".parquet"):
        return pd.read_parquet(path)
    if path.endswith(".tsv"):
        return pd.read_csv(path, sep="\t")
    raise ValueError(f"Unsupported table format: {path}")


def write_table(df: pd.DataFrame, path: str) -> None:
    if path.endswith(".pq") or path.endswith(".parquet"):
        df.to_parquet(path, index=False)
    elif path.endswith(".tsv"):
        df.to_csv(path, sep="\t", index=False)
    else:
        raise ValueError(f"Unsupported table format: {path}")


def _strip_version(transcript_id: str) -> str:
    return transcript_id.split(".")[0]


def tag_priority(df: pd.DataFrame, curated: dict[str, str]) -> pd.DataFrame:
    def tier_for_row(row) -> int:
        gene = str(row.get("SYMBOL", "")).strip().upper()
        feature = str(row.get("Feature", "")).strip()
        curated_transcript = curated.get(gene)
        if curated_transcript and feature and _strip_version(feature) == _strip_version(curated_transcript):
            return 1
        mane = row.get("MANE_SELECT", "")
        if pd.notna(mane) and str(mane).strip() not in ("", "."):
            return 2
        return 3

    df = df.copy()
    df["TRANSCRIPT_PRIORITY_TIER"] = df.apply(tier_for_row, axis=1)
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--gene-transcript-mapping", default="resources/gene_transcript_mapping.txt")
    args = parser.parse_args()

    df = read_table(args.input)
    curated = load_curated_mapping(args.gene_transcript_mapping)
    print(f"Loaded {len(curated)} curated gene->transcript mappings from {args.gene_transcript_mapping}")

    df = tag_priority(df, curated)
    tier_counts = df["TRANSCRIPT_PRIORITY_TIER"].value_counts().sort_index()
    print(f"Tier distribution: {dict(tier_counts)}")

    write_table(df, args.output)
    print(f"Wrote {len(df)} rows to {args.output}")


if __name__ == "__main__":
    sys.exit(main())
