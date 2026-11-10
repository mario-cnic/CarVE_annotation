#!/usr/bin/env python3
"""Compare UTRAnnotator output between two runs of the same input, and against two older runs.

Checks (walkthrough/20261110_utrannotator_deterministic_fork.md):
  - every row whose old runs disagreed has a non-empty 5UTR_consequence and `consequence=`-tagged entries;
  - 5UTR_consequence, 5UTR_annotation and UTRAnnotator_status are identical between run1 and run2;
  - on each such row the new 5UTR_consequence contains both classes the old runs reported.
Exit code 0 when all pass.
"""
import argparse
import sys

import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds

KEY = ["Locus", "Feature"]
UTR = ["5UTR_consequence", "5UTR_annotation", "UTRAnnotator_status"]
EMPTY = {"", "-", ".", "nan", "None"}


def load(path: str, loci: set[str] | None) -> pd.DataFrame:
    d = ds.dataset(path)
    cols = [c for c in KEY + UTR if c in d.schema.names]
    flt = pc.field("Locus").isin(sorted(loci)) if loci else pc.field("5UTR_consequence").is_valid()
    df = d.to_table(columns=cols, filter=flt).to_pandas()
    return df.astype(object).where(df.notna(), "").astype(str)


def terms(value: str) -> set[str]:
    return {t for t in value.replace(",", "&").split("&") if t.strip() and t not in EMPTY}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run1", required=True, help="table (.pq) of the first new run")
    ap.add_argument("--run2", required=True, help="table (.pq) of the second new run")
    ap.add_argument("--old1", required=True, help="table (.pq) of an older run of the same variants")
    ap.add_argument("--old2", required=True, help="table (.pq) of another older run")
    a = ap.parse_args()

    old = load(a.old1, None).merge(load(a.old2, None), on=KEY, suffixes=("_o1", "_o2"))
    old = old[old["5UTR_consequence_o1"] != old["5UTR_consequence_o2"]]
    loci = set(old.Locus)
    print(f"old runs disagree on {len(old)} rows, {len(loci)} variants")

    r1, r2 = load(a.run1, loci), load(a.run2, loci)
    ok = True
    both = r1.merge(r2, on=KEY, how="outer", suffixes=("_1", "_2"), indicator=True)
    if (both._merge != "both").any():
        print(f"FAIL: {int((both._merge != 'both').sum())} rows only in one run")
        ok = False
    for col in UTR:
        diff = both[(both._merge == "both") & (both[f"{col}_1"] != both[f"{col}_2"])]
        print(f"{col}: {len(diff)} rows differ between run1 and run2")
        ok &= diff.empty

    cmp = old.merge(r1, on=KEY, how="left")
    missing = cmp["5UTR_consequence"].isna() | cmp["5UTR_consequence"].isin(EMPTY)
    untagged = ~cmp["5UTR_annotation"].fillna("").str.startswith("consequence=")
    lost = [
        (row.Locus, row.Feature) for _, row in cmp.iterrows()
        if not (terms(row["5UTR_consequence_o1"]) | terms(row["5UTR_consequence_o2"])) <= terms(str(row["5UTR_consequence"]))
    ]
    print(f"rows without 5UTR_consequence in run1: {int(missing.sum())}; without consequence= tag: {int(untagged.sum())}")
    print(f"rows where run1 lacks a class an old run reported: {len(lost)}")
    for locus, feature in lost[:10]:
        print(f"  {locus} {feature}")
    print(f"UTRAnnotator_status in run1 on these rows: {cmp.UTRAnnotator_status.value_counts().to_dict()}")
    ok &= not missing.any() and not untagged.any() and not lost
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
