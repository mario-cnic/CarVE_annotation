#!/usr/bin/env python3
"""Count how many rows of a parsed table survive successive candidate flags.

Reads five columns only, one row group at a time. Rows are variant x transcript,
not variants. Thresholds are illustrative, not validated for any disease.
"""
import argparse
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

AF_COLUMN = "gnomADv4_AF_grpmax_joint"
AF_MAX = 0.001
IMPACTS = ["HIGH", "MODERATE"]
MAX_TIER = 2
COLUMNS = ["Locus", "TRANSCRIPT_PRIORITY_TIER", AF_COLUMN, "IMPACT"]


def count(mask):
    """Number of true values in a boolean array (nulls count as false)."""
    return int(pc.sum(pc.cast(pc.fill_null(mask, False), pa.int64())).as_py() or 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("parquet", help="parsed.clean.pq table")
    args = ap.parse_args()

    pf = pq.ParquetFile(args.parquet)
    steps = ["all_rows", "tier_le_max", "tier_af", "tier_af_impact"]
    totals = dict.fromkeys(steps, 0)
    loci = []
    for rg in range(pf.num_row_groups):
        t = pf.read_row_group(rg, columns=COLUMNS)
        af = t[AF_COLUMN]
        tier_ok = pc.less_equal(t["TRANSCRIPT_PRIORITY_TIER"], MAX_TIER)
        # a missing frequency is kept: absent from gnomAD is not the same as common
        af_ok = pc.or_(pc.is_null(af), pc.less_equal(af, AF_MAX))
        impact_ok = pc.is_in(t["IMPACT"], value_set=pa.array(IMPACTS))
        step_masks = {
            "tier_le_max": tier_ok,
            "tier_af": pc.and_(tier_ok, af_ok),
            "tier_af_impact": pc.and_(pc.and_(tier_ok, af_ok), impact_ok),
        }
        totals["all_rows"] += t.num_rows
        for name, mask in step_masks.items():
            totals[name] += count(mask)
        loci.append(pc.unique(t["Locus"]))

    n_loci = len(pc.unique(pa.chunked_array(loci)))
    print(f"distinct_Locus\t{n_loci}")
    for name in steps:
        print(f"{name}\t{totals[name]}\t{100 * totals[name] / totals['all_rows']:.3f}%")


if __name__ == "__main__":
    main()
