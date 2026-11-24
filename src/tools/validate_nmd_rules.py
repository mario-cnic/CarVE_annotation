#!/usr/bin/env python3
"""Check parse_nmd against the VEP NMD plugin on an existing parsed table, without re-running VEP.

1. Reproduces the plugin's rules at the variant position from the transcript exon table and compares the result
   with the table's `NMD` flag (stop_gained/frameshift rows). Disagreement points at a coordinate or model bug.
   Plugin rules (NMD.pm): intronless transcript; variant CDS end <= 101; variant genomic end in the last exon or
   within 52 nt (inclusive window) of the penultimate exon's 3' end.
2. Counts how NMD_STATUS moves from the table's values to the new parse_nmd result, by consequence and NMD_RULE.

Reads only the needed columns, one row group at a time. Writes aggregate counts (JSON) to --out; with --rows, also
the per-row comparison (contains sample variants: keep it out of the repository).
"""
import argparse, json, sys
from collections import Counter
from pathlib import Path

import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from parsing.modules.predictors import _first_int, load_transcript_exons, parse_nmd  # noqa: E402

COLUMNS = ["Feature", "Consequence", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "STRAND", "VARIANT_CLASS",
           "NMD", "NMD_STATUS", "MANE_SELECT"]
SCOPE = r"stop_gained|frameshift_variant|splice_donor_variant|splice_acceptor_variant"


def plugin_flag(row, exons) -> str:
    """'escape', 'none' or why the plugin's rules cannot be reproduced for this row."""
    model = exons.get(str(row.Feature).split(".")[0])
    if model is None:
        return "no_model"
    ends = model[1]
    cds_txt, cdna_txt = str(row.CDS_position), str(row.cDNA_position)
    if len(ends) == 1:
        return "escape"
    cds_end = _first_int(cds_txt.split("-")[-1]) if cds_txt not in ("", "None", "nan") else None
    cds_start = _first_int(cds_txt)
    if cds_start is None or cds_end is None:
        return "none"  # plugin returns nothing without a CDS position
    if cds_end <= 101:
        return "escape"
    parts = [_first_int(x) for x in cdna_txt.split("-")]
    if any(p is None for p in parts):
        return "cdna_unknown"
    a, b = parts[0], parts[-1]
    insertion = str(row.VARIANT_CLASS) == "insertion"
    plus = float(row.STRAND) == 1.0
    # cDNA of the variant's genomic end (Ensembl insertions have end = start - 1)
    vf_end = (a if insertion else b) if plus else (b if insertion else a)
    pen = ends[-2]
    if vf_end > pen or pen - 51 <= vf_end <= pen:
        return "escape"
    return "none"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="parsed .pq table with the VEP NMD column")
    ap.add_argument("--exons", required=True, help="transcript exon table (build_transcript_exon_table.py)")
    ap.add_argument("--out", required=True, help="output JSON with aggregate counts")
    ap.add_argument("--rows", help="optional TSV with the per-row comparison")
    args = ap.parse_args(argv)

    exons = load_transcript_exons(args.exons)
    pf = pq.ParquetFile(args.table)
    parts = []
    for i in range(pf.num_row_groups):
        t = pf.read_row_group(i, columns=COLUMNS)
        t = t.filter(pc.match_substring_regex(pc.cast(t["Consequence"], "string"), SCOPE))
        parts.append(t.to_pandas())
    d = pd.concat(parts, ignore_index=True)
    old = d["NMD_STATUS"].astype(str)
    new = parse_nmd(d.drop(columns=["NMD_STATUS"]).copy(), exons)
    d["NMD_STATUS_new"], d["NMD_RULE"], d["NMD_PTC_BASIS"] = new.NMD_STATUS, new.NMD_RULE, new.NMD_PTC_BASIS
    trunc = d.Consequence.astype(str).str.contains("stop_gained|frameshift_variant")
    d["kind"] = d.Consequence.astype(str).map(
        lambda c: "frameshift" if "frameshift_variant" in c else "stop_gained" if "stop_gained" in c else "splice")
    d["plugin_table"] = d.NMD.fillna("").astype(str).str.contains("NMD_escaping_variant").map({True: "escape", False: "none"})
    d.loc[trunc, "plugin_repro"] = [plugin_flag(r, exons) for r in d.loc[trunc].itertuples(index=False)]

    def counts(frame, cols):
        return [dict(zip(cols + ["n"], [*(k if isinstance(k, tuple) else (k,)), int(v)]))
                for k, v in frame.groupby(cols, dropna=False).size().items()]

    d["NMD_STATUS_old"] = old
    d["mane_select"] = d.MANE_SELECT.notna() & (d.MANE_SELECT.astype(str) != "")
    summary = {
        "table": args.table,
        "exons": args.exons,
        "rows_in_scope": int(len(d)),
        "plugin_reproduction": counts(d.loc[trunc], ["kind", "plugin_table", "plugin_repro"]),
        "status_moves": counts(d, ["mane_select", "kind", "NMD_STATUS_old", "NMD_STATUS_new", "NMD_RULE"]),
    }
    Path(args.out).write_text(json.dumps(summary, indent=1))
    if args.rows:
        d.to_csv(args.rows, sep="\t", index=False)
    agree = d.loc[trunc & d.plugin_repro.isin(["escape", "none"])]
    print(f"rows in scope {len(d)}; plugin reproduced on {len(agree)} rows, agreement "
          f"{(agree.plugin_table == agree.plugin_repro).mean():.4%}")
    print(pd.crosstab([d.kind, d.NMD_STATUS_old], [d.NMD_STATUS_new]).to_string())


if __name__ == "__main__":
    main()
