#!/usr/bin/env python3
"""Audit gene/transcript attribution of predictor scores in a parsed annotation table.

Question: when a variant overlaps (or lies near) more than one gene, does the score placed on
a row (one row = one variant x one VEP transcript) come from that row's own gene/transcript?

Read-only. Reads a few columns of a `*.parsed.clean.pq` in row-group batches and writes
aggregate counts only (no loci, no identifiers) as JSON.

Mirrors the parser logic in shared/utils/src/modules/splicing.py (read 2026-10-01):
  - parse_spliceai_custom: uses only the FIRST entry of the SpliceAI string.
  - parse_pangolin:        max |score| over every gene entry; gene id never read.
  - parse_spip:            first entry whose gene symbol == row SYMBOL, else first entry.
Multi-entry values are joined with '&' in the table (verified: pipes-per-string of 18, 27, ...).

Thresholds are reporting thresholds only (0.20 = the splicing floor used by
generate_clinical_prioritization_report.py:833 and parse_spip's mechanism rule).
"""
import argparse
import json
import re
import sys
from collections import Counter

import pyarrow.parquet as pq

FLOOR = 0.20
SPLIT = re.compile(r"[&]")


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def spliceai_entries(s):
    out = []
    for e in SPLIT.split(s):
        p = e.split("|")
        if len(p) < 10:
            continue
        sc = [v for v in (_f(p[2]), _f(p[3]), _f(p[4]), _f(p[5])) if v is not None]
        out.append((p[1], max(sc) if sc else None))
    return out


def pangolin_entries(s):
    out = []
    for e in SPLIT.split(s):
        p = e.split("|")
        gene = p[0].split(",")[0]
        sc = []
        for x in p[1:]:
            if ":" in x:
                a = x.split(":")
                if len(a) == 2 and _f(a[1]) is not None and 0 <= abs(_f(a[1])) <= 1:
                    sc.append(abs(_f(a[1])))
        out.append((gene, max(sc) if sc else None))
    return out


def spip_entries(s):
    out = []
    for e in re.split(r"[,&]", s):
        p = e.split("|")
        if len(p) > 12:
            out.append((p[11].strip().split(".")[0], p[12].strip(), _f(p[4].strip())))
    return out


def crosses(a, b):
    """True when two scores fall on different sides of the reporting floor."""
    if a is None or b is None:
        return a != b
    return (a >= FLOOR) != (b >= FLOOR)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--pq", required=True)
    ap.add_argument("--row-groups", default="all", help="comma list or 'all'")
    ap.add_argument("--batch", type=int, default=100_000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    f = pq.ParquetFile(a.pq)
    rgs = range(f.num_row_groups) if a.row_groups == "all" else [int(x) for x in a.row_groups.split(",")]
    cols = ["SYMBOL", "Gene", "Feature", "MANE_SELECT", "SpliceAI", "Pangolin", "SPiP"]
    c = Counter()
    for rg in rgs:
        for b in f.iter_batches(batch_size=a.batch, row_groups=[rg], columns=cols):
            d = b.to_pydict()
            for sym, gene, _ft, mane, sai, pan, spip in zip(*(d[k] for k in cols)):
                c["rows"] += 1
                # ---------------- SpliceAI ----------------
                if sai:
                    ent = spliceai_entries(sai)
                    if ent:
                        c["sai_rows"] += 1
                        used = ent[0]
                        if len(ent) > 1:
                            c["sai_multi_entry_rows"] += 1
                        symset = {g for g, _ in ent}
                        if sym is None:
                            c["sai_row_no_symbol"] += 1
                        elif used[0] == sym:
                            c["sai_used_entry_is_row_gene"] += 1
                        elif sym in symset:
                            c["sai_WRONG_entry_used__correct_entry_available"] += 1
                            right = next(s for g, s in ent if g == sym)
                            if crosses(used[1], right):
                                c["sai_WRONG_entry__floor_flip"] += 1
                        else:
                            c["sai_row_gene_not_scored__neighbour_score_attached"] += 1
                            if (used[1] or 0) >= FLOOR:
                                c["sai_neighbour_score_ge_floor"] += 1
                # ---------------- Pangolin ----------------
                if pan:
                    ent = pangolin_entries(pan)
                    c["pan_rows"] += 1
                    if len(ent) > 1:
                        c["pan_multi_gene_rows"] += 1
                        sc = [s for _, s in ent if s is not None]
                        if len(sc) > 1 and crosses(max(sc), min(sc)):
                            c["pan_multi_gene__floor_depends_on_which_gene"] += 1
                    m = re.match(r"ENS([A-Z]*)G", ent[0][0])
                    c["pan_ids_" + ("mouse" if m and m.group(1) == "MUS" else "human" if m else "other")] += 1
                # ---------------- SPiP ----------------
                if spip and "caused an error" not in spip:
                    ent = spip_entries(spip)
                    if ent:
                        c["spip_rows"] += 1
                        genes = {g for _, g, _ in ent}
                        mine = [e for e in ent if e[1] == sym] if sym else []
                        if len(genes) > 1:
                            c["spip_multi_gene_rows"] += 1
                        if sym is not None and not mine:
                            c["spip_row_gene_absent__first_entry_of_other_gene_used"] += 1
                            if (ent[0][2] or 0) >= FLOOR:
                                c["spip_fallback_score_ge_floor"] += 1
                        if len({t for t, _, _ in mine}) > 1:
                            c["spip_row_gene_has_multiple_transcripts"] += 1
                            used = mine[0]
                            if mane:
                                m = str(mane).split(".")[0]
                                c["spip_multi_tx_rows_with_mane"] += 1
                                hit = [e for e in mine if e[0] == m]
                                if hit and hit[0][0] != used[0]:
                                    c["spip_WRONG_transcript_used__MANE_available"] += 1
                                    if crosses(used[2], hit[0][2]):
                                        c["spip_WRONG_transcript__floor_flip"] += 1
                                elif not hit:
                                    c["spip_mane_not_among_spip_entries"] += 1
                elif spip:
                    c["spip_error_rows"] += 1
        print(f"row group {rg} done", file=sys.stderr, flush=True)

    json.dump({"pq": a.pq, "row_groups": list(rgs), "floor": FLOOR, "counts": dict(sorted(c.items()))},
              open(a.out, "w"), indent=2)
    print(json.dumps(dict(sorted(c.items())), indent=1))


if __name__ == "__main__":
    main()
