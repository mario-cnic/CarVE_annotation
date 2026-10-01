#!/usr/bin/env python3
"""Derive which dbNSFP columns are lists aligned position-by-position with `Ensembl_transcriptid`.

dbNSFP/VEP write per-transcript values as '&'-joined lists attached to EVERY transcript row of a variant.
A column is called transcript-aligned when, on rows where it is multi-valued, its item count equals the
item count of `Ensembl_transcriptid` on >= --min-fraction of those rows (and there are >= --min-rows such
rows). Output: one column name per line (resources/dbnsfp_transcript_aligned_columns.txt) + manifest JSON.
Read-only on the input table.
"""
import argparse, hashlib, json, sys
from datetime import datetime, timezone
import pandas as pd, pyarrow.parquet as pq

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--pq", required=True); ap.add_argument("--row-groups", default="0")
ap.add_argument("--min-rows", type=int, default=1000); ap.add_argument("--min-fraction", type=float, default=0.999)
ap.add_argument("--out", required=True)
a = ap.parse_args()

f = pq.ParquetFile(a.pq)
frames = [f.read_row_group(int(g)).to_pandas() for g in a.row_groups.split(",")]
t = pd.concat(frames, ignore_index=True)
db = t["Ensembl_transcriptid"].notna()
n_t = t.loc[db, "Ensembl_transcriptid"].astype(str).str.count("&") + 1
aligned, report = [], {}
for c in t.columns:
    if c == "Ensembl_transcriptid" or not t[c].notna().any() or (t[c].notna() & ~db).any():
        continue
    s = t.loc[db, c]; nn = s.notna()
    cnt = s[nn].astype(str).str.count("&") + 1
    multi = cnt > 1
    if multi.sum() < a.min_rows:
        continue
    frac = float((cnt[multi] == n_t[nn][multi]).mean())
    report[c] = {"multi_rows": int(multi.sum()), "fraction_len_equals_transcripts": round(frac, 4)}
    if frac >= a.min_fraction:
        aligned.append(c)
aligned = sorted(["Ensembl_transcriptid"] + aligned)
open(a.out, "w").write("\n".join(aligned) + "\n")
sha = hashlib.sha256(open(a.pq, "rb").read(1 << 20)).hexdigest()
json.dump({"derived": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
           "input": a.pq, "row_groups": a.row_groups, "dbnsfp_rows": int(db.sum()), "min_rows": a.min_rows,
           "min_fraction": a.min_fraction, "n_aligned_columns": len(aligned),
           "multi_valued_but_not_aligned": sorted(c for c, r in report.items() if c not in aligned),
           "argv": sys.argv}, open(a.out + ".manifest.json", "w"), indent=1)
print(len(aligned), "aligned columns;", "multi-valued but NOT aligned:", sorted(c for c in report if c not in aligned))
