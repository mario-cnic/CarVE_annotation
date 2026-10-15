#!/usr/bin/env python3
"""Compare two Parquet tables: schema (types), footer metadata, row groups, column order, values, row order.

Usage: compare_tables.py A.pq B.pq [--ignore-footer-key carve.provenance]
Exit 0 when the tables match, 1 otherwise. The provenance stamp differs between launches, so ignore its key
when comparing two runs.
"""
import argparse
import sys

import pyarrow.parquet as pq

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("a")
ap.add_argument("b")
ap.add_argument("--ignore-footer-key", action="append", default=[])
args = ap.parse_args()

fa, fb = pq.ParquetFile(args.a), pq.ParquetFile(args.b)
skip = {k.encode() for k in args.ignore_footer_key}
# The serialized ARROW:schema embeds the other keys; with keys ignored, compare the decoded schema metadata.
if skip:
    skip.add(b"ARROW:schema")
meta_a = {k: v for k, v in (fa.metadata.metadata or {}).items() if k not in skip}
meta_b = {k: v for k, v in (fb.metadata.metadata or {}).items() if k not in skip}
meta_a["schema"] = {k: v for k, v in (fa.schema_arrow.metadata or {}).items() if k not in skip}
meta_b["schema"] = {k: v for k, v in (fb.schema_arrow.metadata or {}).items() if k not in skip}
groups = lambda f: [f.metadata.row_group(i).num_rows for i in range(f.metadata.num_row_groups)]
da, db = fa.read().to_pandas(), fb.read().to_pandas()
checks = {
    "schema": fa.schema_arrow.remove_metadata().equals(fb.schema_arrow.remove_metadata()),
    "footer": meta_a == meta_b,
    "row_groups": groups(fa) == groups(fb),
    "column_order": list(da.columns) == list(db.columns),
    "rows": len(da) == len(db),
}
differing = [c for c in da.columns if c in db.columns and not da[c].equals(db[c])] if checks["rows"] else []
checks["values"] = checks["rows"] and not differing
for name, ok in checks.items():
    print(f"{name:13s} {'same' if ok else 'DIFFERENT'}")
if differing:
    print("differing columns:", ", ".join(differing[:20]))
sys.exit(0 if all(checks.values()) else 1)
