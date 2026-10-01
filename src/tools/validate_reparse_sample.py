#!/usr/bin/env python3
"""Re-parse the raw predictor strings of a parsed table with the gene/transcript-matched parsers and
compare with the table's existing (old-parser) columns. Read-only on the input; writes aggregates only.

Usage: validate_reparse_sample.py --pq T.parsed.clean.pq --row-group 16 --hgnc H.txt --symbol-map M.tsv \
          --enst-spip-map resources/enst_to_spip_nm.grch38.gencode45.tsv --out result.json
"""
import argparse, json, os, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python", "parsing"))
from modules import splicing as sp, gene_identity as gi  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--pq", required=True); ap.add_argument("--row-group", type=int, required=True)
ap.add_argument("--hgnc", required=True); ap.add_argument("--symbol-map", required=True)
ap.add_argument("--enst-spip-map", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()

need = ["Locus", "Gene", "SYMBOL", "Feature", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SpliceAI", "SPiP", "Pangolin",
        "spliceai_custom_MAX", "SPiP_prediction", "Pangolin_max_score"]
df = pq.ParquetFile(a.pq).read_row_group(a.row_group, columns=need).to_pandas()
old = df[["spliceai_custom_MAX", "SPiP_prediction", "Pangolin_max_score"]].copy()
ident = gi.GeneIdentity(a.symbol_map, a.hgnc, a.enst_spip_map)
df = df.drop(columns=["spliceai_custom_MAX", "SPiP_prediction", "Pangolin_max_score"])
df = sp.parse_spliceai_custom(df, identity=ident)
df = sp.parse_spip(df, identity=ident)
df = sp.parse_pangolin(df)

FLOOR = 0.20
res = {"rows": len(df)}
for name, mcol, ncol, ocol in (("spliceai", "spliceai_custom_match", "spliceai_custom_MAX", "spliceai_custom_MAX"),
                               ("spip", "SPiP_match", "SPiP_prediction", "SPiP_prediction"),
                               ("pangolin", "Pangolin_match", "Pangolin_max_score", "Pangolin_max_score")):
    new, o = df[ncol], old[ocol]
    present_old = o.notna()
    res[name] = {
        "match_status": df[mcol].value_counts().to_dict(),
        "old_present": int(present_old.sum()), "new_present": int(new.notna().sum()),
        "old_present_new_absent": int((present_old & new.isna()).sum()),
        "old_absent_new_present": int((o.isna() & new.notna()).sum()),
        "value_changed_both_present": int(((o != new) & present_old & new.notna()).sum()),
        f"floor_{FLOOR}_flip_both_present": int((((o >= FLOOR) != (new >= FLOOR)) & present_old & new.notna()).sum()),
        "old_ge_floor_new_absent": int(((o >= FLOOR) & new.isna()).sum()),
    }
# independent check: every alias_resolved SpliceAI row must agree with the coordinate table
m = pd.read_csv(a.symbol_map, sep="\t").set_index("symbol")["ensg"]
al = df[df["spliceai_custom_match"] == gi.ALIAS_RESOLVED]
res["spliceai_alias_rows_checked"] = int(len(al))
res["spliceai_alias_rows_coordinate_disagree"] = int((al["spliceai_custom_SYMBOL"].map(m) != al["Gene"].str.split(".").str[0]).sum())
json.dump(res, open(a.out, "w"), indent=1)
print(json.dumps(res, indent=1))
