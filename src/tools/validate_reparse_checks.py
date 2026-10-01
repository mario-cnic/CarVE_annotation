#!/usr/bin/env python3
"""Extra checks of the gene/transcript-matched parsers on real row groups (read-only, aggregates only):
 (a) match status by transcript tier (1 = curated clinical ENST, 2 = row has MANE_SELECT, 3 = other);
 (b) independent check of SpliceAI rows accepted as `matched` / `alias_resolved` against the coordinate table;
 (c) SPiP `no_entry_for_row_gene` recount with the allele filter off;
 (d) old SPiP >= 0.20 removed, by tier.
"""
import argparse, json, os, sys
from collections import Counter
import numpy as np, pandas as pd, pyarrow.parquet as pq
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python", "parsing"))
from modules import splicing as sp, gene_identity as gi  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--pq", required=True); ap.add_argument("--row-groups", required=True)
ap.add_argument("--hgnc", required=True); ap.add_argument("--symbol-map", required=True)
ap.add_argument("--curated", required=True, help="resources/gene_transcript_mapping.txt (defines tier 1 only)")
ap.add_argument("--enst-spip-map", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
ident = gi.GeneIdentity(a.symbol_map, a.hgnc, a.enst_spip_map)
cur = {ln.split(',')[2].strip().split('.')[0] for ln in list(open(a.curated))[1:] if ln.count(',') >= 2}  # tier 1 = curated ENST
smap = pd.read_csv(a.symbol_map, sep="\t").set_index("symbol")["ensg"]
need = ["Locus", "Gene", "SYMBOL", "Feature", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SpliceAI", "SPiP", "SPiP_prediction"]
C = Counter()
for rg in (int(x) for x in a.row_groups.split(",")):
    df = pq.ParquetFile(a.pq).read_row_group(rg, columns=need).to_pandas()
    old_spip = df.pop("SPiP_prediction")
    tier = np.where(df.Feature.astype(str).str.split(".").str[0].isin(cur), 1, np.where(df.MANE_SELECT.notna(), 2, 3))
    df = sp.parse_spliceai_custom(df, identity=ident)
    df = sp.parse_spip(df, identity=ident)
    for t, m in zip(tier, df.SPiP_match): C[f"spip_tier{t}:{m}"] += 1
    for t, m in zip(tier, df.spliceai_custom_match): C[f"spliceai_tier{t}:{m}"] += 1
    rem = old_spip.ge(0.2) & df.SPiP_prediction.isna()
    for t in tier[rem.values]: C[f"spip_old_ge0.2_removed_tier{t}"] += 1
    for st in (gi.MATCHED, gi.ALIAS_RESOLVED):
        s = df[df.spliceai_custom_match == st]
        insym = s.spliceai_custom_SYMBOL.isin(smap.index)
        bad = insym & (s.spliceai_custom_SYMBOL.map(smap) != s.Gene.astype(str).str.split(".").str[0])
        C[f"spliceai_{st}_rows"] += len(s); C[f"spliceai_{st}_symbol_not_in_coord_map"] += int((~insym).sum())
        C[f"spliceai_{st}_coord_map_ENSG_differs_from_row_Gene"] += int(bad.sum())
    ne = df[df.SPiP_match == gi.NO_ENTRY].copy()
    C["spip_no_entry_rows"] += len(ne)
    if len(ne):
        ne2 = sp.parse_spip(ne[["SPiP", "Feature", "MANE_SELECT", "MANE_PLUS_CLINICAL"]].assign(Locus=None), identity=ident)
        C["spip_no_entry_would_match_without_allele_filter"] += int((ne2.SPiP_match == gi.MATCHED).sum())
    print("row group", rg, "done", file=sys.stderr, flush=True)
json.dump(dict(sorted(C.items())), open(a.out, "w"), indent=1)
print(json.dumps(dict(sorted(C.items())), indent=1))
