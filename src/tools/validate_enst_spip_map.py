#!/usr/bin/env python3
"""Independent checks of resources/enst_to_spip_nm.*.tsv (Ensembl transcript -> SPiP-database RefSeq transcript by exon structure).

 (1) Agreement with known pairs: every MANE (ENST, NM) pair whose NM accession exists in SPiP's database. Does the
     structural map give that ENST the same NM, a different NM only, or nothing? Only "different NM only" is a mis-attribution.
 (2) Cross-gene pairs: ENST's GENCODE gene versus the SPiP transcript's gene (symbol resolved through HGNC).
Read-only; aggregates plus the list of cross-gene pairs.
"""
import argparse, csv, gzip, json, re
from collections import Counter, defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--map", required=True); ap.add_argument("--mane-gtf", required=True)
ap.add_argument("--spip-db-tsv", required=True); ap.add_argument("--hgnc", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()

spip_nm = {r["nm"] for r in csv.DictReader(open(a.spip_db_tsv), delimiter="\t")}
pairs = defaultdict(set); rows = []
for r in csv.DictReader(open(a.map), delimiter="\t"):
    pairs[r["enst"]].add(r["spip_nm"]); rows.append(r)

mane = {}
for ln in gzip.open(a.mane_gtf, "rt"):
    if "\ttranscript\t" in ln:
        t = re.search(r'transcript_id "([^"]+)"', ln).group(1).split(".")[0]
        x = re.search(r'db_xref "RefSeq:([^"]+)"', ln)
        if x: mane[t] = x.group(1).split(".")[0]
c1 = Counter(); diff = []
for enst, nm in mane.items():
    if nm not in spip_nm: c1["MANE NM not in SPiP database (cannot test)"] += 1; continue
    got = pairs.get(enst, set())
    if nm in got: c1["same NM"] += 1
    elif got: c1["different NM only (mis-attribution?)"] += 1; diff.append((enst, nm, sorted(got)))
    else: c1["no structural counterpart"] += 1

cur, other = {}, defaultdict(set)
for r in csv.DictReader(open(a.hgnc), delimiter="\t"):
    e = r["ensembl_gene_id"]
    if not e: continue
    cur.setdefault(r["symbol"], set()).add(e)
    for k in ("prev_symbol", "alias_symbol"):
        for s in r[k].split("|"):
            if s: other[s].add(e)
c2 = Counter(); cross = []
for r in rows:
    g = r["gene_id"].split(".")[0]
    if r["gene_name"] == r["spip_symbol"]: c2["same symbol"] += 1
    elif g in (cur.get(r["spip_symbol"]) or other.get(r["spip_symbol"]) or set()): c2["same gene via HGNC (renamed)"] += 1
    elif not (cur.get(r["spip_symbol"]) or other.get(r["spip_symbol"])): c2["SPiP symbol unknown to HGNC (cannot test)"] += 1
    else: c2["DIFFERENT gene"] += 1; cross.append((r["enst"], r["gene_name"], r["spip_nm"], r["spip_symbol"], r["match_type"]))
res = {"known_pairs": dict(c1), "gene_check": dict(c2), "different_nm_only_examples": diff[:15], "cross_gene_pairs": cross}
json.dump(res, open(a.out, "w"), indent=1)
print(json.dumps({"known_pairs": res["known_pairs"], "gene_check": res["gene_check"]}, indent=1)); print("cross-gene examples:", cross[:12]); print("different-NM-only examples:", diff[:6])
