#!/usr/bin/env python3
"""Can symbol-alias mismatches be told apart from true neighbour-gene attribution?

For rows of a parsed table where the row's own gene (VEP SYMBOL / Gene) is absent from the
SpliceAI entries, resolve each SpliceAI symbol to an Ensembl gene ID WITHOUT any alias list:
by chrom + strand + exon-overlap between SpliceAI's own annotation (spliceai/annotations/
grch38.txt) and the GENCODE gene models (gffutils db built for Pangolin, GENCODE 45).
Then classify each (row gene, SpliceAI symbol) pair:

  same_symbol      SpliceAI symbol == row SYMBOL                 (not a mismatch)
  alias            different symbol, resolves to the row's own ENSG  -> harmless rename
  neighbour        resolves to a different ENSG                  -> true cross-gene attribution
  unresolved       no overlapping GENCODE gene found

Read-only; aggregates only (counts of rows and of distinct symbol pairs).
"""
import argparse
import json
import sqlite3
from collections import Counter, defaultdict

import pyarrow.parquet as pq


def load_spliceai_annotation(path):
    """symbol -> list of (chrom, strand, [(exon_start, exon_end), ...]) ; coords 1-based as in file."""
    ann = defaultdict(list)
    for line in open(path):
        if line.startswith("#"):
            continue
        name, chrom, strand, _s, _e, es, ee = line.rstrip("\n").split("\t")
        ex = list(zip(map(int, es.strip(",").split(",")), map(int, ee.strip(",").split(","))))
        ann[name].append((chrom, strand, ex))
    return ann


def load_gencode_exons(db):
    """(chrom, strand) -> list of (start, end, ensg_unversioned) exons, from the gffutils sqlite."""
    c = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    genes = {}
    for gid, in c.execute("select id from features where featuretype='gene'"):
        genes[gid] = gid.split(".")[0]
    ex = defaultdict(list)
    for seqid, s, e, strand, attrs in c.execute(
            "select seqid,start,end,strand,attributes from features where featuretype='exon'"):
        i = attrs.find('"gene_id":["') + len('"gene_id":["')
        ens = attrs[i:attrs.find('"', i)].split(".")[0]
        ex[(seqid.replace("chr", ""), strand)].append((s, e, ens))
    for k in ex:
        ex[k].sort()
    return ex


def resolve(symbol, ann, gex):
    """Best ENSG for a SpliceAI symbol by summed exon overlap; None if nothing overlaps."""
    score = Counter()
    for chrom, strand, exons in ann.get(symbol, []):
        cand = gex.get((chrom.replace("chr", ""), strand), [])
        lo = exons[0][0]
        hi = exons[-1][1]
        for s, e, ens in cand:
            if e < lo or s > hi:
                continue
            for xs, xe in exons:
                ov = min(e, xe) - max(s, xs) + 1
                if ov > 0:
                    score[ens] += ov
    return score.most_common(1)[0][0] if score else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pq", required=True)
    ap.add_argument("--row-groups", default="all")
    ap.add_argument("--spliceai-annotation", required=True)
    ap.add_argument("--gencode-db", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    f = pq.ParquetFile(a.pq)
    rgs = range(f.num_row_groups) if a.row_groups == "all" else [int(x) for x in a.row_groups.split(",")]
    pairs = Counter()  # (row ENSG, row SYMBOL, spliceai first symbol) -> rows
    for rg in rgs:
        for b in f.iter_batches(batch_size=200_000, row_groups=[rg], columns=["SYMBOL", "Gene", "SpliceAI"]):
            d = b.to_pydict()
            for sym, gene, sai in zip(d["SYMBOL"], d["Gene"], d["SpliceAI"]):
                if not sai or sym is None or gene is None:
                    continue
                syms = [e.split("|")[1] for e in sai.split("&") if e.count("|") >= 9]
                if not syms or sym in syms:
                    continue
                pairs[(gene.split(".")[0], sym, syms[0])] += 1

    ann = load_spliceai_annotation(a.spliceai_annotation)
    gex = load_gencode_exons(a.gencode_db)
    cache = {}
    rows, npairs = Counter(), Counter()
    for (ens, sym, ssym), n in pairs.items():
        if ssym not in cache:
            cache[ssym] = resolve(ssym, ann, gex)
        r = cache[ssym]
        cls = "unresolved" if r is None else ("alias" if r == ens else "neighbour")
        rows[cls] += n
        npairs[cls] += 1
    json.dump({"row_groups": list(rgs), "rows": dict(rows), "distinct_pairs": dict(npairs),
               "distinct_spliceai_symbols": len(cache),
               "unresolved_symbols_in_spliceai_annotation": sum(1 for s, r in cache.items() if r is None and s in ann),
               "symbols_missing_from_spliceai_annotation": sum(1 for s in cache if s not in ann)},
              open(a.out, "w"), indent=2)
    print(open(a.out).read())


if __name__ == "__main__":
    main()
