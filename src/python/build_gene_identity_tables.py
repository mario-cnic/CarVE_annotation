#!/usr/bin/env python3
"""Build resources/spliceai_symbol_to_ensg.grch38.gencode45.tsv + a build manifest.

Maps each gene symbol in SpliceAI's bundled GRCh38 annotation (the `grch38` annotation that
src/python/annotate_spliceai.py passes to spliceai.utils.Annotator) to the Ensembl gene ID whose
GENCODE 45 exons overlap it most on the same chromosome and strand. No alias list is used.

Input gene models: the gffutils db the Pangolin step reads (GENCODE 45 = the gene models of the VEP
111 cache), built with the Ensembl_canonical filter (canonical-transcript exons only); symbols the
coordinate route cannot place are left to the HGNC route at run time.

Coordinates: exon intervals overlap-tested in the same 1-based inclusive convention both files use;
a +-1 offset between the two sources does not change the best overlap.
"""
import argparse
import hashlib
import json
import sqlite3
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_spliceai_annotation(path):
    ann = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            name, chrom, strand, _s, _e, es, ee = line.rstrip("\n").split("\t")
            exons = list(zip(map(int, es.strip(",").split(",")), map(int, ee.strip(",").split(","))))
            ann[name].append((chrom.replace("chr", ""), strand, exons))
    return ann


def load_gencode_exons(db_path):
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    ex = defaultdict(list)
    for seqid, s, e, strand, attrs in con.execute(
            "select seqid,start,end,strand,attributes from features where featuretype='exon'"):
        i = attrs.find('"gene_id":["') + len('"gene_id":["')
        ens = attrs[i:attrs.find('"', i)].split(".")[0]
        ex[(seqid.replace("chr", ""), strand)].append((s, e, ens))
    for k in ex:
        ex[k].sort()
    return ex


def best_gene(symbol, ann, gex):
    score = Counter()
    for chrom, strand, exons in ann[symbol]:
        lo, hi = exons[0][0], exons[-1][1]
        for s, e, ens in gex.get((chrom, strand), []):
            if e < lo or s > hi:
                continue
            for xs, xe in exons:
                ov = min(e, xe) - max(s, xs) + 1
                if ov > 0:
                    score[ens] += ov
    top = score.most_common(2)
    return top


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--spliceai-annotation", required=True)
    ap.add_argument("--gencode-db", required=True)
    ap.add_argument("--out-tsv", required=True)
    a = ap.parse_args()

    ann = load_spliceai_annotation(a.spliceai_annotation)
    gex = load_gencode_exons(a.gencode_db)
    n_unplaced = 0
    with open(a.out_tsv, "w") as out:
        out.write("symbol\tensg\toverlap_bp\trunner_up_ensg\trunner_up_bp\n")
        for sym in sorted(ann):
            top = best_gene(sym, ann, gex)
            if not top:
                n_unplaced += 1
                continue
            ru = top[1] if len(top) > 1 else ("", 0)
            out.write(f"{sym}\t{top[0][0]}\t{top[0][1]}\t{ru[0]}\t{ru[1]}\n")
    manifest = {
        "built": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "spliceai_annotation": a.spliceai_annotation, "spliceai_annotation_sha256": sha256(a.spliceai_annotation),
        "gencode_db": a.gencode_db, "gencode_db_sha256": sha256(a.gencode_db),
        "symbols_in_spliceai_annotation": len(ann), "symbols_not_placed_by_coordinates": n_unplaced,
        "output": a.out_tsv, "output_sha256": sha256(a.out_tsv),
        "argv": sys.argv,
    }
    with open(a.out_tsv + ".manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
