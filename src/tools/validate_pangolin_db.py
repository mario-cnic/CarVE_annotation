#!/usr/bin/env python3
"""Validates a Pangolin gffutils annotation db against the GTF it was built from.

Checks: species/release directives, human contig set, gene ID namespace, one Ensembl_canonical
transcript per gene, gene/transcript/exon counts equal to an independent parse of the GTF, and
spot genes whose coordinates, strand and canonical exons must match the GTF exactly. Spot
positions are then queried through Pangolin's own `get_genes()` (extracted from pangolin.py by
AST, so torch is not imported) to confirm each gene is found with a non-empty exon list.

Needs gffutils (same version as the reader, see BUG_TRACKER.md MISC-12). Exit 1 on any failure.
"""

import argparse
import ast
import collections
import gzip
import re
import sys

import gffutils

ATTR = re.compile(r'(\S+) "([^"]*)"')
EXPECTED_CONTIGS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
SPOT_GENES = ["MYBPC3", "TTN", "MYH7", "LMNA", "SCN5A", "EMD", "DMD", "SRY"]


def parse_gtf(path, tag="Ensembl_canonical"):
    """Independent GTF parse mirroring create_db.py's default filter.

    Returns (counts, genes, exons): counts by featuretype after filtering; genes maps gene_name to
    (gene_id, chrom, start, end, strand); exons maps gene_id to a sorted list of (start, end) of the
    tagged transcript's exons.
    """
    counts = collections.Counter()
    genes, exons = {}, collections.defaultdict(list)
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            ft = p[2]
            if ft not in ("gene", "transcript", "exon"):
                continue
            a = collections.defaultdict(list)
            for k, v in ATTR.findall(p[8]):
                a[k].append(v)
            if ft != "gene" and tag not in a.get("tag", []):
                continue
            counts[ft] += 1
            gid = a["gene_id"][0]
            if ft == "gene":
                genes.setdefault(a["gene_name"][0], (gid, p[0], int(p[3]), int(p[4]), p[6]))
            elif ft == "exon":
                exons[gid].append((int(p[3]), int(p[4])))
    return counts, genes, {k: sorted(v) for k, v in exons.items()}


def load_get_genes(pangolin_py):
    """Compiles Pangolin's get_genes() from source without executing the module's imports."""
    tree = ast.parse(open(pangolin_py).read())
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "get_genes")
    ns = {}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), pangolin_py, "exec"), ns)
    return ns["get_genes"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True)
    ap.add_argument("--gtf", required=True, help="GTF the db was built from (.gtf.gz)")
    ap.add_argument("--pangolin-py", required=True, help="Path to pangolin/pangolin.py")
    a = ap.parse_args()

    fails = []

    def check(ok, msg):
        print(f"[{'PASS' if ok else 'FAIL'}] {msg}")
        if not ok:
            fails.append(msg)

    db = gffutils.FeatureDB(a.db)
    con = db.conn
    directives = [r[0] for r in con.execute("select directive from directives")]
    print("directives:", directives)
    check(any("human genome (GRCh38)" in d for d in directives), "directives declare human GRCh38")

    ft_counts = dict(con.execute("select featuretype, count(*) from features group by featuretype"))
    gtf_counts, gtf_genes, gtf_exons = parse_gtf(a.gtf)
    print("db featuretype counts:", ft_counts, "| GTF (filtered):", dict(gtf_counts))
    for ft in ("gene", "transcript", "exon"):
        check(ft_counts.get(ft, 0) == gtf_counts[ft], f"{ft} count db={ft_counts.get(ft, 0)} GTF={gtf_counts[ft]}")

    contigs = {r[0] for r in con.execute("select distinct seqid from features")}
    print("contigs:", sorted(contigs))
    check(EXPECTED_CONTIGS <= contigs, "chr1-chr22, chrX, chrY all present")
    check(contigs <= EXPECTED_CONTIGS | {"chrM"}, "no contigs beyond chr1-22/X/Y/M")

    n_ensg = con.execute("select count(*) from features where featuretype='gene' and id like 'ENSG%'").fetchone()[0]
    n_mouse = con.execute("select count(*) from features where attributes like '%ENSMUS%'").fetchone()[0]
    check(n_ensg == ft_counts.get("gene"), f"all {n_ensg} gene IDs are ENSG")
    check(n_mouse == 0, f"no ENSMUS identifiers ({n_mouse})")

    multi = con.execute(
        "select count(*) from (select parent from relations r join features f on f.id=r.child "
        "where f.featuretype='transcript' and r.level=1 group by parent having count(*)>1)"
    ).fetchone()[0]
    check(multi == 0, f"one transcript per gene ({multi} genes with >1)")
    not_canon = sum(1 for t in db.features_of_type("transcript") if "Ensembl_canonical" not in t.attributes.get("tag", []))
    check(not_canon == 0, f"every transcript tagged Ensembl_canonical ({not_canon} not)")

    get_genes = load_get_genes(a.pangolin_py)
    for name in SPOT_GENES:
        gid, chrom, start, end, strand = gtf_genes[name]
        g = db[gid]
        check((g.seqid, g.start, g.end, g.strand) == (chrom, start, end, strand),
              f"{name} {gid} {chrom}:{start}-{end}({strand}) matches GTF")
        db_exons = sorted((e.start, e.end) for e in db.children(g, featuretype="exon"))
        check(db_exons == gtf_exons[gid], f"{name} canonical exons db={len(db_exons)} GTF={len(gtf_exons[gid])}")
        ex_mid = (gtf_exons[gid][len(gtf_exons[gid]) // 2][0] + gtf_exons[gid][len(gtf_exons[gid]) // 2][1]) // 2
        for chrom_q in (chrom, chrom[3:]):
            pos_g, neg_g = get_genes(chrom_q, ex_mid, db)
            hit = (pos_g if strand == "+" else neg_g).get(gid)
            check(bool(hit) and len(hit) == 2 * len(gtf_exons[gid]),
                  f"get_genes({chrom_q}:{ex_mid}) finds {name} on {strand} with {len(hit or [])} exon bounds")

    pos_g, neg_g = get_genes("chr1", 5000, db)
    check(not pos_g and not neg_g, "get_genes(chr1:5000) intergenic returns no genes")

    print(f"\n{len(fails)} failure(s)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
