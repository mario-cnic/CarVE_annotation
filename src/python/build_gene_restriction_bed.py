#!/usr/bin/env python3
"""
Builds resources/v5_genes_loc.bed: gene-body coordinates for every gene in
Complete_gene_list_V5, for the gene-restricted tier's (SpliceAI/Pangolin/SPiP)
region-restriction filter (`bcftools view -R`).

Purely local/offline: no cluster, no network. Joins Complete_gene_list_V5's gene
symbols against the local MANE GTF's `gene` feature lines by gene_name, pads each
interval by a splice-prediction window margin, and logs (does not fail on) any gene
that doesn't match by symbol.

Name column (col 4) is a bare gene symbol, e.g. "MYBPC3" — NOT the same format as
resources/cardio_genes_loc.bed's "GENE_NM_ENST" (that file's name column matters to
main.sh:gene_coords.sh's own "already have this gene" grep check; this file's does
not, since it only needs to answer "which genomic interval," not "which transcript" —
no transcript resolution happens here by design).
"""

import argparse
import csv
import gzip
import sys
from collections import OrderedDict


def read_gene_list(path: str) -> list[str]:
    """Reads Complete_gene_list_V5.csv, returns deduplicated gene symbols in first-seen order."""
    genes = OrderedDict()
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gene = row["GENE"].strip()
            if gene:
                genes[gene] = True
    return list(genes.keys())


def _attr(attrs: str, key: str):
    for attr in attrs.split(";"):
        attr = attr.strip()
        if attr.startswith(key + " "):
            return attr.split('"')[1]
    return None


def read_gtf_genes(path: str) -> dict[str, tuple[str, int, int, str, str | None]]:
    """{gene_name: (chrom, start, end, strand, gene_id without version)} from a GTF's `gene` lines.
    A name seen on a primary contig keeps that record over one on an alt/fix/random contig."""
    genes = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            name = _attr(fields[8], "gene_name")
            if not name:
                continue
            gid = _attr(fields[8], "gene_id")
            rec = (fields[0], int(fields[3]), int(fields[4]), fields[6], gid.split(".")[0] if gid else None)
            if name not in genes or (is_primary_contig(rec[0]) and not is_primary_contig(genes[name][0])):
                genes[name] = rec
    return genes


def is_primary_contig(chrom: str) -> bool:
    """chr1-22, X, Y, M (with or without `chr`); not alt / fix / random / chrUn / EBV."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return c in {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}


def resolve_genes(genes: list[str], mane: dict, primary: dict | None = None) -> tuple[dict, dict]:
    """{listed gene: [(chrom, start, end, strand, how), ...]} and {listed gene: reason} for the rest (MISC-4).

    A listed name is looked up in MANE by exact symbol, then case-insensitively (`C10ORF71` ->
    `C10orf71`, only if one MANE symbol fits), then as an Ensembl gene ID (`ENSG...`, version ignored).
    A MANE record is always kept. When it sits only on an alt/fix contig, or the name is not in MANE,
    the primary-assembly GTF (`primary`, e.g. GENCODE) is tried the same way and its record is added
    if it is on a primary contig. Genes left on alt/fix contigs only are reported (`mane_alt_contig_only`).
    """
    def lookup(name, table):
        if not table:
            return None, None
        if name in table:
            return table[name], "symbol"
        folded = [k for k in table if k.upper() == name.upper()]
        if len(folded) == 1:
            return table[folded[0]], "symbol_case"
        if name.upper().startswith("ENSG"):
            hits = [v for v in table.values() if v[4] == name.split(".")[0].upper()]
            if len(hits) == 1:
                return hits[0], "ensembl_id"
        return None, None

    found, missing = {}, {}
    for name in genes:
        recs = []
        rec, how = lookup(name, mane)
        if rec:
            recs.append((*rec[:4], f"mane_{how}"))
        if not rec or not is_primary_contig(rec[0]):
            prec, phow = lookup(name, primary)
            if prec and is_primary_contig(prec[0]):
                recs.append((*prec[:4], f"primary_gtf_{phow}"))
        if recs:
            found[name] = recs
        if not any(is_primary_contig(r[0]) for r in recs):
            missing[name] = "mane_alt_contig_only" if recs else "no_match"
    return found, missing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gene-list",
        default="/home/mruizp/data_lab_PGP/shared/utils/data/Complete_gene_list_V5.csv",
        help="Complete_gene_list_V5.csv path",
    )
    parser.add_argument(
        "--mane-gtf",
        default="/home/mruizp/data_lab_PGP/shared/utils/data/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz",
        help="MANE GRCh38 GTF path (gzipped or plain)",
    )
    parser.add_argument(
        "--primary-gtf",
        default=None,
        help="Optional primary-assembly GTF (e.g. GENCODE 45) for genes MANE lacks or has on alt/fix contigs only",
    )
    parser.add_argument("--output", default="resources/v5_genes_loc.bed", help="Output BED path")
    parser.add_argument(
        "--margin",
        type=int,
        default=10000,
        help="Padding added to each side of the gene body (matches SpliceAI/Pangolin's -d window)",
    )
    parser.add_argument(
        "--unmatched-output",
        default=None,
        help="Optional path to write the genes without a primary-contig region (gene<TAB>reason)",
    )
    args = parser.parse_args()

    genes = read_gene_list(args.gene_list)
    print(f"Read {len(genes)} unique genes from {args.gene_list}")

    mane = read_gtf_genes(args.mane_gtf)
    print(f"Parsed {len(mane)} gene coordinates from {args.mane_gtf}")
    primary = read_gtf_genes(args.primary_gtf) if args.primary_gtf else None
    if primary is not None:
        print(f"Parsed {len(primary)} gene coordinates from {args.primary_gtf}")

    found, unmatched = resolve_genes(genes, mane, primary)
    how = {}
    for recs in found.values():
        for r in recs:
            how[r[4]] = how.get(r[4], 0) + 1
    print(f"Genes with a primary-contig region: {len(genes) - len(unmatched)}/{len(genes)}; records by match: {how}")
    if unmatched:
        print(f"Without a primary-contig region ({len(unmatched)}): "
              + ", ".join(f"{g} ({r})" for g, r in list(unmatched.items())[:10]) + (", ..." if len(unmatched) > 10 else ""))

    # Sort by (chrom, start) for a well-formed BED; chrom sort is lexicographic (matches
    # this repo's existing cardio_genes_loc.bed, which is not karyotype-sorted either).
    #
    # Emit BOTH contig-naming conventions per region (bare "11" AND "chr11"): different VCF
    # sources use different conventions, and `bcftools view -R` needs an exact string match with
    # no synonym-file support (unlike VEP's --synonyms). Harmless for whichever convention doesn't
    # match a given input — those rows just never overlap anything.
    rows = set()
    for gene, recs in found.items():
        for chrom, start, end, strand, _ in recs:
            bed_start = max(0, (start - 1) - args.margin)
            bed_end = end + args.margin
            bare_chrom = chrom[3:] if chrom.startswith("chr") else chrom
            chr_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
            for c in {bare_chrom, chr_chrom}:
                rows.add((c, bed_start, bed_end, gene, ".", strand))
    rows = sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3]))

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        for row in rows:
            writer.writerow(row)
    print(f"Wrote {len(rows)} regions to {args.output}")

    if args.unmatched_output:
        with open(args.unmatched_output, "w") as f:
            f.write("".join(f"{g}\t{r}\n" for g, r in unmatched.items()))
        print(f"Wrote {len(unmatched)} genes without a primary-contig region to {args.unmatched_output}")


if __name__ == "__main__":
    sys.exit(main())
