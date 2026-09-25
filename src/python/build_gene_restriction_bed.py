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


def read_mane_gene_coords(path: str) -> dict[str, tuple[str, int, int, str]]:
    """Parses the MANE GTF's `gene` feature lines into {gene_name: (chrom, start, end, strand)}.

    GTF coordinates are 1-based inclusive; BED is 0-based half-open, so start is
    converted (gtf_start - 1) when writing the BED, not here.
    """
    coords = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, _, _, start, end, _, strand = fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6]
            attrs = fields[8]
            gene_name = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name "):
                    gene_name = attr.split('"')[1]
                    break
            if gene_name:
                coords[gene_name] = (chrom, int(start), int(end), strand)
    return coords


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
        help="Optional path to write the list of genes that had no MANE GTF match",
    )
    args = parser.parse_args()

    genes = read_gene_list(args.gene_list)
    print(f"Read {len(genes)} unique genes from {args.gene_list}")

    mane_coords = read_mane_gene_coords(args.mane_gtf)
    print(f"Parsed {len(mane_coords)} gene coordinates from {args.mane_gtf}")

    matched = []
    unmatched = []
    for gene in genes:
        if gene in mane_coords:
            matched.append(gene)
        else:
            unmatched.append(gene)

    print(f"Matched {len(matched)}/{len(genes)} genes ({100 * len(matched) / len(genes):.1f}%)")
    if unmatched:
        print(f"Unmatched ({len(unmatched)}): {', '.join(unmatched[:10])}"
              + (", ..." if len(unmatched) > 10 else ""))

    # Sort by (chrom, start) for a well-formed BED; chrom sort is lexicographic (matches
    # this repo's existing cardio_genes_loc.bed, which is not karyotype-sorted either).
    #
    # Emit BOTH contig-naming conventions per region (bare "11" AND "chr11"): different VCF
    # sources use different conventions, and `bcftools view -R` needs an exact string match with
    # no synonym-file support (unlike VEP's --synonyms). Harmless for whichever convention doesn't
    # match a given input — those rows just never overlap anything.
    rows = []
    for gene in matched:
        chrom, start, end, strand = mane_coords[gene]
        bed_start = max(0, (start - 1) - args.margin)
        bed_end = end + args.margin
        bare_chrom = chrom[3:] if chrom.startswith("chr") else chrom
        chr_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
        for c in {bare_chrom, chr_chrom}:
            rows.append((c, bed_start, bed_end, gene, ".", strand))
    rows.sort(key=lambda r: (r[0], r[1]))

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        for row in rows:
            writer.writerow(row)
    print(f"Wrote {len(rows)} regions to {args.output}")

    if args.unmatched_output:
        with open(args.unmatched_output, "w") as f:
            f.write("\n".join(unmatched) + "\n")
        print(f"Wrote {len(unmatched)} unmatched gene names to {args.unmatched_output}")


if __name__ == "__main__":
    sys.exit(main())
