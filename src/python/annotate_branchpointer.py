#!/usr/bin/env python3
"""
LaBranchoR branchpoint lookup (GRCh38): a positional annotation, not an allele-aware prediction.

The LaBranchoR BED holds the top predicted branchpoint per 3' splice site (GENCODE v19 introns, lifted to
GRCh38; built by src/tools/build_labranchor_grch38.py). A variant within 2 bp of such a branchpoint gets:
- LaBranchoR_score: LaBranchoR's probability for that branchpoint (0.0 - 1.0), a property of the site, not
  of the variant's alleles.
- LaBranchoR_acc_dist: distance in bp from the variant to that intron's 3' splice site.

No REF/ALT sequence is scored (PRED-18): use these fields to see that a variant sits on a predicted
branchpoint, not as evidence that it disrupts it.
"""

import argparse
import os

import pysam

DEFAULT_LABRANCHOR_BED = "/home/mruizp/data_lab_PGP/resources/annotation/labranchor/labranchor_grch38_top.bed.gz"


def nearby_branchpoint(tbx, chrom: str, pos: int):
    """(score, distance to the 3'SS) of the first BED branchpoint within 2 bp of `pos`, or None.

    BED columns: chrom, bp_start, bp_end, acceptor, score, strand. A row qualifies when the variant lies in
    the branchpoint interval, or 18-44 bp from the acceptor (LaBranchoR's usual branchpoint range)."""
    query_chrom = chrom if chrom.startswith("chr") else "chr" + chrom
    try:
        rows = tbx.fetch(query_chrom, max(0, pos - 2), pos + 2)
    except ValueError:  # contig not in the BED index (e.g. alt contigs, chrM)
        return None
    for row in rows:
        parts = row.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        bp_start, bp_end, acc_pos, score = int(parts[1]), int(parts[2]), int(parts[3]), float(parts[4])
        dist = abs(acc_pos - pos)
        if bp_start <= pos <= bp_end or 18 <= dist <= 44:
            return score, dist
    return None


def annotate_branchpoint_vcf(in_vcf: str, out_vcf: str, labranchor_bed: str = DEFAULT_LABRANCHOR_BED):
    if not os.path.exists(labranchor_bed):
        raise FileNotFoundError(f"LaBranchoR BED not found: {labranchor_bed}")
    tbx = pysam.TabixFile(labranchor_bed)
    in_v = pysam.VariantFile(in_vcf)
    header = in_v.header
    if "LaBranchoR_score" not in header.info:
        header.info.add("LaBranchoR_score", "1", "Float",
                        "LaBranchoR probability of the predicted branchpoint within 2 bp of the variant (site property, alleles not scored)")
    if "LaBranchoR_acc_dist" not in header.info:
        header.info.add("LaBranchoR_acc_dist", "1", "Integer", "Distance in bp from the variant to that branchpoint's 3' splice site")
    out_v = pysam.VariantFile(out_vcf, "w", header=header)

    total = annotated = 0
    for rec in in_v:
        total += 1
        hit = nearby_branchpoint(tbx, rec.chrom, rec.pos)
        if hit:
            rec.info["LaBranchoR_score"] = round(hit[0], 4)
            rec.info["LaBranchoR_acc_dist"] = hit[1]
            annotated += 1
        out_v.write(rec)
    in_v.close()
    out_v.close()
    tbx.close()
    print(f"LaBranchoR lookup completed. Total: {total}, near a predicted branchpoint: {annotated}")


def main():
    parser = argparse.ArgumentParser(description="LaBranchoR branchpoint lookup (positional, alleles not scored)")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_vcf", help="Output VCF file")
    parser.add_argument("--labranchor-bed", default=DEFAULT_LABRANCHOR_BED, help="Path to tabix-indexed GRCh38 LaBranchoR top branchpoint BED file")
    args = parser.parse_args()
    annotate_branchpoint_vcf(args.input_vcf, args.output_vcf, labranchor_bed=args.labranchor_bed)


if __name__ == "__main__":
    main()
