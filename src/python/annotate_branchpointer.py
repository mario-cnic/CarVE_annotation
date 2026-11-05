#!/usr/bin/env python3
"""
LaBranchoR branchpoint lookup (GRCh38): a positional annotation, not an allele-aware prediction.

The LaBranchoR BED holds the top predicted branchpoint per 3' splice site (GENCODE v19 introns, lifted to
GRCh38; built by src/tools/build_labranchor_grch38.py). A variant within 2 bp of such a branchpoint gets:
- LaBranchoR_score: LaBranchoR's probability for that branchpoint (0.0 - 1.0), a property of the site, not
  of the variant's alleles.
- LaBranchoR_acc_dist: distance in bp from the variant to that intron's 3' splice site.

These two fields score no REF/ALT sequence (PRED-18): they show that a variant sits on a predicted
branchpoint, not that it disrupts it. The allele-aware fields come from LaBranchoR's in silico mutagenesis
table (--labranchor-ism; built by src/tools/build_labranchor_ism_grch38.py): the model scored on REF and ALT
for every SNV in the 70 bases upstream of a GENCODE v19 3' splice site. For an SNV in several such windows
the 3'SS with the lowest d_best_bp is reported:
- LaBranchoR_ISM_d_best: change of the reference top branchpoint's probability (P_alt - P_ref max).
- LaBranchoR_ISM_d_max: change of the best branchpoint probability anywhere in the window.
- LaBranchoR_ISM_dist: bases from the SNV to that 3' splice site (1-70, upstream).
- LaBranchoR_ISM_n_3pss: number of 3' splice sites whose window holds the SNV.
Informational: no calibrated threshold, nothing here changes NEW_IMPACT.
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


ISM_INFO = [
    ("LaBranchoR_ISM_d_best", "Float", "LaBranchoR ISM: change of the reference top branchpoint probability (ALT - REF), worst 3'SS"),
    ("LaBranchoR_ISM_d_max", "Float", "LaBranchoR ISM: change of the best branchpoint probability in the window (ALT - REF), same 3'SS"),
    ("LaBranchoR_ISM_dist", "Integer", "LaBranchoR ISM: bases from the SNV to that 3' splice site (1-70, upstream)"),
    ("LaBranchoR_ISM_n_3pss", "Integer", "LaBranchoR ISM: number of 3' splice sites whose 70-base window holds the SNV"),
]


def _fetch(tbx, chrom: str, start0: int, end: int):
    for name in (chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom):
        if name in tbx.contigs:
            return tbx.fetch(name, start0, end)
    return []


def ism_effect(tbx, chrom: str, pos: int, ref: str, alt: str):
    """(d_best, d_max, dist, n_3pss) of a biallelic SNV from the ISM table, or None if no row matches.

    ISM columns: chrom, pos (1-based), ref, alt (+ strand), strand, 3pSS, d_best_bp, d_max_bp. The 3pSS
    coordinate lies after the window on '+' (pos in [3pSS-69, 3pSS]) and before it on '-'
    (pos in [3pSS+1, 3pSS+70]), so the distance is 3pSS + 1 - pos or pos - 3pSS."""
    if len(ref) != 1 or len(alt) != 1:
        return None
    hits = []
    for row in _fetch(tbx, chrom, pos - 1, pos):
        f = row.rstrip("\n").split("\t")
        if int(f[1]) != pos or f[2].upper() != ref.upper() or f[3].upper() != alt.upper():
            continue
        pss = int(f[5])
        dist = pss + 1 - pos if f[4] == "+" else pos - pss
        hits.append((float(f[6]), float(f[7]), dist))
    if not hits:
        return None
    d_best, d_max, dist = min(hits)
    return d_best, d_max, dist, len(hits)


def annotate_branchpoint_vcf(in_vcf: str, out_vcf: str, labranchor_bed: str = DEFAULT_LABRANCHOR_BED,
                             labranchor_ism: str | None = None):
    if not os.path.exists(labranchor_bed):
        raise FileNotFoundError(f"LaBranchoR BED not found: {labranchor_bed}")
    if labranchor_ism and not os.path.exists(labranchor_ism):
        raise FileNotFoundError(f"LaBranchoR ISM table not found: {labranchor_ism}")
    tbx = pysam.TabixFile(labranchor_bed)
    tbx_ism = pysam.TabixFile(labranchor_ism) if labranchor_ism else None
    in_v = pysam.VariantFile(in_vcf)
    header = in_v.header
    if "LaBranchoR_score" not in header.info:
        header.info.add("LaBranchoR_score", "1", "Float",
                        "LaBranchoR probability of the predicted branchpoint within 2 bp of the variant (site property, alleles not scored)")
    if "LaBranchoR_acc_dist" not in header.info:
        header.info.add("LaBranchoR_acc_dist", "1", "Integer", "Distance in bp from the variant to that branchpoint's 3' splice site")
    for tag, typ, desc in ISM_INFO:   # declared even without the table: the merge transfers these tags
        if tag not in header.info:
            header.info.add(tag, "1", typ, desc)
    out_v = pysam.VariantFile(out_vcf, "w", header=header)

    total = annotated = ism_scored = 0
    for rec in in_v:
        total += 1
        hit = nearby_branchpoint(tbx, rec.chrom, rec.pos)
        if hit:
            rec.info["LaBranchoR_score"] = round(hit[0], 4)
            rec.info["LaBranchoR_acc_dist"] = hit[1]
            annotated += 1
        if tbx_ism is not None and rec.alts and len(rec.alts) == 1:
            eff = ism_effect(tbx_ism, rec.chrom, rec.pos, rec.ref, rec.alts[0])
            if eff:
                for (tag, _, _), value in zip(ISM_INFO, (round(eff[0], 6), round(eff[1], 6), eff[2], eff[3])):
                    rec.info[tag] = value
                ism_scored += 1
        out_v.write(rec)
    in_v.close()
    out_v.close()
    tbx.close()
    if tbx_ism is not None:
        tbx_ism.close()
    print(f"LaBranchoR lookup completed. Total: {total}, near a predicted branchpoint: {annotated}, ISM-scored SNVs: {ism_scored}")


def main():
    parser = argparse.ArgumentParser(description="LaBranchoR branchpoint lookup (positional, alleles not scored)")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_vcf", help="Output VCF file")
    parser.add_argument("--labranchor-bed", default=DEFAULT_LABRANCHOR_BED, help="Path to tabix-indexed GRCh38 LaBranchoR top branchpoint BED file")
    parser.add_argument("--labranchor-ism", default=None, help="Path to tabix-indexed GRCh38 LaBranchoR ISM table (allele-aware SNV scores)")
    args = parser.parse_args()
    annotate_branchpoint_vcf(args.input_vcf, args.output_vcf, labranchor_bed=args.labranchor_bed,
                             labranchor_ism=args.labranchor_ism)


if __name__ == "__main__":
    main()
