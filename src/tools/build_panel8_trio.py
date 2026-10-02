#!/usr/bin/env python3
"""Builds the synthetic test panel `panel8_trio`: the 72 variants of `panel7_refcorrected.vcf.gz`
plus 6 variants in an 8th gene (EMD, chrX), with a synthetic trio (PROBAND_P8 = daughter,
FATHER_P8, MOTHER_P8). No patient data. Written for testing the genotype outputs.

Outputs (deterministic for a given seed):
  <out-prefix>.vcf.gz (+ .tbi)   FORMAT GT:GQ:DP:AD, contig names as in panel7 (1..22, X, bare)
  <out-prefix>.expected.tsv      Locus, sample_id, gt_raw, gt_status expected from the scenario that
                                 generated each genotype (not from the parser under test)

Check a genotype table against the expectation:
  build_panel8_trio.py --compare-genotypes <run>.genotypes.pq --expected <out-prefix>.expected.tsv

Scenarios (record -> PROBAND / FATHER / MOTHER): inherited from mother or father, de novo,
homozygous alternative with carrier parents, parent-only variants, all hom-ref, compound-het pair
(two DSP variants, one from each parent), mother no-call, proband partial no-call, phased 1|0,
unphased 1/0 with no reference reads (split-multiallelic style), missing GQ/DP, and chrX variants
with a haploid father (hemizygous alt, ref, no-call).
"""

import argparse
import os
import random
import sys

import pysam

PANEL7 = "test_data/raw_vcfs/panel7_refcorrected.vcf.gz"
FASTA = "/home/mruizp/data_references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta"
SAMPLES = ("PROBAND_P8", "FATHER_P8", "MOTHER_P8")
EMD_POSITIONS = [154379400, 154379900, 154380300, 154380800, 154381100, 154381400]  # inside the gene body
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}

# scenario -> (proband, father, mother) GT strings
BULK = {
    "from_mother": ("0/1", "0/0", "0/1"),
    "from_father": ("0/1", "0/1", "0/0"),
    "de_novo": ("0/1", "0/0", "0/0"),
    "hom_alt": ("1/1", "0/1", "0/1"),
    "father_only": ("0/0", "0/1", "0/0"),
    "mother_only": ("0/0", "0/0", "0/1"),
    "all_hom_ref": ("0/0", "0/0", "0/0"),
}
BULK_WEIGHTS = [("from_mother", 40), ("from_father", 22), ("de_novo", 6), ("hom_alt", 8),
                ("father_only", 10), ("mother_only", 10), ("all_hom_ref", 4)]


def status_of(gt):
    """Expected gt_status for a GT string, from the scenario table, written independently of the parser."""
    table = {"0/0": "hom_ref", "0/1": "het", "1/1": "hom_alt", "./.": "no_call", ".": "no_call", "./1": "partial_no_call",
             "1|0": "het", "1/0": "alt_plus_other_allele", "1": "haploid_alt", "0": "haploid_ref"}
    return table[gt]


def format_call(gt, rng, gq=None, dp=None):
    """GT:GQ:DP:AD string with read counts consistent with the genotype."""
    if gt in ("./.", "."):
        return f"{gt}:.:.:."
    depth = rng.randint(18, 60) if dp is None else dp
    alt_frac = {"0/0": 0.0, "0/1": 0.5, "1|0": 0.5, "1/1": 1.0, "./1": 0.5, "1": 1.0, "0": 0.0}.get(gt)
    if gt == "1/0":   # alt plus another allele: no reference reads, some reads for the other allele
        alt = depth - 6
        return f"{gt}:{gq or 99}:{depth}:0,{alt}"
    alt = round(depth * alt_frac) if gt not in ("0/1", "1|0", "./1") else max(1, min(depth - 1, depth // 2 + rng.randint(-4, 4)))
    ad = f"{depth - alt},{alt}"
    return f"{gt}:{99 if gq is None else gq}:{depth}:{ad}"


def emd_variants(fasta):
    """Six EMD (chrX) records with REF taken from the FASTA: 4 SNVs, 1 deletion, 1 insertion."""
    out = []
    for i, pos in enumerate(EMD_POSITIONS):
        seq = fasta.fetch("chrX", pos - 1, pos + 2).upper()
        if i < 4:
            out.append(("X", pos, seq[0], TRANSITION[seq[0]]))
        elif i == 4:
            out.append(("X", pos, seq[:2], seq[0]))          # deletion of the base after the anchor
        else:
            out.append(("X", pos, seq[0], seq[0] + ("A" if seq[0] != "A" else "C")))  # insertion
    return out


def build(panel7, fasta_path, out_prefix, seed):
    rng = random.Random(seed)
    fasta = pysam.FastaFile(fasta_path)
    src = pysam.VariantFile(panel7)
    header = src.header.copy()
    header.add_meta("FORMAT", items=[("ID", "GT"), ("Number", "1"), ("Type", "String"), ("Description", "Genotype")])
    header.add_meta("FORMAT", items=[("ID", "GQ"), ("Number", "1"), ("Type", "Integer"), ("Description", "Genotype Quality")])
    header.add_meta("FORMAT", items=[("ID", "DP"), ("Number", "1"), ("Type", "Integer"), ("Description", "Read depth")])
    header.add_meta("FORMAT", items=[("ID", "AD"), ("Number", "R"), ("Type", "Integer"), ("Description", "Allelic depths")])
    lengths = dict(zip(fasta.references, fasta.lengths))
    for bare in ("X",):
        header.contigs.add(bare, length=lengths["chr" + bare])
        header.contigs.add("chr" + bare, length=lengths["chr" + bare])
    for s in SAMPLES:
        header.add_sample(s)

    records = []   # (chrom, pos, ref, alt, info/qual/id source, gts)
    for r in src:
        records.append(["panel7", r, None])
    for chrom, pos, ref, alt in emd_variants(fasta):
        records.append(["emd", (chrom, pos, ref, alt), None])

    # fixed scenarios by position in the list
    dsp = [i for i, (kind, r, _) in enumerate(records) if kind == "panel7" and r.info.get("gene") == "DSP"]
    panel7_idx = [i for i, (kind, _, _) in enumerate(records) if kind == "panel7"]
    special = {
        dsp[0]: ("0/1", "0/1", "0/0"),            # compound-het pair, variant 1 from the father
        dsp[1]: ("0/1", "0/0", "0/1"),            # compound-het pair, variant 2 from the mother
        panel7_idx[3]: ("0/1", "0/0", "./."),     # mother no-call
        panel7_idx[5]: ("./1", "0/0", "0/1"),     # proband partial no-call
        panel7_idx[7]: ("1|0", "0/0", "0/1"),     # phased
        panel7_idx[9]: ("1/0", "0/0", "0/0"),     # unphased 1/0, no reference reads
    }
    emd0 = len(panel7_idx)
    special.update({
        emd0 + 0: ("1/1", "1", "0/1"),            # X-linked: hemizygous father, carrier mother, affected daughter
        emd0 + 1: ("0/1", "0", "0/1"),            # haploid ref father
        emd0 + 2: ("0/0", ".", "0/0"),            # haploid no-call father
        emd0 + 3: ("0/1", "1", "0/0"),            # father hemizygous alt, mother ref, daughter het
        emd0 + 4: ("0/1", "0", "0/1"),            # deletion, haploid ref father
        emd0 + 5: ("0/0", "0", "0/0"),            # insertion, all reference
    })
    names, weights = zip(*BULK_WEIGHTS)
    out_vcf = out_prefix + ".vcf"
    expected = []
    with pysam.VariantFile(out_vcf, "w", header=header) as out:
        for i, (kind, r, _) in enumerate(records):
            gts = special.get(i) or BULK[rng.choices(names, weights)[0]]
            if kind == "panel7":
                chrom, pos, ref, alt = r.chrom, r.pos, r.ref, r.alts[0]
                rec = out.new_record(contig=chrom, start=r.start, stop=r.stop, alleles=r.alleles, id=r.id,
                                     qual=r.qual, filter=list(r.filter.keys()))
                for k, v in r.info.items():
                    rec.info[k] = v
            else:
                chrom, pos, ref, alt = r
                rec = out.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt), qual=1000.0, filter=["PASS"])
                rec.info["gene"] = "EMD"
                rec.info["variant_key"] = f"{chrom}:{pos}-{ref}-{alt}"
            calls = []
            for sample, gt in zip(SAMPLES, gts):
                gq = 40 if (i % 11 == 0 and sample == SAMPLES[1] and gt not in ("./.", ".")) else None
                calls.append(format_call(gt, rng, gq=gq))
                expected.append((f"{chrom}:{pos}-{ref}-{alt}", sample, gt, status_of(gt)))
            for sample, call in zip(SAMPLES, calls):
                f = dict(zip(("GT", "GQ", "DP", "AD"), call.split(":")))
                gt = f["GT"]
                sep = "|" if "|" in gt else "/"
                alleles = tuple(None if a == "." else int(a) for a in gt.replace("|", "/").split("/"))
                rec.samples[sample]["GT"] = alleles
                rec.samples[sample].phased = sep == "|"
                if f["GQ"] != ".":
                    rec.samples[sample]["GQ"] = int(f["GQ"])
                if f["DP"] != ".":
                    rec.samples[sample]["DP"] = int(f["DP"])
                if f["AD"] != ".":
                    rec.samples[sample]["AD"] = tuple(int(x) for x in f["AD"].split(","))
            out.write(rec)
    pysam.tabix_compress(out_vcf, out_vcf + ".gz", force=True)
    pysam.tabix_index(out_vcf + ".gz", preset="vcf", force=True)
    os.remove(out_vcf)
    with open(out_prefix + ".expected.tsv", "w") as f:
        f.write("Locus\tsample_id\tgt_raw\tgt_status\n")
        f.writelines("\t".join(row) + "\n" for row in expected)
    print(f"{len(records)} records, {len(SAMPLES)} samples -> {out_prefix}.vcf.gz; expectation -> {out_prefix}.expected.tsv")


def compare(genotypes_path, expected_path):
    import pandas as pd
    got = pd.read_parquet(genotypes_path)
    exp = pd.read_csv(expected_path, sep="\t", dtype=str, keep_default_na=False)
    merged = exp.merge(got[["Locus", "sample_id", "gt_raw", "gt_status"]], on=["Locus", "sample_id"], how="outer",
                       suffixes=("_expected", "_got"), indicator=True)
    problems = merged[(merged["_merge"] != "both") | (merged["gt_raw_expected"] != merged["gt_raw_got"])
                      | (merged["gt_status_expected"] != merged["gt_status_got"])]
    print(f"{len(exp)} expected rows, {len(got)} rows in the table, {len(problems)} mismatch(es)")
    if len(problems):
        print(problems.head(20).to_string())
        return 1
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel7", default=PANEL7, help="REF-corrected panel7 VCF (default %(default)s)")
    p.add_argument("--fasta", default=FASTA, help="GRCh38 FASTA the REF alleles come from (default %(default)s)")
    p.add_argument("--out-prefix", default="test_data/raw_vcfs/panel8_trio", help="Output prefix (default %(default)s)")
    p.add_argument("--seed", type=int, default=8)
    p.add_argument("--compare-genotypes", help="Compare this genotype .pq with --expected instead of building")
    p.add_argument("--expected", help="Expected TSV written by the build step")
    a = p.parse_args()
    if a.compare_genotypes:
        sys.exit(compare(a.compare_genotypes, a.expected or a.out_prefix + ".expected.tsv"))
    build(a.panel7, a.fasta, a.out_prefix, a.seed)


if __name__ == "__main__":
    main()
