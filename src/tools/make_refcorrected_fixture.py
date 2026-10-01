#!/usr/bin/env python3
"""Make a reference-consistent copy of a synthetic SNV test VCF.

The original test/raw_vcfs/panel7_test.vcf.gz was generated with arbitrary REF alleles (35 of 72 do not match GRCh38 and
5 have REF == ALT), so the pre-flight assembly check (MISC-14) rightly rejects it. This writes a NEW file where
REF is the reference base at POS; the original ALT is kept unless it equals the new REF, in which case the next base
in A->C->G->T order is used (deterministic). `variant_key` in INFO is rewritten to match. Positions, genes, QUAL, FILTER
and other INFO values are untouched. SNV-only; the original file is never modified.
"""
import argparse, gzip, re, subprocess
import pysam

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--vcf", required=True); ap.add_argument("--fasta", required=True)
ap.add_argument("--out-vcf", required=True, help="plain .vcf path; bgzipped + indexed to .vcf.gz")
ap.add_argument("--bgzip", default="bgzip"); ap.add_argument("--tabix", default="tabix")
a = ap.parse_args()

fa = pysam.FastaFile(a.fasta)
names = set(fa.references)
n = changed_ref = changed_alt = 0
out = open(a.out_vcf, "w")
for ln in gzip.open(a.vcf, "rt"):
    if ln.startswith("##contig") or ln.startswith("##fileformat"):
        out.write(ln); continue
    if ln.startswith("#CHROM"):
        out.write(f"##reference_consistent_copy=source:{a.vcf};reference:{a.fasta};REF set to the reference base, ALT kept unless equal to REF (then next of A,C,G,T)\n"); out.write(ln); continue
    if ln.startswith("#"):
        out.write(ln); continue
    f = ln.rstrip("\n").split("\t")
    chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
    if len(ref) != 1 or len(alt) != 1:
        raise SystemExit(f"SNV-only tool; got {chrom}:{pos} {ref}>{alt}")
    contig = chrom if chrom in names else "chr" + chrom
    base = fa.fetch(contig, pos - 1, pos).upper()
    if base not in "ACGT":
        raise SystemExit(f"non-ACGT reference base at {contig}:{pos}")
    new_alt = alt if alt != base else "ACGT"[("ACGT".index(base) + 1) % 4]
    changed_ref += base != ref; changed_alt += new_alt != alt
    f[3], f[4] = base, new_alt
    f[7] = re.sub(r"variant_key=[^;]*", f"variant_key={chrom}:{pos}-{base}-{new_alt}", f[7])
    out.write("\t".join(f) + "\n"); n += 1
out.close()
subprocess.run([a.bgzip, "-f", a.out_vcf], check=True)
subprocess.run([a.tabix, "-f", "-p", "vcf", a.out_vcf + ".gz"], check=True)
print(f"{n} records; REF changed on {changed_ref}; ALT changed on {changed_alt}; wrote {a.out_vcf}.gz")
