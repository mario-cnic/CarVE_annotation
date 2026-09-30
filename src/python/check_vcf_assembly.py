#!/usr/bin/env python3
"""Pre-flight check that an input VCF matches the reference FASTA the predictors use.

Two independent checks:

1. Header: every ##contig with a length is resolved to a FASTA contig (exact name, then with the
   `chr` prefix toggled, then MT <-> chrM) and its length compared with the FASTA index. Any
   length mismatch fails. A header whose contigs match nothing in the FASTA also fails.
2. REF spot-check: the REF allele of the first N records is compared with the FASTA sequence at
   the same 1-based position. Fails when the mismatch fraction exceeds the threshold.

A VCF with no ##contig lengths cannot pass check 1, so it passes only if check 2 examined at least
`--min-ref-checked` records and stayed under the threshold; otherwise it fails as unverifiable.

Exit status: 0 = pass (warnings possible), 1 = fail. A TSV report is written to `--report`.
"""

import argparse
import sys

import pysam

ACGT = set("ACGT")


def resolve_contig(name, fasta_lengths):
    """Returns the FASTA contig name matching a VCF contig name, or None."""
    if name in fasta_lengths:
        return name
    alt = name[3:] if name.startswith("chr") else f"chr{name}"
    if alt in fasta_lengths:
        return alt
    mito = {"MT": "chrM", "chrM": "MT", "M": "chrM", "chrMT": "chrM"}.get(name)
    if mito in fasta_lengths:
        return mito
    return None


def check_header(header_contigs, fasta_lengths):
    """Compares header ##contig lengths with the FASTA index.

    Returns (n_with_length, matched, mismatched, unmatched), where mismatched is a list of
    (vcf_name, fasta_name, vcf_length, fasta_length) and unmatched a list of VCF contig names.
    """
    with_length = {n: l for n, l in header_contigs.items() if l is not None}
    matched, mismatched, unmatched = [], [], []
    for name, length in with_length.items():
        fname = resolve_contig(name, fasta_lengths)
        if fname is None:
            unmatched.append(name)
        elif fasta_lengths[fname] != length:
            mismatched.append((name, fname, length, fasta_lengths[fname]))
        else:
            matched.append(name)
    return len(with_length), matched, mismatched, unmatched


def check_ref_alleles(vcf, fasta, fasta_lengths, max_records):
    """Compares REF alleles of the first `max_records` records with the FASTA.

    Records with a non-ACGT REF (symbolic, N) or on a contig absent from the FASTA are skipped.
    A FASTA base outside ACGT (N or IUPAC ambiguity code) is treated as compatible.
    Returns (checked, skipped, mismatches) where mismatches lists (contig, pos, vcf_ref, fasta_ref).
    """
    checked, skipped, mismatches = 0, 0, []
    for i, rec in enumerate(vcf):
        if i >= max_records:
            break
        ref = (rec.ref or "").upper()
        fname = resolve_contig(rec.chrom, fasta_lengths)
        if not ref or set(ref) - ACGT or fname is None or rec.stop > fasta_lengths[fname]:
            skipped += 1
            continue
        fasta_ref = fasta.fetch(fname, rec.start, rec.stop).upper()
        checked += 1
        if any(v != f and f in ACGT for v, f in zip(ref, fasta_ref)):
            mismatches.append((rec.chrom, rec.pos, ref, fasta_ref))
    return checked, skipped, mismatches


def run(vcf_path, fasta_path, report_path, ref_check_records, max_ref_mismatch_frac, min_ref_checked):
    fasta = pysam.FastaFile(fasta_path)
    fasta_lengths = dict(zip(fasta.references, fasta.lengths))
    vcf = pysam.VariantFile(vcf_path)
    header_contigs = {name: c.length for name, c in vcf.header.contigs.items()}

    failures, warnings = [], []

    n_len, matched, mismatched, unmatched = check_header(header_contigs, fasta_lengths)
    if n_len == 0:
        warnings.append("no ##contig lines with a length; relying on the REF spot-check alone")
    else:
        if mismatched:
            shown = "; ".join(f"{v} ({vl}) vs FASTA {f} ({fl})" for v, f, vl, fl in mismatched[:10])
            failures.append(f"{len(mismatched)} ##contig length(s) differ from the FASTA: {shown}")
        if not matched and not mismatched:
            failures.append(f"none of the {n_len} ##contig names resolve to a FASTA contig")
        if unmatched:
            warnings.append(f"{len(unmatched)} ##contig name(s) not in the FASTA: {', '.join(unmatched[:10])}")

    checked, skipped, ref_mm = check_ref_alleles(vcf, fasta, fasta_lengths, ref_check_records)
    frac = len(ref_mm) / checked if checked else 0.0
    if checked and frac > max_ref_mismatch_frac:
        shown = "; ".join(f"{c}:{p} VCF={r} FASTA={f}" for c, p, r, f in ref_mm[:10])
        failures.append(
            f"REF mismatch in {len(ref_mm)}/{checked} records ({frac:.2%} > {max_ref_mismatch_frac:.2%}): {shown}"
        )
    elif ref_mm:
        warnings.append(f"REF mismatch in {len(ref_mm)}/{checked} records (within threshold)")
    if n_len == 0 and checked < min_ref_checked:
        failures.append(
            f"assembly unverifiable: no ##contig lengths and only {checked} REF allele(s) checkable "
            f"(need >= {min_ref_checked}); add ##contig lines, e.g. `bcftools reheader --fai`"
        )

    status = "FAIL" if failures else "PASS"
    with open(report_path, "w") as out:
        out.write("key\tvalue\n")
        rows = [
            ("status", status),
            ("vcf", vcf_path),
            ("fasta", fasta_path),
            ("header_contigs_with_length", n_len),
            ("header_contigs_matched", len(matched)),
            ("header_contigs_length_mismatch", len(mismatched)),
            ("header_contigs_not_in_fasta", len(unmatched)),
            ("ref_records_checked", checked),
            ("ref_records_skipped", skipped),
            ("ref_mismatches", len(ref_mm)),
            ("max_ref_mismatch_frac", max_ref_mismatch_frac),
        ]
        rows += [("failure", f) for f in failures] + [("warning", w) for w in warnings]
        for k, v in rows:
            out.write(f"{k}\t{v}\n")

    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    for f in failures:
        print(f"ERROR: {f}", file=sys.stderr)
    print(
        f"Assembly check {status}: {len(matched)}/{n_len} header contigs match, "
        f"{len(ref_mm)}/{checked} REF mismatches"
    )
    return 1 if failures else 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="Input VCF/BCF (bgzipped or plain)")
    p.add_argument("--fasta", required=True, help="Reference FASTA with a .fai index")
    p.add_argument("--report", required=True, help="Output TSV report path")
    p.add_argument("--ref-check-records", type=int, default=1000, help="Records to REF-check (default 1000)")
    p.add_argument("--max-ref-mismatch-frac", type=float, default=0.01,
                   help="Fail above this REF mismatch fraction (default 0.01)")
    p.add_argument("--min-ref-checked", type=int, default=20,
                   help="Minimum REF-checked records when ##contig lengths are absent (default 20)")
    a = p.parse_args()
    sys.exit(run(a.vcf, a.fasta, a.report, a.ref_check_records, a.max_ref_mismatch_frac, a.min_ref_checked))


if __name__ == "__main__":
    main()
