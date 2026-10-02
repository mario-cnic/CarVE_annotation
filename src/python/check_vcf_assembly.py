#!/usr/bin/env python3
"""Pre-flight check of an input VCF: reference assembly, record normalisation and sample columns.

Assembly (two independent checks):

1. Header: every ##contig with a length is resolved to a FASTA contig (exact name, then with the
   `chr` prefix toggled, then MT <-> chrM) and its length compared with the FASTA index. Any
   length mismatch fails. A header whose contigs match nothing in the FASTA also fails.
2. REF spot-check: the REF allele of the first N records is compared with the FASTA sequence at
   the same 1-based position. Fails when the mismatch fraction exceeds the threshold.

A VCF with no ##contig lengths cannot pass check 1, so it passes only if check 2 examined at least
`--min-ref-checked` records and stayed under the threshold; otherwise it fails as unverifiable.

Records (one pass over the whole file, no FASTA needed):

3. Multiallelic records (more than one ALT) fail above `--max-multiallelic`.
4. Biallelic indels/MNVs that are not left-aligned or not parsimonious fail above
   `--max-unnormalized`: a record is not normalised when REF and ALT end with the same base, or
   both are longer than one base and start with the same base.
5. gVCF input (`<NON_REF>` / `<*>` ALT, or `##GVCFBlock` header lines) always fails.
6. Duplicate chrom:pos:ref:alt records, and records whose ALT equals REF, only warn.

Samples: duplicate sample names fail. The report records the genotype mode (`no_samples`,
`samples_without_GT`, `single_sample`, `multi_sample`), the sample names and the declared FORMAT
fields; nothing is changed downstream by them.

Exit status: 0 = pass (warnings possible), 1 = fail. A TSV report is written to `--report`.
"""

import argparse
import gzip
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


def header_sample_names(vcf_path):
    """Sample names from the #CHROM line, read without htslib (which refuses duplicate names)."""
    opener = gzip.open if vcf_path.endswith((".gz", ".bgz")) else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    return []


def is_normalised(ref, alt):
    """True when a biallelic ACGT record is trimmed and left-aligned (VCF anchor-base convention)."""
    if ref[-1] == alt[-1]:
        return False
    return not (len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0])


def scan_records(vcf):
    """One pass over all records. Returns a dict of counts."""
    c = dict(records=0, multiallelic=0, not_normalised=0, symbolic=0, gvcf=0, no_alt=0, duplicates=0, alt_is_ref=0)
    group, seen = None, set()
    for rec in vcf:
        c["records"] += 1
        alts = rec.alts or ()
        key = (rec.chrom, rec.pos)
        if key != group:
            group, seen = key, set()
        if not alts:
            c["no_alt"] += 1
            continue
        if any(a in ("<NON_REF>", "<*>") for a in alts):
            c["gvcf"] += 1
        sig = (rec.ref, alts)
        if sig in seen:
            c["duplicates"] += 1
        seen.add(sig)
        if len(alts) > 1:
            c["multiallelic"] += 1
            continue
        ref, alt = rec.ref.upper(), alts[0].upper()
        if set(ref) - ACGT or set(alt) - ACGT:
            c["symbolic"] += 1
        elif ref == alt:
            c["alt_is_ref"] += 1
        elif not is_normalised(ref, alt):
            c["not_normalised"] += 1
    return c


def genotype_mode(samples, format_ids):
    """Genotype mode of the input from its header."""
    if not samples:
        return "no_samples"
    if "GT" not in format_ids:
        return "samples_without_GT"
    return "single_sample" if len(samples) == 1 else "multi_sample"


def write_report(report_path, rows, failures, warnings):
    with open(report_path, "w") as out:
        out.write("key\tvalue\n")
        for k, v in rows + [("failure", f) for f in failures] + [("warning", w) for w in warnings]:
            out.write(f"{k}\t{v}\n")


def run(vcf_path, fasta_path, report_path, ref_check_records, max_ref_mismatch_frac, min_ref_checked,
        max_multiallelic=0, max_unnormalised=0):
    fasta = pysam.FastaFile(fasta_path)
    fasta_lengths = dict(zip(fasta.references, fasta.lengths))
    try:
        vcf = pysam.VariantFile(vcf_path)
    except (ValueError, OSError) as exc:
        names = header_sample_names(vcf_path)
        dups = sorted({n for n in names if names.count(n) > 1})
        msg = (f"duplicate sample name(s) in the #CHROM line: {', '.join(dups)}" if dups
               else f"VCF header could not be read: {exc}")
        write_report(report_path, [("status", "FAIL"), ("vcf", vcf_path), ("fasta", fasta_path)], [msg], [])
        print(f"ERROR: {msg}", file=sys.stderr)
        return 1
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

    # Record scan reopens the file: the REF spot-check above consumed the first records.
    scan = scan_records(pysam.VariantFile(vcf_path))
    fix = f"bcftools norm -m -any -f {fasta_path}"
    if scan["gvcf"] or any(k.startswith("GVCFBlock") for k in (r.key for r in vcf.header.records)):
        failures.append("input is a gVCF (<NON_REF>/<*> ALT or ##GVCFBlock header); use the genotyped VCF")
    if scan["multiallelic"] > max_multiallelic:
        failures.append(f"{scan['multiallelic']} multiallelic record(s) (allowed {max_multiallelic}); split with `{fix}`")
    if scan["not_normalised"] > max_unnormalised:
        failures.append(
            f"{scan['not_normalised']} indel/MNV record(s) not left-aligned or not parsimonious "
            f"(allowed {max_unnormalised}); normalise with `{fix}`"
        )
    if scan["duplicates"]:
        warnings.append(f"{scan['duplicates']} duplicate chrom:pos:ref:alt record(s)")
    if scan["alt_is_ref"]:
        warnings.append(f"{scan['alt_is_ref']} record(s) with ALT identical to REF")
    if scan["no_alt"]:
        warnings.append(f"{scan['no_alt']} record(s) without an ALT allele")

    samples = list(vcf.header.samples)
    dup_names = sorted({n for n in samples if samples.count(n) > 1})
    if dup_names:
        failures.append(f"duplicate sample name(s): {', '.join(dup_names)}")
    format_ids = list(vcf.header.formats.keys())
    mode = genotype_mode(samples, format_ids)
    if mode == "samples_without_GT":
        warnings.append("sample columns present but no GT FORMAT field declared")

    status = "FAIL" if failures else "PASS"
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
        ("records_scanned", scan["records"]),
        ("multiallelic_records", scan["multiallelic"]),
        ("not_normalised_records", scan["not_normalised"]),
        ("symbolic_or_non_acgt_records", scan["symbolic"]),
        ("duplicate_site_records", scan["duplicates"]),
        ("alt_equals_ref_records", scan["alt_is_ref"]),
        ("genotype_mode", mode),
        ("n_samples", len(samples)),
        ("sample_names", ",".join(samples)),
        ("format_fields_declared", ",".join(f for f in ("GT", "GQ", "DP", "AD", "PL", "PS", "FT") if f in format_ids)),
    ]
    write_report(report_path, rows, failures, warnings)

    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    for f in failures:
        print(f"ERROR: {f}", file=sys.stderr)
    print(
        f"Input check {status}: {len(matched)}/{n_len} header contigs match, "
        f"{len(ref_mm)}/{checked} REF mismatches, {scan['records']} records "
        f"({scan['multiallelic']} multiallelic, {scan['not_normalised']} not normalised), "
        f"genotype mode {mode} ({len(samples)} sample(s))"
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
    p.add_argument("--max-multiallelic", type=int, default=0,
                   help="Fail above this many multiallelic records (default 0)")
    p.add_argument("--max-unnormalised", type=int, default=0,
                   help="Fail above this many not-left-aligned / non-parsimonious records (default 0)")
    a = p.parse_args()
    sys.exit(run(a.vcf, a.fasta, a.report, a.ref_check_records, a.max_ref_mismatch_frac, a.min_ref_checked,
                 a.max_multiallelic, a.max_unnormalised))


if __name__ == "__main__":
    main()
