import os
import random
import sys

import pysam
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

from check_vcf_assembly import run, resolve_contig

CHR1_LEN, CHR2_LEN = 3000, 2000


@pytest.fixture
def fasta(tmp_path):
    rng = random.Random(7)
    seqs = {"chr1": "".join(rng.choice("ACGT") for _ in range(CHR1_LEN)),
            "chr2": "".join(rng.choice("ACGT") for _ in range(CHR2_LEN)),
            "chrM": "ACGT" * 50}
    seqs["chr2"] = seqs["chr2"][:99] + "R" + seqs["chr2"][100:]  # IUPAC code at 1-based pos 100
    path = tmp_path / "ref.fa"
    with open(path, "w") as f:
        for name, s in seqs.items():
            f.write(f">{name}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")
    pysam.faidx(str(path))
    return str(path), seqs


def write_vcf(path, contig_lines, records):
    """records: list of (chrom, 1-based pos, ref, alt). Writes bgzipped + tabix-indexed VCF."""
    txt = str(path)[:-3]
    with open(txt, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in contig_lines:
            f.write(c + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in records:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    pysam.tabix_index(txt, preset="vcf", force=True)
    return txt + ".gz"


def true_records(seqs, chrom="chr1", n=30, label=None):
    """SNVs whose REF is the FASTA base at that position."""
    recs = []
    for i in range(n):
        pos = 100 + i * 50
        ref = seqs[chrom][pos - 1]
        alt = "A" if ref != "A" else "C"
        recs.append((label or chrom, pos, ref, alt))
    return recs


def wrong_records(seqs, chrom="chr1", n=30):
    """SNVs whose REF disagrees with the FASTA."""
    return [(c, p, {"A": "C", "C": "G", "G": "T", "T": "A"}[r], r) for c, p, r, _ in true_records(seqs, chrom, n)]


def report(tmp_path):
    return dict(line.rstrip("\n").split("\t", 1) for line in open(tmp_path / "r.tsv") if not line.startswith("key\t"))


def go(tmp_path, vcf, fasta_path, **kw):
    args = dict(ref_check_records=1000, max_ref_mismatch_frac=0.01, min_ref_checked=20)
    args.update(kw)
    return run(vcf, fasta_path, str(tmp_path / "r.tsv"), **args)


GRCH38_LIKE = [f"##contig=<ID=chr1,length={CHR1_LEN}>", f"##contig=<ID=chr2,length={CHR2_LEN}>"]


def test_resolve_contig_prefix_and_mito():
    lengths = {"chr1": 1, "chrM": 1}
    assert resolve_contig("chr1", lengths) == "chr1"
    assert resolve_contig("1", lengths) == "chr1"
    assert resolve_contig("MT", lengths) == "chrM"
    assert resolve_contig("GL000192.1", lengths) is None


def test_matching_header_and_refs_pass(tmp_path, fasta):
    fa, seqs = fasta
    vcf = write_vcf(tmp_path / "ok.vcf.gz", GRCH38_LIKE, true_records(seqs))
    assert go(tmp_path, vcf, fa) == 0
    r = report(tmp_path)
    assert r["status"] == "PASS" and r["header_contigs_matched"] == "2" and r["ref_mismatches"] == "0"


def test_bare_contig_names_pass(tmp_path, fasta):
    fa, seqs = fasta
    header = [f"##contig=<ID=1,length={CHR1_LEN}>", f"##contig=<ID=2,length={CHR2_LEN}>"]
    vcf = write_vcf(tmp_path / "bare.vcf.gz", header, true_records(seqs, label="1"))
    assert go(tmp_path, vcf, fa) == 0


def test_grch37_like_lengths_fail(tmp_path, fasta):
    fa, seqs = fasta
    header = [f"##contig=<ID=chr1,length={CHR1_LEN + 294199}>", f"##contig=<ID=chr2,length={CHR2_LEN}>"]
    vcf = write_vcf(tmp_path / "b37.vcf.gz", header, true_records(seqs))
    assert go(tmp_path, vcf, fa) == 1
    r = report(tmp_path)
    assert r["status"] == "FAIL" and r["header_contigs_length_mismatch"] == "1"


def test_no_contig_names_resolve_fail(tmp_path, fasta):
    fa, _ = fasta
    header = ["##contig=<ID=NC_000001.11,length=248956422>"]
    vcf = write_vcf(tmp_path / "refseq.vcf.gz", header, [])
    assert go(tmp_path, vcf, fa) == 1


def test_ref_mismatch_with_good_header_fails(tmp_path, fasta):
    fa, seqs = fasta
    vcf = write_vcf(tmp_path / "badref.vcf.gz", GRCH38_LIKE, wrong_records(seqs))
    assert go(tmp_path, vcf, fa) == 1
    assert report(tmp_path)["ref_mismatches"] == "30"


def test_single_ref_mismatch_under_threshold_passes_with_warning(tmp_path, fasta):
    fa, seqs = fasta
    # 1 mismatch in 79 checked records (chr2:100 is the IUPAC position and is skipped as a REF).
    recs = true_records(seqs, n=50) + true_records(seqs, chrom="chr2", n=30)[1:]
    recs[0] = wrong_records(seqs, n=1)[0]
    vcf = write_vcf(tmp_path / "one.vcf.gz", GRCH38_LIKE, recs)
    assert go(tmp_path, vcf, fa, max_ref_mismatch_frac=0.02) == 0
    r = report(tmp_path)
    assert r["ref_mismatches"] == "1" and "warning" in r


def test_missing_contig_lines_pass_when_refs_match(tmp_path, fasta):
    fa, seqs = fasta
    vcf = write_vcf(tmp_path / "nohdr.vcf.gz", [], true_records(seqs))
    assert go(tmp_path, vcf, fa) == 0
    r = report(tmp_path)
    assert r["header_contigs_with_length"] == "0" and r["ref_records_checked"] == "30"


def test_missing_contig_lines_fail_when_refs_mismatch(tmp_path, fasta):
    fa, seqs = fasta
    vcf = write_vcf(tmp_path / "nohdr_bad.vcf.gz", [], wrong_records(seqs))
    assert go(tmp_path, vcf, fa) == 1


def test_missing_contig_lines_fail_when_too_few_records(tmp_path, fasta):
    fa, seqs = fasta
    vcf = write_vcf(tmp_path / "nohdr_few.vcf.gz", [], true_records(seqs, n=5))
    assert go(tmp_path, vcf, fa) == 1
    assert "unverifiable" in open(tmp_path / "r.tsv").read()


def test_iupac_fasta_base_is_not_a_mismatch(tmp_path, fasta):
    fa, seqs = fasta
    recs = true_records(seqs) + [("chr2", 100, "A", "G")]
    vcf = write_vcf(tmp_path / "iupac.vcf.gz", GRCH38_LIKE, recs)
    assert go(tmp_path, vcf, fa, max_ref_mismatch_frac=0.0) == 0


def test_symbolic_and_out_of_range_records_skipped(tmp_path, fasta):
    fa, seqs = fasta
    recs = true_records(seqs) + [("chr2", 10, "N", "<DEL>"), ("chr2", CHR2_LEN + 5, "A", "G")]
    vcf = write_vcf(tmp_path / "skip.vcf.gz", GRCH38_LIKE, recs)
    assert go(tmp_path, vcf, fa) == 0
    assert report(tmp_path)["ref_records_skipped"] == "2"


def write_vcf_samples(path, records, samples=(), format_ids=("GT",), gt="0/1"):
    """records: (chrom, pos, ref, alt-string). One sample column per name; FORMAT = format_ids."""
    records = sorted(records, key=lambda r: (r[0], r[1]))
    txt = str(path)[:-3]
    with open(txt, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in GRCH38_LIKE:
            f.write(c + "\n")
        for i in format_ids:
            f.write(f'##FORMAT=<ID={i},Number=1,Type=String,Description="{i}">\n')
        cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        cols += ("\tFORMAT\t" + "\t".join(samples)) if samples else ""
        f.write(cols + "\n")
        for chrom, pos, ref, alt in records:
            line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t."
            if samples:
                line += "\t" + ":".join(format_ids) + "".join("\t" + (gt if i == 0 else "0") for i in range(len(samples)))
            f.write(line + "\n")
    pysam.tabix_index(txt, preset="vcf", force=True)
    return txt + ".gz"


def indel_records(seqs, normalised, n=5):
    """Deletions anchored on true FASTA bases; REF last base differs from ALT only when normalised."""
    recs = []
    for i in range(n):
        pos = 100 + i * 50
        b1, b2 = seqs["chr1"][pos - 1], seqs["chr1"][pos]
        if (b1 != b2) == normalised:
            recs.append(("chr1", pos, b1 + b2, b1 if normalised else b2))
    return recs


def test_multiallelic_record_fails_and_can_be_allowed(tmp_path, fasta):
    fa, seqs = fasta
    ref = seqs["chr1"][99]
    alts = ",".join(b for b in "ACGT" if b != ref)
    vcf = write_vcf_samples(tmp_path / "multi.vcf.gz", true_records(seqs, n=29) + [("chr1", 100, ref, alts)])
    assert go(tmp_path, vcf, fa) == 1
    r = report(tmp_path)
    assert r["multiallelic_records"] == "1" and "bcftools norm -m -any" in open(tmp_path / "r.tsv").read()
    assert go(tmp_path, vcf, fa, max_multiallelic=1) == 0


def test_non_left_aligned_and_non_parsimonious_records_fail(tmp_path, fasta):
    fa, seqs = fasta
    b = seqs["chr1"][1199:1202]
    recs = true_records(seqs, n=29) + [("chr1", 1200, b[:2], b[1]),            # REF/ALT end with the same base
                                       ("chr1", 1300, seqs["chr1"][1299:1302], seqs["chr1"][1299:1301] + "A")]
    if recs[-1][3] == recs[-1][2]:
        recs[-1] = (recs[-1][0], recs[-1][1], recs[-1][2], recs[-1][3][:2] + "C")  # shared first base, differing last
    vcf = write_vcf_samples(tmp_path / "unnorm.vcf.gz", recs)
    assert go(tmp_path, vcf, fa) == 1
    assert report(tmp_path)["not_normalised_records"] == "2"


def test_normalised_indel_and_mnv_pass(tmp_path, fasta):
    fa, seqs = fasta
    recs = true_records(seqs, n=20) + indel_records(seqs, normalised=True)
    pos = 1400
    ref = seqs["chr1"][pos - 1:pos + 1]
    recs.append(("chr1", pos, ref, ("T" if ref[0] != "T" else "G") + ("A" if ref[1] != "A" else "C")))  # MNV
    vcf = write_vcf_samples(tmp_path / "norm.vcf.gz", recs)
    assert go(tmp_path, vcf, fa) == 0
    assert report(tmp_path)["not_normalised_records"] == "0"


def test_gvcf_fails(tmp_path, fasta):
    fa, seqs = fasta
    ref = seqs["chr1"][99]
    vcf = write_vcf_samples(tmp_path / "g.vcf.gz", true_records(seqs, n=29) + [("chr1", 100, ref, "<NON_REF>")])
    assert go(tmp_path, vcf, fa) == 1
    assert "gVCF" in open(tmp_path / "r.tsv").read()


def test_duplicate_sites_and_alt_equal_ref_only_warn(tmp_path, fasta):
    fa, seqs = fasta
    recs = true_records(seqs)
    recs.append(recs[0])
    recs.append(("chr2", 300, seqs["chr2"][299], seqs["chr2"][299]))
    vcf = write_vcf(tmp_path / "dup.vcf.gz", GRCH38_LIKE, sorted(recs, key=lambda r: (r[0], r[1])))
    assert go(tmp_path, vcf, fa) == 0
    r = report(tmp_path)
    assert r["duplicate_site_records"] == "1" and r["alt_equals_ref_records"] == "1" and r["not_normalised_records"] == "0"


def test_duplicate_sample_names_fail_with_clear_message(tmp_path, fasta):
    fa, seqs = fasta
    vcf = write_vcf_samples(tmp_path / "dupname.vcf.gz", true_records(seqs), samples=("A", "A"))
    assert go(tmp_path, vcf, fa) == 1
    assert "duplicate sample name(s) in the #CHROM line: A" in open(tmp_path / "r.tsv").read()


@pytest.mark.parametrize("samples,format_ids,mode,n", [
    ((), ("GT",), "no_samples", "0"),
    (("S1",), ("GT", "GQ"), "single_sample", "1"),
    (("MOTHER", "INDEX", "FATHER"), ("GT", "DP"), "multi_sample", "3"),
    (("S1",), ("DP",), "samples_without_GT", "1"),
])
def test_genotype_mode_and_samples_reported(tmp_path, fasta, samples, format_ids, mode, n):
    fa, seqs = fasta
    vcf = write_vcf_samples(tmp_path / "gt.vcf.gz", true_records(seqs), samples=samples, format_ids=format_ids)
    assert go(tmp_path, vcf, fa) == 0
    r = report(tmp_path)
    assert r["genotype_mode"] == mode and r["n_samples"] == n and r["sample_names"] == ",".join(samples)
    assert r["format_fields_declared"] == ",".join(f for f in ("GT", "GQ", "DP") if f in format_ids)
