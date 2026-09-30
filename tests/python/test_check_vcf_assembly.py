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
