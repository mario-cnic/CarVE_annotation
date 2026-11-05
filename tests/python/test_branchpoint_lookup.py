"""annotate_branchpointer.py: positional LaBranchoR lookup only (PRED-18)."""
import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python")))
import annotate_branchpointer as ab  # noqa: E402

# chrom, bp_start, bp_end, acceptor (0-based exon first base), score, strand. Branchpoints at 1-based 1001
# (+, exon starts at 1031) and 5001 (-, exon starts at 4971): both 30 bases from the 3'SS. The 8001 one is
# 100 bases away, so it matches only on the branchpoint itself.
BED = [("chr1", 1000, 1001, 1030, 0.91, "+"), ("chr1", 5000, 5001, 4970, 0.35, "-"),
       ("chr1", 8000, 8001, 8100, 0.50, "+")]


@pytest.fixture
def bed(tmp_path):
    plain = tmp_path / "lab.bed"
    plain.write_text("".join("\t".join(map(str, r)) + "\n" for r in BED))
    pysam.tabix_compress(str(plain), str(plain) + ".gz", force=True)
    pysam.tabix_index(str(plain) + ".gz", preset="bed", force=True)
    return str(plain) + ".gz"


def write_vcf(path, records, contigs=("chr1", "chr2")):
    header = pysam.VariantHeader()
    for c in contigs:
        header.add_line(f"##contig=<ID={c},length=100000>")
    with pysam.VariantFile(str(path), "w", header=header) as out:
        for chrom, pos in records:
            out.write(out.new_record(contig=chrom, start=pos - 1, alleles=("A", "G")))


def test_lookup_fields_only(tmp_path, bed):
    vcf, out = tmp_path / "in.vcf", tmp_path / "out.vcf"
    write_vcf(vcf, [("chr1", 1001), ("chr1", 5001), ("chr1", 3000), ("chr1", 8000), ("chr1", 8001), ("chr2", 1001)])
    ab.annotate_branchpoint_vcf(str(vcf), str(out), bed)
    with pysam.VariantFile(str(out)) as vf:
        assert set(vf.header.info) == {"LaBranchoR_score", "LaBranchoR_acc_dist"} | {t for t, _, _ in ab.ISM_INFO}
        got = [(r.pos, r.info.get("LaBranchoR_score"), r.info.get("LaBranchoR_acc_dist")) for r in vf]
    # on each predicted branchpoint; off any; the base before an out-of-range branchpoint; on it; chr2 not in the BED index
    assert got == [(1001, pytest.approx(0.91), 30), (5001, pytest.approx(0.35), 30), (3000, None, None),
                   (8000, None, None), (8001, pytest.approx(0.50), 100), (1001, None, None)]


def test_acc_dist_matches_ism_dist(tmp_path, bed, ism_table):
    # the + BED branchpoint and the + ISM row share the 3'SS (1030): both distances count the same bases
    vcf, out = tmp_path / "in.vcf", tmp_path / "out.vcf"
    write_alleles(vcf, [("chr1", 1001, "A", "T")])
    ab.annotate_branchpoint_vcf(str(vcf), str(out), bed, ism_table)
    with pysam.VariantFile(str(out)) as vf:
        r = next(iter(vf))
        assert r.info["LaBranchoR_acc_dist"] == r.info["LaBranchoR_ISM_dist"] == 30


def test_missing_bed_fails(tmp_path):
    vcf = tmp_path / "in.vcf"
    write_vcf(vcf, [("chr1", 1001)])
    with pytest.raises(FileNotFoundError):
        ab.annotate_branchpoint_vcf(str(vcf), str(tmp_path / "out.vcf"), str(tmp_path / "absent.bed.gz"))


ISM = [  # chrom, pos, ref, alt, strand, 3pSS, d_best_bp, d_max_bp
    ("chr1", 1001, "A", "G", "+", 1030, -0.80, -0.30),     # in the window of two 3'SS (+ and -)
    ("chr1", 1001, "A", "G", "-", 990, -0.10, 0.20),
    ("chr1", 1001, "A", "T", "+", 1030, -0.05, 0.00),
    ("chr1", 2000, "C", "T", "+", 2010, 0.00, 0.40),
]


@pytest.fixture
def ism_table(tmp_path):
    plain = tmp_path / "ism.tsv"
    plain.write_text("#chrom\tpos\tref\talt\tstrand\t3pSS\td_best_bp\td_max_bp\n"
                     + "".join("\t".join(map(str, r)) + "\n" for r in ISM))
    pysam.tabix_compress(str(plain), str(plain) + ".gz", force=True)
    pysam.tabix_index(str(plain) + ".gz", seq_col=0, start_col=1, end_col=1, meta_char="#", zerobased=False, force=True)
    return str(plain) + ".gz"


def write_alleles(path, records):
    header = pysam.VariantHeader()
    for c in ("chr1", "chr2"):
        header.add_line(f"##contig=<ID={c},length=100000>")
    with pysam.VariantFile(str(path), "w", header=header) as out:
        for chrom, pos, ref, alt in records:
            out.write(out.new_record(contig=chrom, start=pos - 1, alleles=(ref, alt)))


def test_ism_scores_are_allele_matched(tmp_path, bed, ism_table):
    vcf, out = tmp_path / "in.vcf", tmp_path / "out.vcf"
    write_alleles(vcf, [("chr1", 1001, "A", "G"), ("chr1", 1001, "A", "C"), ("chr1", 2000, "C", "T"),
                        ("chr1", 2000, "C", "CT"), ("chr2", 1001, "A", "G")])
    ab.annotate_branchpoint_vcf(str(vcf), str(out), bed, ism_table)
    tags = [t for t, _, _ in ab.ISM_INFO]
    with pysam.VariantFile(str(out)) as vf:
        assert set(tags) <= set(vf.header.info)
        got = [tuple(r.info.get(t) for t in tags) for r in vf]
    # A>G: two 3'SS, the one with the lowest d_best (the '+' one, 30 bases upstream) is reported
    assert got[0][2:] == (30, 2) and got[0][0] == pytest.approx(-0.80) and got[0][1] == pytest.approx(-0.30)
    assert got[1] == (None,) * 4                                     # A>C: no ISM row for that ALT
    assert got[2][2:] == (11, 1) and got[2][1] == pytest.approx(0.40)
    assert got[3] == (None,) * 4 and got[4] == (None,) * 4          # indel; contig absent from the table


def test_ism_headers_declared_without_the_table(tmp_path, bed):
    vcf, out = tmp_path / "in.vcf", tmp_path / "out.vcf"
    write_alleles(vcf, [("chr1", 1001, "A", "G")])
    ab.annotate_branchpoint_vcf(str(vcf), str(out), bed)
    with pysam.VariantFile(str(out)) as vf:
        assert {t for t, _, _ in ab.ISM_INFO} <= set(vf.header.info)


def test_missing_ism_table_fails(tmp_path, bed):
    vcf = tmp_path / "in.vcf"
    write_alleles(vcf, [("chr1", 1001, "A", "G")])
    with pytest.raises(FileNotFoundError):
        ab.annotate_branchpoint_vcf(str(vcf), str(tmp_path / "out.vcf"), bed, str(tmp_path / "absent.tsv.gz"))
