"""annotate_branchpointer.py: positional LaBranchoR lookup only (PRED-18)."""
import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python")))
import annotate_branchpointer as ab  # noqa: E402

# chrom, bp_start, bp_end, acceptor, score, strand
BED = [("chr1", 1000, 1001, 1030, 0.91, "+"), ("chr1", 5000, 5001, 4970, 0.35, "-")]


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
    write_vcf(vcf, [("chr1", 1001), ("chr1", 5001), ("chr1", 3000), ("chr2", 1001)])
    ab.annotate_branchpoint_vcf(str(vcf), str(out), bed)
    with pysam.VariantFile(str(out)) as vf:
        assert set(vf.header.info) == {"LaBranchoR_score", "LaBranchoR_acc_dist"}
        got = [(r.pos, r.info.get("LaBranchoR_score"), r.info.get("LaBranchoR_acc_dist")) for r in vf]
    # on each predicted branchpoint; off any; chr2 not in the BED index
    assert got == [(1001, pytest.approx(0.91), 29), (5001, pytest.approx(0.35), 31), (3000, None, None), (1001, None, None)]


def test_missing_bed_fails(tmp_path):
    vcf = tmp_path / "in.vcf"
    write_vcf(vcf, [("chr1", 1001)])
    with pytest.raises(FileNotFoundError):
        ab.annotate_branchpoint_vcf(str(vcf), str(tmp_path / "out.vcf"), str(tmp_path / "absent.bed.gz"))
