"""hg19 -> GRCh38 liftover in universal_variant_converter.py: strand, REF check, FORMAT kept (LIFT-1, LIFT-2).

A fake converter replaces pyliftover (which downloads UCSC chains): hg19 chr1:p maps to GRCh38 chr1:p+10 on
'+'; hg19 chr2:p maps to GRCh38 chr1:(80-p) on '-' (an inverted region); hg19 chr3 is unmapped.
"""
import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python")))
import universal_variant_converter as uvc  # noqa: E402

SEQ = "ACGTTGCAAC" * 10   # GRCh38 chr1, 100 bp


class FakeLiftOver:
    def convert_coordinate(self, chrom, pos0):
        if chrom == "chr1":
            return [("chr1", pos0 + 10, "+", 1)]
        if chrom == "chr2":
            return [("chr1", 80 - pos0, "-", 1)]
        return []


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_text(">chr1\n" + SEQ + "\n")
    pysam.faidx(str(path))
    return pysam.FastaFile(str(path))


def g38(pos0, n=1):
    return SEQ[pos0:pos0 + n]


def test_plus_strand_snv_and_ref_check(fasta):
    ref = g38(14)                                                    # hg19 chr1:5 -> GRCh38 0-based 14
    assert uvc.lift_variant(FakeLiftOver(), fasta, "1", 5, ref, "A" if ref != "A" else "C")[0][:3] == ("1", 15, ref)
    wrong = "A" if ref != "A" else "C"
    assert uvc.lift_variant(FakeLiftOver(), fasta, "1", 5, wrong, "G") == (None, "ref_mismatch_after_liftover")


def test_minus_strand_snv_and_mnp_are_reverse_complemented(fasta):
    # hg19 chr2:11 -> GRCh38 0-based 70; on the hg19 strand the REF base is the complement
    ref19 = uvc.revcomp(g38(70))
    lifted, why = uvc.lift_variant(FakeLiftOver(), fasta, "2", 11, ref19, "A")
    assert why is None and lifted == ("1", 71, g38(70), "T")
    # two-base MNP at hg19 chr2:11-12 -> GRCh38 69..70, reverse-complemented
    lifted, why = uvc.lift_variant(FakeLiftOver(), fasta, "2", 11, uvc.revcomp(g38(69, 2)), "GA")
    assert why is None and lifted == ("1", 70, g38(69, 2), "TC")


def test_refused_cases(fasta):
    assert uvc.lift_variant(FakeLiftOver(), fasta, "2", 11, uvc.revcomp(g38(70)), "AT")[1] == "liftover_minus_strand_indel"
    assert uvc.lift_variant(FakeLiftOver(), fasta, "3", 5, "A", "G") == (None, "liftover_failed_unmapped")


def test_vcf_liftover_keeps_samples_and_rewrites_the_header(tmp_path, fasta):
    inp = tmp_path / "in.hg19.vcf"
    inp.write_text(
        "##fileformat=VCFv4.2\n##reference=hg19.fa\n##contig=<ID=1,length=249250621>\n##contig=<ID=2,length=243199373>\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        f"2\t11\tminus\t{uvc.revcomp(g38(70))}\tA\t50\tPASS\tDP=7\tGT\t0/1\t1/1\n"
        f"1\t5\tplus\t{g38(14)}\tT\t60\tPASS\t.\tGT\t0/0\t0/1\n"
        "3\t5\tunmapped\tA\tG\t10\tPASS\t.\tGT\t0/1\t0/0\n")
    out = tmp_path / "out.vcf"
    n, errors = uvc.liftover_vcf(str(inp), str(out), FakeLiftOver(), fasta)
    assert n == 2 and [e[2] for e in errors] == ["liftover_failed_unmapped"]
    text = out.read_text()
    assert "length=249250621" not in text and "##reference=hg19.fa" not in text and "##contig=<ID=1,length=248956422>" in text
    assert '##INFO=<ID=DP,' in text and '##INFO=<ID=ORIG_HG19_POS,' in text and "\tFORMAT\tS1\tS2\n" in text
    body = [l.split("\t") for l in text.splitlines() if not l.startswith("#")]
    assert [r[2] for r in body] == ["plus", "minus"]                  # GRCh38 position order
    assert body[0][:5] == ["1", "15", "plus", g38(14), "T"] and body[0][7] == "ORIG_HG19_POS=5"
    assert body[1][:5] == ["1", "71", "minus", g38(70), "T"] and body[1][7] == "ORIG_HG19_POS=11;DP=7"
    assert body[1][8:] == ["GT", "0/1", "1/1"]                        # FORMAT and samples kept (LIFT-1)
    assert (tmp_path / "out.vcf.err").read_text().splitlines()[1].endswith("liftover_failed_unmapped")


def test_mitochondrial_query_uses_the_ucsc_name(fasta):
    seen = []

    class Spy:
        def convert_coordinate(self, chrom, pos0):
            seen.append(chrom)
            return []
    uvc.lift_variant(Spy(), fasta, "MT", 100, "A", "G")
    assert seen == ["chrM"]
