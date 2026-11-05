"""LaBranchoR ISM hg19 -> GRCh38 builder: clean lifts only, 1-based in and out, REF checked."""
import gzip
import os
import shutil
import sys
from collections import Counter

import pysam
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/tools")))
import build_labranchor_ism_grch38 as ism  # noqa: E402

pytestmark = pytest.mark.skipif(not shutil.which("sort"), reason="needs coreutils sort")
SEQ = ("ACGTTGCAACGGTA" * 1000)[:12000]   # GRCh38 chr1


class FakeLiftOver:
    """0-based like pyliftover. chr1: +100 below 5000, gap 5000-5009, +105 from 6000; chr2 inverted; chr3 two hits."""
    def convert_coordinate(self, chrom, p0):
        if chrom == "chr1":
            if 5000 <= p0 < 5010:
                return []
            return [("chr1", p0 + (105 if p0 >= 6000 else 100), "+", 1)]
        if chrom == "chr2":
            return [("chr1", 9000 - p0, "-", 1)]
        if chrom == "chr3":
            return [("chr1", p0, "+", 2), ("chr1", p0 + 1, "+", 1)]
        return None


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "ref.fa"
    path.write_text(">chr1\n" + SEQ + "\n")
    pysam.faidx(str(path))
    return pysam.FastaFile(str(path))


def g38(pos1):
    return SEQ[pos1 - 1]


def other(base):
    return "A" if base != "A" else "C"


def test_builder(tmp_path, fasta):
    rows = [
        ("chr1", 3000, g38(3100), "T", "+", 3030, "-0.5", "0.1"),        # clean, listed first: output is sorted
        ("chr1", 1000, g38(1100), "G", "+", 1030, "-0.9", "-0.2"),       # clean
        ("chr1", 1001, other(g38(1101)), "G", "+", 1030, "0", "0"),       # REF differs in GRCh38
        ("chr1", 5005, "A", "G", "+", 4980, "0", "0"),                    # pos in a chain gap
        ("chr1", 5995, g38(6095), "T", "-", 6010, "0", "0"),              # window crosses a chain change
        ("chr2", 100, "A", "G", "+", 130, "0", "0"),                      # inverted mapping
        ("chr3", 100, "A", "G", "+", 130, "0", "0"),                      # multi-mapped
    ]
    src = tmp_path / "ism.hg19.tsv.gz"
    with gzip.open(src, "wt") as fh:
        fh.write("chrom\tpos\tref\talt\tstrand\t3pSS\td_best_bp\td_max_bp\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    out = tmp_path / "ism.grch38.tsv.gz"
    counts = ism.build(str(src), FakeLiftOver(), fasta, str(out))
    assert counts == Counter({"rows": 7, "lifted": 2, "ref_mismatch_grch38": 1, "pos_unmapped": 1,
                              "distance_to_3pSS_changed": 1, "3pSS_minus_strand_mapping": 1, "3pSS_multi_mapped": 1})
    with gzip.open(out, "rt") as fh:
        lines = fh.read().splitlines()
    assert lines[0] == ism.HEADER.rstrip("\n")
    assert [l.split("\t")[:6] for l in lines[1:]] == [
        ["chr1", "1100", g38(1100), "G", "+", "1130"], ["chr1", "3100", g38(3100), "T", "+", "3130"]]
    tbx = pysam.TabixFile(str(out))
    assert [r.split("\t")[1] for r in tbx.fetch("chr1", 3099, 3100)] == ["3100"]       # 1-based index
