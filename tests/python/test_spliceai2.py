"""annotate_spliceai2.py lookup and parse_spliceai2 status (PR-03)."""
import os
import sys

import numpy as np
import pandas as pd
import pysam
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, os.path.join(ROOT, "src/python"))
import annotate_spliceai2 as sa2  # noqa: E402

sys.path.insert(0, os.path.join(ROOT, "src/python/parsing"))
from filter_variants import parse_spliceai2  # noqa: E402

# Real layout: chrom pos ref alt strand, then the source columns the lookup reads, then the summary score.
COLS = ["chrom", "pos", "ref", "alt", "strand", "donor_gain_delta_score_0", "donor_gain_dist_0",
        "donor_loss_delta_score_0", "donor_loss_dist_0", "acceptor_gain_delta_score_0", "acceptor_gain_dist_0",
        "acceptor_loss_delta_score_0", "acceptor_loss_dist_0", "spliceai2_summary_score"]
SNV = [("chr1", 100, "A", "G", "+", 0.0, -5, 0.0, 7, 0.81, 1, 0.57, 5, 0.81),
       ("chr1", 200, "C", "T", "+", 0.0, -5, 0.0, 7, 0.0, 1, 0.0, 5, 0.0),   # genes on both strands:
       ("chr1", 200, "C", "T", "-", 0.02, 9, 0.0, 7, 0.0, 1, 0.0, 5, 0.02)]  # the higher row wins
INDEL = [("chr1", 300, "GT", "G", "+", 0.1, 3, 0.4, -2, 0.0, 1, 0.0, 5, 0.4)]


def write_table(path, rows):
    plain = str(path)[:-3]
    with open(plain, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    pysam.tabix_compress(plain, str(path), force=True)
    pysam.tabix_index(str(path), seq_col=0, start_col=1, end_col=1, line_skip=1, force=True)


@pytest.fixture
def scores(tmp_path):
    for kind, rows in (("snv", SNV), ("indel", INDEL)):
        (tmp_path / kind).mkdir()
        write_table(tmp_path / kind / "chr1.tsv.gz", rows)
    return str(tmp_path)


def test_lookup_exact_allele_and_best_strand(tmp_path, scores):
    vcf, out = tmp_path / "in.vcf", tmp_path / "out.vcf"
    header = pysam.VariantHeader()
    for c in ("chr1", "chr2"):
        header.add_line(f"##contig=<ID={c},length=100000>")
    records = [("chr1", 100, ("A", "G")), ("chr1", 100, ("A", "C")), ("chr1", 200, ("C", "T")),
               ("chr1", 300, ("GT", "G")), ("chr1", 300, ("G", "GA")), ("chr2", 100, ("A", "G")),
               ("chr1", 100, ("A", "G", "T"))]
    with pysam.VariantFile(str(vcf), "w", header=header) as fh:
        for chrom, pos, alleles in records:
            fh.write(fh.new_record(contig=chrom, start=pos - 1, alleles=alleles))
    sa2.annotate_spliceai2_vcf(str(vcf), str(out), scores)
    with pysam.VariantFile(str(out)) as fh:
        got = [(r.info.get("SpliceAI2_score"), r.info.get("SpliceAI2_strand"), r.info.get("SpliceAI2_DP_AG")) for r in fh]
    assert got[0] == (pytest.approx(0.81), "+", 1)
    assert got[1] == (None, None, None)               # other ALT at a scored position
    assert got[2][:2] == (pytest.approx(0.02), "-")
    assert got[3][0] == pytest.approx(0.4)            # indel table
    assert got[4] == (None, None, None)               # indel not in the population set
    assert got[5] == (None, None, None)               # chromosome without a file
    assert got[6] == (None, None, None)               # multiallelic record


def test_missing_scores_dir_fails(tmp_path):
    with pytest.raises(FileNotFoundError):
        sa2.SpliceAI2Lookup(str(tmp_path / "nope"))


def test_parse_status_never_zero_for_not_covered():
    df = pd.DataFrame({"SpliceAI2_score": ["0.81", ".", np.nan, "0"], "SpliceAI2_DP_AG": ["1", ".", np.nan, "-3"]})
    res = parse_spliceai2(df.copy())
    assert res["SpliceAI2_status"].tolist() == ["scored", "not_covered", "not_covered", "scored"]
    assert res["SpliceAI2_score"].iloc[3] == 0.0 and np.isnan(res["SpliceAI2_score"].iloc[1])


def test_parse_skips_without_columns():
    res = parse_spliceai2(pd.DataFrame({"Locus": ["chr1:1-A-G"]}))
    assert list(res.columns) == ["Locus"]
