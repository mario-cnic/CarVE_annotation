import os
import shutil
import sys

import pandas as pd
import pysam
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

from extract_genotypes import parse_gt, run
from genotypes_to_wide import run as wide_run

BCFTOOLS = shutil.which("bcftools") or "/data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools"
pytestmark = pytest.mark.skipif(not os.access(BCFTOOLS, os.X_OK), reason="bcftools not available")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chrX,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="PS">
##FORMAT=<ID=FT,Number=1,Type=String,Description="FT">
"""


def write_vcf(path, samples, records, fmt="GT:GQ:DP:AD:PS:FT", header=HEADER):
    """records: (chrom, pos, ref, alt, [per-sample FORMAT strings])."""
    txt = str(path)[:-3]
    with open(txt, "w") as f:
        f.write(header)
        cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        f.write(cols + (f"\tFORMAT\t" + "\t".join(samples) if samples else "") + "\n")
        for chrom, pos, ref, alt, calls in records:
            line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t."
            if samples:
                line += f"\t{fmt}\t" + "\t".join(calls)
            f.write(line + "\n")
    pysam.tabix_index(txt, preset="vcf", force=True)
    return txt + ".gz"


TRIO = ("PROBAND-1", "FATHER_1", "MOTHER_1")
RECORDS = [
    ("chr1", 100, "A", "C", ["0/1:99:30:15,15:.:PASS", "0/0:90:28:28,0:.:PASS", "1|0:80:25:12,13:100:LowQual"]),
    ("chr1", 200, "G", "T", ["./.:.:.:.:.:.", "1/1:50:10:0,10:.:PASS", "./1:20:5:2,3:.:.", ]),
    ("chr1", 300, "A", "C,G", ["1/2:60:20:0,10,10:.:PASS", "0/1:60:20:10,10,0:.:PASS", "0/0:60:20:20,0,0:.:PASS"]),
    ("chrX", 400, "T", "A", ["1:40:12:0,12:.:PASS", "0:40:12:12,0:.:PASS", ".:.:.:.:.:."]),
]


def go(tmp_path, vcf, name="g.pq"):
    out = str(tmp_path / name)
    rc = run(vcf, out, BCFTOOLS, batch_rows=5)
    return rc, out


@pytest.mark.parametrize("raw,expected", [
    ("0/0", ("hom_ref", 0, 0, False, 2)),
    ("0/1", ("het", 0, 1, False, 2)),
    ("1/0", ("alt_plus_other_allele", 1, 0, False, 2)),
    ("2/0", ("alt_plus_other_allele", 2, 0, False, 2)),
    ("0/2", ("het", 0, 2, False, 2)),
    ("1|0", ("het", 1, 0, True, 2)),
    ("1/1", ("hom_alt", 1, 1, False, 2)),
    ("./.", ("no_call", None, None, False, 2)),
    (".", ("no_call", None, None, False, 1)),
    ("./1", ("partial_no_call", None, 1, False, 2)),
    ("1/2", ("multiallelic_other", 1, 2, False, 2)),
    ("1", ("haploid_alt", 1, None, False, 1)),
    ("0", ("haploid_ref", 0, None, False, 1)),
    ("0/0/1", ("het", 0, 0, False, 3)),
    ("1|1", ("hom_alt", 1, 1, True, 2)),
])
def test_parse_gt(raw, expected):
    assert parse_gt(raw) == expected


def test_trio_long_table_has_one_row_per_variant_and_sample(tmp_path):
    vcf = write_vcf(tmp_path / "trio.vcf.gz", TRIO, RECORDS)
    rc, out = go(tmp_path, vcf)
    assert rc == 0
    df = pd.read_parquet(out)
    assert len(df) == len(RECORDS) * 3
    assert list(df.sample_id.unique()) == list(TRIO)
    assert df.Locus.iloc[0] == "chr1:100-A-C"
    row = df[(df.Locus == "chr1:100-A-C") & (df.sample_id == "MOTHER_1")].iloc[0]
    assert (row.gt_raw, row.gt_status, row.phased, row.gq, row.dp, row.ad_ref, row.ad_alt, row.ps, row.ft) == \
        ("1|0", "het", True, 80, 25, 12, 13, "100", "LowQual")


def test_statuses_and_missing_values_stay_null(tmp_path):
    vcf = write_vcf(tmp_path / "trio.vcf.gz", TRIO, RECORDS)
    _, out = go(tmp_path, vcf)
    df = pd.read_parquet(out).set_index(["Locus", "sample_id"])
    assert df.loc[("chr1:200-G-T", "PROBAND-1"), "gt_status"] == "no_call"
    assert pd.isna(df.loc[("chr1:200-G-T", "PROBAND-1"), "gq"]) and pd.isna(df.loc[("chr1:200-G-T", "PROBAND-1"), "ad_ref"])
    assert df.loc[("chr1:200-G-T", "MOTHER_1"), "gt_status"] == "partial_no_call"
    assert df.loc[("chr1:300-A-C,G", "PROBAND-1"), "gt_status"] == "multiallelic_other"
    assert df.loc[("chrX:400-T-A", "PROBAND-1"), "gt_status"] == "haploid_alt"
    assert df.loc[("chrX:400-T-A", "FATHER_1"), "gt_status"] == "haploid_ref"
    assert df.loc[("chrX:400-T-A", "MOTHER_1"), "gt_status"] == "no_call"
    assert df.loc[("chr1:100-A-C", "FATHER_1"), "gt_status"] == "hom_ref"
    assert pd.isna(df.loc[("chr1:100-A-C", "PROBAND-1"), "ps"])
    assert pd.isna(df.loc[("chr1:200-G-T", "PROBAND-1"), "ad"])


def test_single_sample_without_optional_format_fields(tmp_path):
    header = HEADER.split("##FORMAT=<ID=PS")[0]
    header = "\n".join(l for l in header.splitlines() if "ID=AD" not in l) + "\n"
    vcf = write_vcf(tmp_path / "one.vcf.gz", ("S1",), [("chr1", 100, "A", "C", ["0/1:50:10"])], fmt="GT:GQ:DP", header=header)
    rc, out = go(tmp_path, vcf)
    df = pd.read_parquet(out)
    assert rc == 0 and len(df) == 1 and df.gq[0] == 50 and df.ad.isna().all()
    assert "ps" not in df.columns and "ft" not in df.columns


def test_sites_only_vcf_writes_no_file(tmp_path):
    vcf = write_vcf(tmp_path / "sites.vcf.gz", (), [("chr1", 100, "A", "C", [])])
    rc, out = go(tmp_path, vcf)
    assert rc == 0 and not os.path.exists(out)


def test_samples_without_gt_field_writes_no_file(tmp_path):
    header = "\n".join(l for l in HEADER.splitlines() if "ID=GT" not in l) + "\n"
    vcf = write_vcf(tmp_path / "nogt.vcf.gz", ("S1",), [("chr1", 100, "A", "C", ["50"])], fmt="GQ", header=header)
    rc, out = go(tmp_path, vcf)
    assert rc == 0 and not os.path.exists(out)


def make_main_table(tmp_path, loci, rows_per_locus=3):
    rows = [{"Locus": l, "Feature": f"ENST{i}"} for l in loci for i in range(rows_per_locus)]
    path = tmp_path / "main.pq"
    pd.DataFrame(rows).to_parquet(path)
    return str(path)


def test_wide_view_maps_statuses_and_keeps_row_count(tmp_path):
    vcf = write_vcf(tmp_path / "trio.vcf.gz", TRIO, RECORDS)
    _, gt = go(tmp_path, vcf)
    main = make_main_table(tmp_path, [r[0] + ":" + str(r[1]) + "-" + r[2] + "-" + r[3] for r in RECORDS])
    out = str(tmp_path / "wide.pq")
    assert wide_run(main, gt, out, max_samples=10) == 0
    w = pd.read_parquet(out)
    assert len(w) == len(RECORDS) * 3
    first = w[w.Locus == "chr1:100-A-C"].iloc[0]
    assert (first["GT_PROBAND-1"], first["GT_FATHER_1"], first["GT_MOTHER_1"]) == ("HET", "HOMREF", "HET")  # 1|0 phased het
    assert first["GQ_PROBAND-1"] == 99
    second = w[w.Locus == "chr1:200-G-T"].iloc[0]
    assert (second["GT_PROBAND-1"], second["GT_FATHER_1"], second["GT_MOTHER_1"]) == ("MISSING", "HOMALT", "MISSING")
    assert w[w.Locus == "chr1:300-A-C,G"].iloc[0]["GT_PROBAND-1"] == "MISSING"
    assert w[w.Locus == "chrX:400-T-A"].iloc[0]["GT_PROBAND-1"] == "HOMALT"
    assert pd.read_parquet(main).shape[1] == 2   # main table untouched


def test_wide_view_refuses_too_many_samples(tmp_path):
    vcf = write_vcf(tmp_path / "trio.vcf.gz", TRIO, RECORDS)
    _, gt = go(tmp_path, vcf)
    main = make_main_table(tmp_path, ["chr1:100-A-C"])
    assert wide_run(main, gt, str(tmp_path / "w.pq"), max_samples=2) == 1
    assert not os.path.exists(tmp_path / "w.pq")


def test_split_multiallelic_1_0_is_not_reported_as_het(tmp_path):
    recs = [("chr1", 100, "A", "C", ["1/0:60:20:0,12:.:PASS", "0/1:60:20:10,10:.:PASS", "1|0:60:20:10,10:.:PASS"])]
    vcf = write_vcf(tmp_path / "split.vcf.gz", TRIO, recs)
    _, out = go(tmp_path, vcf)
    df = pd.read_parquet(out).set_index("sample_id")
    assert df.loc["PROBAND-1", "gt_status"] == "alt_plus_other_allele"
    assert df.loc["FATHER_1", "gt_status"] == "het" and df.loc["MOTHER_1", "gt_status"] == "het"


def test_wide_view_batches_give_the_same_result(tmp_path):
    vcf = write_vcf(tmp_path / "trio.vcf.gz", TRIO, RECORDS)
    _, gt = go(tmp_path, vcf)
    main = make_main_table(tmp_path, [r[0] + ":" + str(r[1]) + "-" + r[2] + "-" + r[3] for r in RECORDS])
    a, b = str(tmp_path / "a.pq"), str(tmp_path / "b.pq")
    assert wide_run(main, gt, a, 10, batch_rows=1000) == 0 and wide_run(main, gt, b, 10, batch_rows=2) == 0
    pd.testing.assert_frame_equal(pd.read_parquet(a), pd.read_parquet(b))


def test_wide_view_fails_when_a_table_row_has_no_genotype(tmp_path):
    vcf = write_vcf(tmp_path / "trio.vcf.gz", TRIO, RECORDS)
    _, gt = go(tmp_path, vcf)
    main = make_main_table(tmp_path, ["chr1:100-A-C", "chr9:999-A-T"])
    assert wide_run(main, gt, str(tmp_path / "w.pq"), 10) == 1


def test_wide_view_rejects_non_parquet_table(tmp_path):
    vcf = write_vcf(tmp_path / "trio.vcf.gz", TRIO, RECORDS)
    _, gt = go(tmp_path, vcf)
    (tmp_path / "t.tsv").write_text("Locus\nchr1:100-A-C\n")
    assert wide_run(str(tmp_path / "t.tsv"), gt, str(tmp_path / "w.pq"), 10) == 1
