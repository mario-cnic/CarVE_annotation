"""vcf_parser_pysam.py (VCF_TO_TSV): every record reaches the TSV, genotype columns hold their own field."""
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("vcf_parser")
SCRIPT = Path(__file__).resolve().parents[2] / "src/python/parsing/vcf_parser_pysam.py"
VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|SYMBOL|Feature">
##INFO=<ID=LaBranchoR_score,Number=1,Type=Float,Description="x">
##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tG\t50\tPASS\tCSQ=G|missense_variant|GENEA|ENST1,G|intron_variant|GENEB|ENST2;LaBranchoR_score=0.5\tGT:GQ:DP:AD\t0/1:60:30:12,18
chr1\t200\t.\tC\tT\t50\tPASS\tLaBranchoR_score=0.7\tGT:GQ:DP:AD\t0/1:60:40:20,20
chr1\t300\t.\tG\tA\t50\tPASS\tCSQ=A|synonymous_variant|GENEA|ENST1\tGT:GQ:DP:AD\t1/1:60:25:0,25
"""


def run(tmp_path, *flags, cols=True):
    (tmp_path / "in.vcf").write_text(VCF)
    (tmp_path / "cols.txt").write_text("Allele\nConsequence\nSYMBOL\nFeature\n")
    args = [sys.executable, str(SCRIPT), "--input", "in.vcf", "--output", "out.tsv", "--overwrite", *flags]
    if cols:
        args += ["--vep_columns", "cols.txt"]
    res = subprocess.run(args, cwd=tmp_path, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr[-2000:]
    lines = (tmp_path / "out.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def test_record_without_csq_keeps_a_row(tmp_path):
    rows = run(tmp_path, "--add_info", "--add_vep")          # VCF_TO_TSV's flags
    assert [(r["POS"], r["Feature"]) for r in rows] == [("100", "ENST1"), ("100", "ENST2"), ("200", "."), ("300", "ENST1")]
    assert rows[2]["LaBranchoR_score"] == "0.7" and rows[2]["SYMBOL"] == "."


def test_gene_filter_still_drops_records_without_csq(tmp_path):
    rows = run(tmp_path, "--add_vep", "--gene_set", "GENEA")
    assert [r["POS"] for r in rows] == ["100", "300"]


def test_an_error_on_a_record_fails_instead_of_dropping_it():
    # MISC-29: the exception was logged and the record skipped, so the run passed with a record missing.
    import io
    sys.path.insert(0, str(SCRIPT.parent))
    import vcf_parser_pysam as vp
    bad = {"variant_id": "chr1_1_A_G", "CHROM": "chr1", "vep_info": {"G": 5}}   # not a list of entries
    with pytest.raises(TypeError):
        vp.write_values(vcf=[bad], file=io.StringIO(), base_header=["CHROM"], vep_cols=["SYMBOL"], info_cols=[],
                        include_gt=False, include_gq=False, include_dp=False, include_ad=False,
                        selected_individuals=[], filter_by="SYMBOL")


def test_depth_column_holds_depth(tmp_path):
    # MISC-30: include_dp / include_ad were passed swapped, so DP_<sample> got the AD field.
    rows = run(tmp_path, "--add_vep", "--add_dp")
    assert [r["DP_S1"] for r in rows] == ["30", "30", "40", "25"]


def test_allele_depth_column_holds_ref_and_alt_depths(tmp_path):
    # MISC-31: AD_<sample> was str(Genotype.allele_depth), a flag of the vcf_parser package: always "False".
    rows = run(tmp_path, "--add_vep", "--add_dp", "--add_ad")
    assert [(r["DP_S1"], r["AD_S1"]) for r in rows] == [("30", "12,18"), ("30", "12,18"), ("40", "20,20"), ("25", "0,25")]
