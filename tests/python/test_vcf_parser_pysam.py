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
