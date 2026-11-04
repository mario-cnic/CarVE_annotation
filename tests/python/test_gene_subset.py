"""GENE_SUBSET's real script block (modules/local/gene_subset.nf), rendered as Nextflow would and run
with bcftools: records outside the gene BED are removed (MISC-9), symbolic / breakend / `*` records are
left out (PRED-17), and the subset follows the header's contig order (MISC-20)."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src/python"))
import annotation_cache as ac  # noqa: E402

TOOLS = {t: shutil.which(t) for t in ("bcftools", "tabix", "bgzip")}
pytestmark = pytest.mark.skipif(not all(TOOLS.values()), reason="bcftools/tabix/bgzip not available")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##contig=<ID=chr2,length=100000>
##contig=<ID=chr10,length=100000>
##ALT=<ID=CNV,Description="copy number">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
RECORDS = [  # (chrom, pos, id, ref, alt)
    ("chr1", 150, "in_snv", "A", "G"),
    ("chr1", 160, "in_del", "AC", "A"),
    ("chr1", 170, "in_cnv", "G", "<CNV>"),
    ("chr1", 180, "in_bnd", "G", "G]chr2:123]"),
    ("chr1", 190, "in_star", "G", "*"),
    ("chr1", 5000, "out_snv", "C", "T"),          # outside the BED
    ("chr2", 250, "in_chr2", "T", "C"),
    ("chr10", 350, "in_chr10", "G", "A"),
]
BED = "chr1\t100\t200\nchr10\t300\t400\nchr2\t200\t300\n"   # lexicographic, as resources/v5_genes_loc.bed


def render(script: str, **values) -> str:
    for k, v in values.items():
        script = script.replace("${" + k + "}", v)
    return script.replace("\\$", "$").replace("\\\\\n", "\\\n")


def test_gene_subset_script(tmp_path):
    vcf = tmp_path / "in.vcf"
    vcf.write_text(HEADER + "".join(f"{c}\t{p}\t{i}\t{r}\t{a}\t.\t.\t.\n" for c, p, i, r, a in RECORDS))
    subprocess.run([TOOLS["bgzip"], "-f", str(vcf)], check=True)
    subprocess.run([TOOLS["tabix"], "-p", "vcf", f"{vcf}.gz"], check=True)
    (tmp_path / "genes.bed").write_text(BED)
    block = ac.script_block((REPO / "modules/local/gene_subset.nf").read_text(), "GENE_SUBSET")
    body = block.split('"""')[1]
    script = render(body, **{"params.bcftools": TOOLS["bcftools"], "params.tabix": TOOLS["tabix"],
                             "gene_bed": "genes.bed", "vcf": "in.vcf.gz"})
    res = subprocess.run(["bash", "-ue", "-c", script], cwd=tmp_path, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert "3 record(s) in the gene regions left out" in res.stderr
    ids = subprocess.run([TOOLS["bcftools"], "query", "-f", "%ID\n", str(tmp_path / "subset.vcf.gz")],
                         capture_output=True, text=True, check=True).stdout.split()
    assert ids == ["in_snv", "in_del", "in_chr2", "in_chr10"]      # header order: chr1, chr2, chr10
    assert (tmp_path / "subset.vcf.gz.tbi").is_file()
