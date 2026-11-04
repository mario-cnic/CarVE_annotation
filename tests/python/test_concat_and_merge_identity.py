"""CONCAT_CHUNKS and MERGE_ANNOTATIONS check CHROM/POS/REF/ALT identity, not only counts (MISC-6).

The real script blocks are rendered as Nextflow interpolates them and run with bcftools on fixtures.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src/python"))
import annotation_cache as ac  # noqa: E402

BCF, TABIX = shutil.which("bcftools"), shutil.which("tabix")
pytestmark = pytest.mark.skipif(not (BCF and TABIX), reason="bcftools/tabix not available")

HEADER = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
          '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="x">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')


def vcf(path, records, info="."):
    path.parent.mkdir(parents=True, exist_ok=True)
    plain = path.with_suffix("")
    plain.write_text(HEADER + "".join(f"chr1\t{p}\t.\t{r}\t{a}\t.\t.\t{info}\n" for p, r, a in records))
    subprocess.run([BCF, "view", str(plain), "-Oz", "-o", str(path)], check=True)
    subprocess.run([TABIX, "-f", "-p", "vcf", str(path)], check=True)
    return path


def render(nf: str, process: str, values: dict) -> str:
    body = ac.script_block((REPO / "modules/local" / nf).read_text(), process).split('"""')[1]
    for k, v in values.items():
        body = body.replace("${" + k + "}", v)
    assert "${" not in body, body
    return body.replace("\\\\", "\0").replace("\\$", "$").replace("\0", "\\")


def concat(tmp_path, chunks):
    paths = [vcf(tmp_path / f"chunk_{i:02d}" / "c.vcf.gz", recs) for i, recs in enumerate(chunks, 1)]
    vcf(tmp_path / "orig.vcf.gz", [(100, "A", "G"), (200, "C", "T"), (300, "G", "A"), (400, "T", "C")])
    rel = [str(p.relative_to(tmp_path)) for p in paths]
    script = render("chunk.nf", "CONCAT_CHUNKS", {
        "params.bcftools": BCF, "params.tabix": TABIX, "chunk_list[0]": rel[0], "chunk_list.join(' ')": " ".join(rel),
        "final_vcf": "run.annX.vcf.gz", "orig_vcf": "orig.vcf.gz", "suffix": "annX"})
    return subprocess.run(["bash", "-ue", "-c", script], cwd=tmp_path, capture_output=True, text=True)


def test_concat_passes_when_sites_are_identical(tmp_path):
    res = concat(tmp_path, [[(300, "G", "A"), (400, "T", "C")], [(100, "A", "G"), (200, "C", "T")]])
    assert res.returncode == 0, res.stderr
    assert "CHROM/POS/REF/ALT identical to the input" in res.stdout


def test_concat_fails_when_a_site_changed_but_the_count_did_not(tmp_path):
    res = concat(tmp_path, [[(100, "A", "G"), (200, "C", "T")], [(300, "G", "T"), (400, "T", "C")]])
    assert res.returncode == 1 and "record count parity confirmed" in res.stdout
    assert "CHROM/POS/REF/ALT of the output differ from the input" in res.stderr


def merge(tmp_path, vep_records):
    vep = vcf(tmp_path / "vep.vcf.gz", vep_records)
    files = {}
    for name in ("bp", "pango", "sai", "spip"):
        files[name] = vcf(tmp_path / f"{name}.vcf.gz", [(100, "A", "G")], info="SpliceAI=x")
    values = {"params.bcftools": BCF, "params.tabix": TABIX, "vep_vcf": vep.name, "meta.partition_id": "run",
              **{f"{k}_vcf": v.name for k, v in files.items()}}
    script = render("merge.nf", "MERGE_ANNOTATIONS", values)
    return subprocess.run(["bash", "-ue", "-c", script], cwd=tmp_path, capture_output=True, text=True)


def _annotate_writes_index():
    usage = subprocess.run([BCF, "annotate"], capture_output=True, text=True).stderr
    return "--write-index" in usage


@pytest.mark.skipif(not (BCF and _annotate_writes_index()), reason="this bcftools annotate has no --write-index (-W); the cluster versions have it")
def test_merge_keeps_veps_sites(tmp_path):
    res = merge(tmp_path, [(100, "A", "G"), (200, "C", "T")])
    assert res.returncode == 0, res.stderr
    out = subprocess.run([BCF, "query", "-f", "%POS %INFO/SpliceAI\n", str(tmp_path / "run.annotated.vcf.gz")],
                         capture_output=True, text=True, check=True).stdout.split("\n")
    assert out[:2] == ["100 x", "200 ."]
