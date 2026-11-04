"""annotate_spliceai.py: a SpliceAI exception fails the run instead of leaving the record unscored (PRED-16).

spliceai/tensorflow are not needed: a stub `spliceai.utils` is put in sys.modules before import.
"""
import os
import sys
import types

import pysam
import pytest

STUB = types.ModuleType("spliceai.utils")
STUB.Annotator = lambda fasta, assembly: None
STUB.get_delta_scores = lambda *a: []  # replaced per test
sys.modules.setdefault("spliceai", types.ModuleType("spliceai"))
sys.modules["spliceai.utils"] = STUB
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python")))

import annotate_spliceai  # noqa: E402


def write_vcf(path, n=3):
    header = pysam.VariantHeader()
    header.add_line("##contig=<ID=chr1,length=100000>")
    with pysam.VariantFile(path, "w", header=header) as out:
        for i in range(n):
            out.write(out.new_record(contig="chr1", start=1000 + 10 * i, alleles=("A", "G")))


def test_scores_and_skips_are_written(tmp_path, monkeypatch):
    src, dst = str(tmp_path / "in.vcf"), str(tmp_path / "out.vcf")
    write_vcf(src)
    # Second record: skipped by spliceai on purpose (returns []), so written without the field.
    scores = {1001: ["G|GENE|0.10|0.00|0.00|0.00|1|2|3|4"], 1011: [], 1021: ["G|GENE|0.50|0.00|0.00|0.00|1|2|3|4"]}
    monkeypatch.setattr(annotate_spliceai, "get_delta_scores", lambda rec, *a: scores[rec.pos])
    annotate_spliceai.annotate_vcf(src, dst, "ref.fa")
    with pysam.VariantFile(dst) as vf:
        got = {r.pos: r.info.get("SpliceAI") for r in vf}
    assert got == {1001: ("G|GENE|0.10|0.00|0.00|0.00|1|2|3|4",), 1011: None,
                   1021: ("G|GENE|0.50|0.00|0.00|0.00|1|2|3|4",)}


def test_exception_fails_and_names_the_record(tmp_path, monkeypatch, capsys):
    src, dst = str(tmp_path / "in.vcf"), str(tmp_path / "out.vcf")
    write_vcf(src)

    def boom(rec, *a):
        if rec.pos == 1011:
            raise MemoryError("OOM")
        return []
    monkeypatch.setattr(annotate_spliceai, "get_delta_scores", boom)
    with pytest.raises(MemoryError):
        annotate_spliceai.annotate_vcf(src, dst, "ref.fa")
    assert "chr1:1011 A>G" in capsys.readouterr().err
