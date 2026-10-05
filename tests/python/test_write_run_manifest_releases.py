"""Database release readers and the per-table provenance file in write_run_manifest.py.

Synthetic fixtures only.
"""
import argparse
import json
import os
import shutil
import sqlite3
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import write_run_manifest as wrm


def make_pangolin_db(path, description):
    con = sqlite3.connect(path)
    con.execute("create table directives (directive text)")
    con.executemany("insert into directives values (?)",
                    [("description: " + description,), ("provider: GENCODE",), ("date: 2023-09-19",)])
    con.commit()
    con.close()


def test_pangolin_release_from_directives_and_build_manifest(tmp_path):
    db = tmp_path / "gencode.v45.test.grch38.db"
    make_pangolin_db(db, "evidence-based annotation of the human genome (GRCh38), version 45 (Ensembl 111)")
    (tmp_path / "gencode.v45.test.grch38.BUILD_MANIFEST.txt").write_text("output_db_sha256: abc\n")
    info = wrm.read_pangolin_db(db)
    assert info["release"].endswith("version 45 (Ensembl 111)")
    assert "output_db_sha256: abc" in info["build_manifest"]


def test_pangolin_unreadable_db_gives_no_release_and_a_reason(tmp_path):
    bad = tmp_path / "x.db"
    bad.write_text("not sqlite")
    info = wrm.read_pangolin_db(bad)
    assert info["release"] is None and info["reason"]


def test_vep_cache_sources(tmp_path):
    cache = tmp_path / "cache" / "homo_sapiens" / "111_GRCh38"
    cache.mkdir(parents=True)
    (cache / "info.txt").write_text("source_ClinVar\t202306\nsource_dbSNP\t156\nsift\tb\nspecies\thomo_sapiens\n")
    info = wrm.read_vep_cache(str(tmp_path), "111")
    assert info["values"] == {"source_ClinVar": "202306", "source_dbSNP": "156", "species": "homo_sapiens"}
    assert wrm.read_vep_cache(str(tmp_path), "999")["release"] is None


def test_dbnsfp_release_and_clinvar_release(tmp_path):
    (tmp_path / "dbNSFP9.9.readme.txt").write_text(
        "dbNSFP version 9.9\n\nRelease:\n\tJune 3, 2024\n\nMajor sources:\n\t\tclinvar release 20240215 (hg38)\n")
    info = wrm.read_dbnsfp(tmp_path / "dbNSFP9.9_grch38.gz")
    assert info["release"] == "dbNSFP version 9.9"
    assert info["release_date"] == "June 3, 2024"
    assert info["clinvar_release"] == "20240215"


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not available")
def test_vcf_header_release_and_missing_release(tmp_path):
    with_date = tmp_path / "a.vcf"
    with_date.write_text("##fileformat=VCFv4.2\n##fileDate=2026-09-28\n##source=ClinVar\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    without = tmp_path / "b.vcf"
    without.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    bcftools = shutil.which("bcftools")
    assert wrm.read_vcf_header(with_date, bcftools)["release"] == "2026-09-28"
    nodate = wrm.read_vcf_header(without, bcftools)
    assert nodate["release"] is None and nodate["reason"]


def test_name_hint_is_separate_from_release():
    assert "v4.1" in wrm.name_hint("/x/gnomAD.v4.1.vcf.gz")
    assert "20260928" in wrm.name_hint("/x/clinvar_20260928.vcf.gz")
    assert wrm.name_hint("/x/plain.txt") == []


def test_database_releases_never_raises_on_missing_files(tmp_path):
    resources = {"pangolin_db": wrm.resource_identity(str(tmp_path / "missing.db")),
                 "clinvar_vcf": wrm.resource_identity(str(tmp_path / "missing.vcf.gz"))}
    out = wrm.database_releases(resources, None, None, None)
    assert out["pangolin_db"]["release"] is None and out["clinvar_vcf"]["release"] is None


def make_manifest(tmp_path, exit_code=0):
    out = tmp_path / "run"
    (out / "pipeline_info").mkdir(parents=True)
    (out / "tables").mkdir()
    (out / "genotypes").mkdir()
    (out / "tables" / "r1.parsed.clean.pq").write_bytes(b"PAR1")
    (out / "tables" / "r1.parsed.clean.withGT.pq").write_bytes(b"PAR1")
    (out / "genotypes" / "r1.genotypes.pq").write_bytes(b"PAR1")
    launch = {"generated_at": "2026-10-05T00:00:00+00:00", "cli_flags": {"run_id": "r1"},
              "pangolin_distance": "4999", "annotation_pipeline_repo": {"commit": "abc", "dirty": False},
              "databases": {"pangolin_db": {"release": "GENCODE 45"}},
              "tool_versions": {"vep_image": {"vep_version": "111.0"}},
              "external_provenance": {"resources": {"pangolin_db": {"path": "/db", "sha256": "ff", "size_bytes": 1}},
                                      "envs": {"pangolin_python": {"binary": "/py", "fingerprint_sha256": "aa",
                                                                   "key_packages": {"pangolin": "pip:pangolin-1.0.2"}}}}}
    (out / "pipeline_info" / "RUN_MANIFEST.json").write_text(json.dumps({"launch": launch, "completed": None}))
    return out


def test_completed_writes_a_provenance_file_next_to_each_table(tmp_path):
    out = make_manifest(tmp_path)
    wrm.cmd_completed(argparse.Namespace(outdir=str(out), exit_code=0, nextflow_log=None))
    for table in ("tables/r1.parsed.clean.pq", "tables/r1.parsed.clean.withGT.pq", "genotypes/r1.genotypes.pq"):
        side = json.loads((out / (table + ".provenance.json")).read_text())
        assert side["kind"] == "table_provenance"
        assert side["table"]["name"] == os.path.basename(table)
        assert side["run"]["run_id"] == "r1" and len(side["run"]["full_manifest_sha256"]) == 64
        assert side["databases"]["pangolin_db"]["release"] == "GENCODE 45"
        assert side["tools"]["vep_image"]["vep_version"] == "111.0"
        assert side["tools"]["envs"]["pangolin_python"]["key_packages"]["pangolin"] == "pip:pangolin-1.0.2"
        assert side["parameters"]["pangolin_distance"] == "4999"


def test_failed_run_writes_no_provenance_file(tmp_path):
    out = make_manifest(tmp_path)
    wrm.cmd_completed(argparse.Namespace(outdir=str(out), exit_code=1, nextflow_log=None))
    assert not list(out.rglob("*.provenance.json"))


def test_existing_provenance_file_is_not_overwritten(tmp_path):
    out = make_manifest(tmp_path)
    side = out / "tables" / "r1.parsed.clean.pq.provenance.json"
    side.write_text("{\"keep\": true}")
    wrm.cmd_completed(argparse.Namespace(outdir=str(out), exit_code=0, nextflow_log=None))
    assert json.loads(side.read_text()) == {"keep": True}
