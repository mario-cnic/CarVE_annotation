"""The provenance stamp in the Parquet file metadata of the three published table writers."""
import json
import os
import shutil
import sys

import pandas as pd
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import write_run_manifest as wrm
from extract_genotypes import run as extract_run
from genotypes_to_wide import run as wide_run
from provenance_stamp import STAMP_KEY, load_stamp, with_stamp
from tag_transcript_priority import write_table
from test_extract_genotypes import HEADER, RECORDS, TRIO, write_vcf

BCFTOOLS = shutil.which("bcftools") or "/data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools"
needs_bcftools = pytest.mark.skipif(not os.access(BCFTOOLS, os.X_OK), reason="bcftools not available")

STAMP = {"kind": "table_stamp", "run_id": "r1", "databases": {"pangolin_db": {"release": "GENCODE 45"}}}


@pytest.fixture
def stamp_file(tmp_path):
    path = tmp_path / "TABLE_STAMP.json"
    path.write_text(json.dumps(STAMP, indent=2))
    return str(path)


def footer(path):
    return {k.decode(): v.decode() for k, v in (pq.ParquetFile(path).metadata.metadata or {}).items()}


def test_load_stamp_is_compact_json_and_none_without_path(stamp_file):
    assert load_stamp(None) is None
    assert json.loads(load_stamp(stamp_file)) == STAMP
    assert "\n" not in load_stamp(stamp_file)


def test_unreadable_stamp_stops_the_writer(tmp_path):
    with pytest.raises(SystemExit):
        load_stamp(str(tmp_path / "missing.json"))
    bad = tmp_path / "bad.json"
    bad.write_text("not json")
    with pytest.raises(SystemExit):
        load_stamp(str(bad))


def test_with_stamp_keeps_existing_entries():
    import pyarrow as pa
    schema = pa.schema([("a", pa.int64())], metadata={b"pandas": b"{}"})
    out = with_stamp(schema, "{}")
    assert out.metadata[b"pandas"] == b"{}" and out.metadata[STAMP_KEY.encode()] == b"{}"
    assert with_stamp(schema, None) is schema


def test_tag_writer_stamps_and_keeps_pandas_metadata_and_data(tmp_path, stamp_file):
    df = pd.DataFrame({"Locus": ["chr1:1-A-C", "chr1:2-G-T"], "x": [1.5, None], "TRANSCRIPT_PRIORITY_TIER": [1, 3]})
    stamped, plain, reference = (str(tmp_path / n) for n in ("stamped.pq", "plain.pq", "reference.pq"))
    write_table(df, stamped, load_stamp(stamp_file))
    write_table(df, plain)
    df.to_parquet(reference, index=False)
    meta = footer(stamped)
    assert json.loads(meta[STAMP_KEY]) == STAMP and "pandas" in meta and "ARROW:schema" in meta
    assert STAMP_KEY not in footer(plain)
    pd.testing.assert_frame_equal(pd.read_parquet(stamped), pd.read_parquet(reference))
    pd.testing.assert_frame_equal(pd.read_parquet(plain), pd.read_parquet(reference))


@needs_bcftools
def test_genotype_table_is_stamped(tmp_path, stamp_file):
    vcf = write_vcf(tmp_path / "in.vcf.gz", TRIO, RECORDS)
    out = str(tmp_path / "g.pq")
    assert extract_run(vcf, out, BCFTOOLS, 500000, load_stamp(stamp_file)) == 0
    assert json.loads(footer(out)[STAMP_KEY]) == STAMP
    out_plain = str(tmp_path / "g_plain.pq")
    assert extract_run(vcf, out_plain, BCFTOOLS, 500000) == 0
    assert STAMP_KEY not in footer(out_plain)


@needs_bcftools
def test_genotype_table_without_records_is_stamped(tmp_path, stamp_file):
    vcf = write_vcf(tmp_path / "empty.vcf.gz", TRIO, [])
    out = str(tmp_path / "g.pq")
    assert extract_run(vcf, out, BCFTOOLS, 500000, load_stamp(stamp_file)) == 0
    assert json.loads(footer(out)[STAMP_KEY]) == STAMP


@needs_bcftools
def test_wide_table_is_stamped(tmp_path, stamp_file):
    vcf = write_vcf(tmp_path / "in.vcf.gz", TRIO, RECORDS)
    genotypes = str(tmp_path / "g.pq")
    assert extract_run(vcf, genotypes, BCFTOOLS, 500000) == 0
    loci = pd.read_parquet(genotypes, columns=["Locus"])["Locus"].drop_duplicates()
    table = str(tmp_path / "main.pq")
    write_table(pd.DataFrame({"Locus": loci, "x": range(len(loci))}), table)
    out = str(tmp_path / "wide.pq")
    assert wide_run(table, genotypes, out, 10, stamp=load_stamp(stamp_file)) == 0
    meta = footer(out)
    assert json.loads(meta[STAMP_KEY]) == STAMP and "pandas" in meta
    assert any(c.startswith("GT_") for c in pq.ParquetFile(out).schema_arrow.names)


def test_compact_provenance_holds_the_identity_fields():
    launch = {"generated_at": "t", "cli_flags": {"run_id": "r1"}, "pangolin_distance": "4999",
              "annotation_pipeline_repo": {"commit": "abc", "dirty": False},
              "databases": {"pangolin_db": {"release": "GENCODE 45"}},
              "external_provenance": {"resources": {"pangolin_db": {"path": "/db", "sha256": "ff"}},
                                      "envs": {"pangolin_python": {"binary": "/py", "fingerprint_sha256": "aa"}}}}
    out = wrm.compact_provenance(launch)
    assert out["run_id"] == "r1" and out["code"]["commit"] == "abc"
    assert out["resources"]["pangolin_db"]["sha256"] == "ff"
    assert out["databases"]["pangolin_db"]["release"] == "GENCODE 45"
    assert out["tools"]["envs"]["pangolin_python"]["fingerprint_sha256"] == "aa"
    assert out["parameters"]["pangolin_distance"] == "4999"


def test_launch_writes_the_stamp_file(tmp_path, monkeypatch):
    import subprocess
    import argparse
    repo = tmp_path / "repo"
    repo.mkdir()
    for cmd in (["init", "-q"], ["add", "-A"]):
        subprocess.run(["git", *cmd], cwd=repo, check=True, capture_output=True)
    (repo / "main.nf").write_text("workflow { }\n")
    subprocess.run(["git", "add", "-A"], cwd=repo, check=True, capture_output=True)
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "i"], cwd=repo,
                   check=True, capture_output=True)
    monkeypatch.setattr(wrm, "REPO_ROOT", repo)
    out = repo / "nf_work" / "out"
    args = argparse.Namespace(
        outdir=str(out), command_line="nextflow run main.nf", nextflow_binary="/bin/true", nextflow_main="main.nf",
        config="nextflow.config", input=str(repo / "main.nf"), run_id="r1", profile="local_dev", spip_tier=None,
        output_format=None, resume=False, bcftools=None, vep_sif=None, vep_cache_version=None,
        spliceai_distance=None, pangolin_distance="4999")
    wrm.cmd_launch(args)
    stamp = json.loads((out / "pipeline_info" / "TABLE_STAMP.json").read_text())
    assert stamp["kind"] == "table_stamp" and stamp["run_id"] == "r1"
    assert stamp["parameters"]["pangolin_distance"] == "4999"
