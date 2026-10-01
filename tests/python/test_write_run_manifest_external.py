"""External-provenance capture in write_run_manifest.py: resources, code repos, conda envs.

Synthetic fixtures only; nothing here touches the real shared/utils or cluster data.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import write_run_manifest as wrm


def git(repo, *args):
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *args],
                          cwd=repo, capture_output=True, check=True, text=True).stdout


def make_env(root, packages):
    (root / "bin").mkdir(parents=True)
    (root / "bin" / "python3").write_text("")
    meta = root / "conda-meta"
    meta.mkdir()
    for p in packages:
        (meta / f"{p}.json").write_text("{}")
    return root / "bin" / "python3"


def test_split_named_requires_name_and_value():
    assert wrm.split_named("fasta=/a/b=c") == ("fasta", "/a/b=c")
    for bad in ("nopath", "=x", "x="):
        with pytest.raises(argparse.ArgumentTypeError):
            wrm.split_named(bad)


def test_resource_small_file_is_hashed(tmp_path):
    f = tmp_path / "db.bin"
    f.write_bytes(b"abc")
    r = wrm.resource_identity(str(f))
    assert r["exists"] and r["size_bytes"] == 3
    assert r["sha256"] == hashlib.sha256(b"abc").hexdigest()
    assert r["hash_skipped"] is None


def test_resource_missing_is_not_guessed(tmp_path):
    r = wrm.resource_identity(str(tmp_path / "absent"))
    assert r["exists"] is False
    assert r["sha256"] is None and r["size_bytes"] is None


def test_resource_over_limit_skips_hash_and_records_fai(tmp_path, monkeypatch):
    monkeypatch.setattr(wrm, "RESOURCE_HASH_MAX_BYTES", 10)
    f = tmp_path / "genome.fasta"
    f.write_bytes(b"A" * 100)
    (tmp_path / "genome.fasta.fai").write_text("chr1\t100\t6\t60\t61\n")
    r = wrm.resource_identity(str(f))
    assert r["sha256"] is None
    assert "size >" in r["hash_skipped"]
    assert r["sidecars"][".fai"]["sha256"] == hashlib.sha256(b"chr1\t100\t6\t60\t61\n").hexdigest()


def test_resource_directory_records_no_hash(tmp_path):
    r = wrm.resource_identity(str(tmp_path))
    assert r["is_dir"] and r["sha256"] is None and r["size_bytes"] is None


def test_code_repo_commit_and_tracked_dirty(tmp_path):
    repo = tmp_path / "shared"
    repo.mkdir()
    git(repo, "init", "-q")
    (repo / "parser.py").write_text("v1\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "init")
    head = git(repo, "rev-parse", "HEAD").strip()

    clean = wrm.code_repo_state(str(repo))
    assert clean["commit"] == head and clean["dirty_tracked"] is False

    (repo / "parser.py").write_text("v2\n")
    dirty = wrm.code_repo_state(str(repo))
    assert dirty["dirty_tracked"] is True and dirty["dirty_files"] == ["parser.py"]


def test_code_repo_ignores_untracked_and_says_so(tmp_path):
    repo = tmp_path / "shared"
    repo.mkdir()
    git(repo, "init", "-q")
    (repo / "a.py").write_text("x\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "init")
    (repo / "untracked_env_file").write_text("y\n")
    s = wrm.code_repo_state(str(repo))
    assert s["dirty_tracked"] is False
    assert s["untracked_not_checked"] is True


def test_code_repo_not_a_git_repo_is_none_not_guess(tmp_path):
    s = wrm.code_repo_state(str(tmp_path))
    assert s["commit"] is None and s["dirty_tracked"] is None


def test_env_fingerprint_changes_with_any_package_change(tmp_path):
    py = make_env(tmp_path / "e1", ["spliceai-1.3.1-pyh_0", "numpy-1.26.4-py312_0"])
    a = wrm.env_identity(str(py))
    assert a["n_packages"] == 2
    assert a["key_packages"]["spliceai"] == "spliceai-1.3.1-pyh_0"
    assert a["key_packages"]["numpy"] == "numpy-1.26.4-py312_0"

    (tmp_path / "e1" / "conda-meta" / "numpy-1.26.4-py312_0.json").unlink()
    (tmp_path / "e1" / "conda-meta" / "numpy-2.0.0-py312_0.json").write_text("{}")
    b = wrm.env_identity(str(py))
    assert b["n_packages"] == 2
    assert b["fingerprint_sha256"] != a["fingerprint_sha256"]


def test_env_fingerprint_is_order_independent_and_stable(tmp_path):
    p1 = make_env(tmp_path / "e1", ["a-1-0", "b-2-0"])
    p2 = make_env(tmp_path / "e2", ["b-2-0", "a-1-0"])
    assert (wrm.env_identity(str(p1))["fingerprint_sha256"]
            == wrm.env_identity(str(p2))["fingerprint_sha256"])


def test_env_pip_packages_are_fingerprinted_and_reported(tmp_path):
    py = make_env(tmp_path / "e1", ["python-3.12.1-0_cpython"])
    before = wrm.env_identity(str(py))
    sp = tmp_path / "e1" / "lib" / "python3.12" / "site-packages"
    (sp / "torch-2.3.1.dist-info").mkdir(parents=True)
    (sp / "pangolin-1.0.2.dist-info").mkdir()
    (sp / "not_key-9.dist-info").mkdir()
    after = wrm.env_identity(str(py))
    assert after["n_pip_packages"] == 3
    assert after["key_packages"]["torch"] == "pip:torch-2.3.1"
    assert after["key_packages"]["pangolin"] == "pip:pangolin-1.0.2"
    assert "not-key" not in after["key_packages"]
    assert after["fingerprint_sha256"] != before["fingerprint_sha256"]


def test_env_unreadable_is_none_not_guess(tmp_path):
    out = wrm.env_identity(str(tmp_path / "nowhere" / "bin" / "python3"))
    assert out["fingerprint_sha256"] is None and out["n_packages"] is None


def test_launch_writes_external_provenance(tmp_path, monkeypatch):
    repo = tmp_path / "shared"
    repo.mkdir()
    git(repo, "init", "-q")
    (repo / "a.py").write_text("x\n")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "init")
    db = tmp_path / "pangolin.db"
    db.write_bytes(b"db-bytes")
    py = make_env(tmp_path / "env", ["pangolin-1.0.2-py_0"])
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    out = tmp_path / "out"

    monkeypatch.setattr(wrm, "REPO_ROOT", tmp_path / "no_git_here")
    (tmp_path / "no_git_here").mkdir()
    args = argparse.Namespace(
        outdir=str(out), command_line="nextflow run x", nextflow_binary="/bin/true",
        nextflow_main="main.nf", config="nextflow.config", input=str(vcf), run_id="r1",
        profile="local_dev", spip_tier=None, output_format=None, resume=False, bcftools=None,
        vep_sif=None, vep_cache_version=None, spliceai_distance=None, pangolin_distance=None,
        resource=[("pangolin_db", str(db))], code_repo=[("shared_utils", str(repo))],
        env=[("pangolin_python", str(py))],
    )
    wrm.cmd_launch(args)
    ext = json.loads((out / "RUN_MANIFEST.json").read_text())["launch"]["external_provenance"]
    assert ext["resources"]["pangolin_db"]["sha256"] == hashlib.sha256(b"db-bytes").hexdigest()
    assert ext["code_repos"]["shared_utils"]["commit"] == git(repo, "rev-parse", "HEAD").strip()
    assert ext["envs"]["pangolin_python"]["key_packages"]["pangolin"] == "pangolin-1.0.2-py_0"


def test_launch_without_external_args_still_works(tmp_path, monkeypatch):
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    (tmp_path / "no_git_here").mkdir()
    monkeypatch.setattr(wrm, "REPO_ROOT", tmp_path / "no_git_here")
    args = argparse.Namespace(
        outdir=str(tmp_path / "out"), command_line="x", nextflow_binary="/bin/true",
        nextflow_main="main.nf", config="nextflow.config", input=str(vcf), run_id="r1",
        profile="local_dev", spip_tier=None, output_format=None, resume=False, bcftools=None,
        vep_sif=None, vep_cache_version=None, spliceai_distance=None, pangolin_distance=None,
    )
    wrm.cmd_launch(args)
    ext = json.loads((tmp_path / "out" / "RUN_MANIFEST.json").read_text())["launch"]["external_provenance"]
    assert ext == {"resources": {}, "code_repos": {}, "envs": {}}
