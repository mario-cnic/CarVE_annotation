import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import write_run_manifest as wrm


def git(repo, *args):
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", *args],
                          cwd=repo, capture_output=True, check=True, text=True).stdout


@pytest.fixture
def repo(tmp_path):
    r = tmp_path / "repo"
    (r / "src").mkdir(parents=True)
    (r / "resources").mkdir()
    git(r, "init", "-q")
    (r / ".gitignore").write_text("nf_work/\n")
    (r / "main.nf").write_text("workflow { }\n")
    (r / "src" / "tool.py").write_text("print('v1')\n")
    (r / "src" / "old.py").write_text("print('to delete')\n")
    (r / "resources" / "blob.bin").write_bytes(bytes(range(256)))
    git(r, "add", "-A")
    git(r, "commit", "-q", "-m", "init")
    return r


def launch_args(outdir, input_path):
    return argparse.Namespace(
        outdir=str(outdir), command_line="nextflow run main.nf", nextflow_binary="/bin/true",
        nextflow_main="main.nf", config="nextflow.config", input=str(input_path), run_id="r1",
        profile="local_dev", spip_tier=None, output_format=None, resume=False, bcftools=None,
        vep_sif=None, vep_cache_version=None, spliceai_distance=None, pangolin_distance=None,
    )


def make_dirty(r):
    (r / "src" / "tool.py").write_text("print('v2, uncommitted')\n")
    (r / "src" / "old.py").unlink()
    (r / "resources" / "blob.bin").write_bytes(bytes(reversed(range(256))))
    (r / "src" / "new_module.py").write_text("NEW = 1\n")
    (r / "resources" / "huge.bed").write_text("x" * 2048)
    (r / "notes.txt").write_text("scratch\n")
    (r / "nf_work").mkdir()
    (r / "nf_work" / "ignored.txt").write_text("run output\n")


def test_clean_launch_writes_no_patch(repo, tmp_path, monkeypatch):
    monkeypatch.setattr(wrm, "REPO_ROOT", repo)
    out = repo / "nf_work" / "out"
    wrm.cmd_launch(launch_args(out, repo / "main.nf"))
    m = json.loads((out / "RUN_MANIFEST.json").read_text())["launch"]["annotation_pipeline_repo"]
    assert m["dirty"] is False and m["dirty_patch"] is None
    assert not list(out.glob("RUN_DIRTY.*.patch"))


def test_dirty_launch_patch_reconstructs_tree(repo, tmp_path, monkeypatch):
    monkeypatch.setattr(wrm, "REPO_ROOT", repo)
    monkeypatch.setattr(wrm, "UNTRACKED_MAX_BYTES", 1024)
    make_dirty(repo)
    status_before = git(repo, "status", "--porcelain")
    head = git(repo, "rev-parse", "HEAD").strip()

    out = repo / "nf_work" / "out"
    wrm.cmd_launch(launch_args(out, repo / "main.nf"))

    # The real index and working tree are untouched.
    assert git(repo, "status", "--porcelain") == status_before
    assert git(repo, "diff", "--cached", "--name-only") == ""

    m = json.loads((out / "RUN_MANIFEST.json").read_text())["launch"]["annotation_pipeline_repo"]
    assert m["dirty"] is True and m["commit"] == head
    dp = m["dirty_patch"]
    patch = out / dp["path"]
    assert patch.is_file() and dp["sha256"] == wrm.sha256_of(patch)
    assert dp["untracked_included"] == ["src/new_module.py"]
    not_captured = {d["path"]: d for d in dp["untracked_not_captured"]}
    assert not_captured["resources/huge.bed"]["sha256"] == wrm.sha256_of(repo / "resources" / "huge.bed")
    assert not_captured["notes.txt"]["reason"] == "outside capture pathspecs"
    assert "nf_work/ignored.txt" not in not_captured

    # Rebuild from the recorded commit + patch in a fresh clone.
    clone = tmp_path / "clone"
    subprocess.run(["git", "clone", "-q", str(repo), str(clone)], check=True)
    git(clone, "checkout", "-q", head)
    git(clone, "apply", str(patch))
    assert (clone / "src" / "tool.py").read_text() == "print('v2, uncommitted')\n"
    assert not (clone / "src" / "old.py").exists()
    assert (clone / "resources" / "blob.bin").read_bytes() == bytes(reversed(range(256)))
    assert (clone / "src" / "new_module.py").read_text() == "NEW = 1\n"
    assert not (clone / "resources" / "huge.bed").exists()


def test_patch_names_are_unique_per_launch(repo, monkeypatch):
    monkeypatch.setattr(wrm, "REPO_ROOT", repo)
    make_dirty(repo)
    out = repo / "nf_work" / "out"
    first = wrm.capture_dirty_tree(repo, out.parent, "20260101T000000Z")["path"]
    second = wrm.capture_dirty_tree(repo, out.parent, "20260101T000001Z")["path"]
    assert first != second
