"""
Run-provenance manifest for annotation_pipeline_new (Module 2) — closes TODO.md's
Priority-0 item and matches the schema Module 1 (sarek_pipeline) already committed
to for its own PLT-012 (`sarek_pipeline/src/write_run_manifest.py`). Each module
writes its own copy per DEC-0002 (repos stay separate) — this is not a shared
library, just a same-schema sibling.

Usage (called by run_annotate_vcf.sh, not normally invoked directly):
    write_run_manifest.py launch --outdir DIR --command-line "..." \
        --nextflow-main PATH --nextflow-binary PATH --config PATH --input PATH \
        --run-id ID [--profile local_dev] [--spip-tier restricted] \
        [--output-format pq] [--resume]
    write_run_manifest.py completed --outdir DIR --exit-code N \
        --nextflow-log PATH

`launch` creates <outdir>/RUN_MANIFEST.json; `completed` updates the same file in
place (fails loudly if the launch record is missing, rather than silently
creating a partial one).
"""
from __future__ import annotations  # cluster's system python3 predates 3.10's `X | None`
                                     # syntax (confirmed by hitting exactly that TypeError on a
                                     # real cluster run) — this defers annotation evaluation
                                     # entirely instead of rewriting every hint, safe on 3.7+.

import argparse
import hashlib
import json
import re
import shutil
import os
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

SCHEMA_VERSION = 1
REPO_ROOT = Path(__file__).resolve().parent.parent.parent  # annotation_pipeline_new/


def sh(cmd: list[str], cwd: Path | None = None) -> str | None:
    try:
        out = subprocess.run(
            cmd, cwd=cwd, capture_output=True, text=True, check=True, timeout=30
        )
        return out.stdout.strip()
    except Exception:
        return None


def sha256_of(path: Path) -> str | None:
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_state(repo_root: Path) -> dict:
    commit = sh(["git", "rev-parse", "HEAD"], cwd=repo_root)
    # NOTE: don't route this through sh()'s blanket .strip() — porcelain's first
    # two columns are meaningful status characters (a leading space is not
    # incidental whitespace), same rationale as sarek_pipeline's copy of this.
    try:
        raw = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=repo_root, capture_output=True, text=True, check=True, timeout=30,
        ).stdout.rstrip("\n")
    except Exception:
        raw = ""
    dirty_files = [l[3:] for l in raw.split("\n") if l] if raw else []
    return {
        "commit": commit,
        "dirty": bool(dirty_files),
        "dirty_files": dirty_files,
    }


# Untracked files under these pathspecs are code/config a run can execute or read, so their
# content goes into the dirty-tree patch; anything else untracked is listed by name only.
UNTRACKED_CAPTURE_PATHSPECS = ["*.nf", "*.config", "*.sh", "modules", "workflows", "src", "config", "resources"]
UNTRACKED_MAX_BYTES = 5 * 1024 * 1024


def _git(args: list[str], repo_root: Path, env: dict | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", *args], cwd=repo_root, capture_output=True, check=True, timeout=300,
        env={**os.environ, **(env or {})},
    )


def capture_dirty_tree(repo_root: Path, outdir: Path, stamp: str) -> dict:
    """Writes <outdir>/RUN_DIRTY.<stamp>.patch, a binary diff from HEAD to the launch-time tree.

    Covers every tracked change (including deletions) plus untracked, non-ignored files under
    UNTRACKED_CAPTURE_PATHSPECS up to UNTRACKED_MAX_BYTES each. Reconstruct with
    `git checkout <commit> && git apply <patch>`. Staging happens in a temporary GIT_INDEX_FILE,
    so the repository's real index is untouched (blobs are written to the object store).
    """
    split = lambda b: [p for p in b.decode().split("\0") if p]
    in_scope = split(_git(["ls-files", "--others", "--exclude-standard", "-z", "--",
                           *UNTRACKED_CAPTURE_PATHSPECS], repo_root).stdout)
    all_untracked = split(_git(["ls-files", "--others", "--exclude-standard", "-z"], repo_root).stdout)

    included, skipped = [], []
    for rel in in_scope:
        size = (repo_root / rel).stat().st_size
        if size <= UNTRACKED_MAX_BYTES:
            included.append(rel)
        else:
            skipped.append({"path": rel, "reason": f"larger than {UNTRACKED_MAX_BYTES} bytes",
                            "sha256": sha256_of(repo_root / rel)})
    skipped += [{"path": rel, "reason": "outside capture pathspecs", "sha256": None}
                for rel in all_untracked if rel not in set(in_scope)]

    patch_path = outdir / f"RUN_DIRTY.{stamp}.patch"
    with tempfile.TemporaryDirectory() as tmp:
        env = {"GIT_INDEX_FILE": str(Path(tmp) / "index")}
        _git(["read-tree", "HEAD"], repo_root, env)
        _git(["add", "-u", "--", "."], repo_root, env)
        if included:
            _git(["add", "--", *included], repo_root, env)
        patch = _git(["diff", "--cached", "--binary", "HEAD"], repo_root, env).stdout
    patch_path.write_bytes(patch)
    return {
        "path": patch_path.name,
        "sha256": sha256_of(patch_path),
        "bytes": len(patch),
        "untracked_included": included,
        "untracked_not_captured": skipped,
        "apply_with": "git checkout <commit> && git apply <path>",
    }


def resolve(path_str: str | None) -> str | None:
    if not path_str:
        return None
    return str(Path(path_str).resolve())


def tool_version(binary: str | None, version_flag: str = "--version") -> str | None:
    """Best-effort: query a tool's own version by actually running it. Returns
    None (not a guess) if the binary isn't executable from here — e.g. every
    predictor's own conda env, which this pipeline's `standard` profile only
    ever reaches on the real cluster, not from a local/dev sandbox."""
    if not binary:
        return None
    out = sh([binary, version_flag])
    if out is None:
        return None
    return out.splitlines()[0] if out else None


def cmd_launch(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "RUN_MANIFEST.json"
    if manifest_path.exists():
        print(
            f"warning: {manifest_path} already exists — this outdir was used before "
            "(expected on a -resume of the same run_id; if this is meant to be a fresh "
            "run, use a different --run_id instead of overwriting history)",
            file=sys.stderr,
        )

    # Resolved by run_annotate_vcf.sh, not re-derived here via shutil.which("nextflow") — this
    # cluster's nextflow lives at /opt/nextflow/nextflow, not on PATH, so a bare `which` lookup
    # from this script would silently record None even on a successful real run. Passing it
    # through explicitly means the manifest reflects the binary that ACTUALLY ran, not a guess.
    nextflow_bin = args.nextflow_binary or shutil.which("nextflow")

    launched_at = datetime.now(timezone.utc)
    repo_state = git_state(REPO_ROOT)
    repo_state["dirty_patch"] = None
    if repo_state["dirty"]:
        try:
            repo_state["dirty_patch"] = capture_dirty_tree(
                REPO_ROOT, outdir, launched_at.strftime("%Y%m%dT%H%M%SZ"))
            print(f"warning: dirty working tree; launch-time changes saved to "
                  f"{outdir / repo_state['dirty_patch']['path']}", file=sys.stderr)
        except Exception as e:
            sys.exit(f"error: working tree is dirty and the dirty-tree patch could not be written ({e}); "
                     "commit or stash first so the run's code can be reconstructed")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "module": "2 (annotation_pipeline_new)",
        "launch": {
            "generated_at": launched_at.isoformat(),
            "command_line": args.command_line,
            "cli_flags": {
                "run_id": args.run_id,
                "profile": args.profile,
                "spip_tier": args.spip_tier,
                "output_format": args.output_format,
                "resume": bool(args.resume),
            },
            "resolved_paths": {
                "nextflow_binary": nextflow_bin,
                "nextflow_main": resolve(args.nextflow_main),
                "config": resolve(args.config),
                "input": resolve(args.input),
                "outdir": resolve(args.outdir),
            },
            "file_hashes_sha256": {
                "input": sha256_of(Path(args.input)) if args.input else None,
                "gene_transcript_mapping": sha256_of(REPO_ROOT / "resources" / "gene_transcript_mapping.txt"),
                "v5_genes_loc_bed": sha256_of(REPO_ROOT / "resources" / "v5_genes_loc.bed"),
            },
            "annotation_pipeline_repo": repo_state,
            "nextflow_version": sh([nextflow_bin, "-v"]) if nextflow_bin else None,
            # Best-effort tool versions, queried at run time rather than hardcoded — closes
            # TODO.md's Priority-0 complaint that no file anywhere records what actually ran.
            # None here means "not queryable from where this manifest was written" (this
            # sandbox has no cluster access to most of these binaries), not "unversioned". VEP
            # itself isn't queried as a local binary — it's containerized (params.vep_sif), so
            # its container's own file hash is the meaningful provenance artifact, not a `vep
            # --help` call against a binary that doesn't exist outside the container.
            "tool_versions": {
                "bcftools": tool_version(args.bcftools, "--version"),
                "vep_cache_version": args.vep_cache_version,
                "vep_sif_sha256": sha256_of(Path(args.vep_sif)) if args.vep_sif else None,
            },
            # SpliceAI/Pangolin's -d/--distance value actually used — closes the specific
            # ambiguity TODO.md named ("the SpliceAI -D question ... there was nothing to
            # check"). Read directly from the resolved config, not retyped here.
            "spliceai_distance": args.spliceai_distance,
            "pangolin_distance": args.pangolin_distance,
            "launched_by": sh(["whoami"]),
            "host": sh(["hostname"]),
            "launch_cwd": str(Path.cwd()),
        },
        "completed": None,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {manifest_path}")


def cmd_completed(args):
    outdir = Path(args.outdir)
    manifest_path = outdir / "RUN_MANIFEST.json"
    if not manifest_path.exists():
        sys.exit(
            f"error: {manifest_path} not found — 'launch' stage must run first "
            "(refusing to write a completion record with no launch record to attach it to)"
        )
    manifest = json.loads(manifest_path.read_text())

    workflow_stats = None
    if args.nextflow_log and Path(args.nextflow_log).is_file():
        text = Path(args.nextflow_log).read_text(errors="replace")
        m = re.search(r"WorkflowStats\[[^\]]*\]", text)
        if m:
            workflow_stats = m.group(0)

    manifest["completed"] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "exit_code": args.exit_code,
        "workflow_stats": workflow_stats,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"updated {manifest_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="stage", required=True)

    p_launch = sub.add_parser("launch")
    p_launch.add_argument("--outdir", required=True)
    p_launch.add_argument("--command-line", required=True)
    p_launch.add_argument("--nextflow-main", required=True)
    p_launch.add_argument("--nextflow-binary")
    p_launch.add_argument("--config", required=True)
    p_launch.add_argument("--input", required=True)
    p_launch.add_argument("--run-id")
    p_launch.add_argument("--profile")
    p_launch.add_argument("--spip-tier")
    p_launch.add_argument("--output-format")
    p_launch.add_argument("--bcftools")
    p_launch.add_argument("--vep-sif")
    p_launch.add_argument("--vep-cache-version")
    p_launch.add_argument("--spliceai-distance")
    p_launch.add_argument("--pangolin-distance")
    p_launch.add_argument("--resume", action="store_true")
    p_launch.set_defaults(func=cmd_launch)

    p_completed = sub.add_parser("completed")
    p_completed.add_argument("--outdir", required=True)
    p_completed.add_argument("--exit-code", required=True, type=int)
    p_completed.add_argument("--nextflow-log")
    p_completed.set_defaults(func=cmd_completed)

    ns = ap.parse_args()
    ns.func(ns)
