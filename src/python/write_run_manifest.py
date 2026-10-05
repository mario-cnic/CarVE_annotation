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

`--outdir` is the run folder. `launch` creates <outdir>/pipeline_info/RUN_MANIFEST.json (runs made
before the sub-folder layout have it directly in <outdir>); `completed` updates the same file in
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
INFO_SUBDIR, REPORTS_SUBDIR, GENOTYPES_SUBDIR = "pipeline_info", "reports", "genotypes"


def run_files(outdir: Path, subdir: str, pattern: str) -> list[Path]:
    """Files matching `pattern` in <outdir>/<subdir>, else in <outdir> itself (runs made before the sub-folder layout)."""
    return sorted((outdir / subdir).glob(pattern)) or sorted(outdir.glob(pattern))


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


# External inputs (databases, models, reference FASTA) are hashed in full up to this size; above it
# only stat() facts plus any sidecar index/checksum are recorded, so a launch never re-reads a
# multi-GB mounted file (the 3 GB FASTA). The Pangolin db (~0.4 GB) is under the limit.
RESOURCE_HASH_MAX_BYTES = 1 << 30
# Packages read from an env's conda-meta (no execution) and reported by name; the fingerprint
# covers the whole env, this list is only for human-readable versions.
KEY_ENV_PACKAGES = ("spliceai", "pangolin", "torch", "pytorch", "tensorflow", "keras", "pysam",
                    "bcftools", "htslib", "samtools", "pandas", "pyarrow", "numpy", "r-base",
                    "gffutils", "pyfaidx", "python")


def split_named(spec: str) -> tuple[str, str]:
    name, sep, value = spec.partition("=")
    if not sep or not name or not value:
        raise argparse.ArgumentTypeError(f"expected NAME=PATH, got {spec!r}")
    return name, value


def resource_identity(path_str: str) -> dict:
    """Identity of an external file/dir. Never guesses: a field that cannot be determined is None."""
    p = Path(path_str)
    out = {"path": str(p), "exists": p.exists(), "is_dir": p.is_dir(), "size_bytes": None,
           "mtime_utc": None, "sha256": None, "hash_skipped": None, "sidecars": {}}
    if not p.exists():
        return out
    st = p.stat()
    out["mtime_utc"] = datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat()
    if p.is_dir():
        return out
    out["size_bytes"] = st.st_size
    if st.st_size <= RESOURCE_HASH_MAX_BYTES:
        out["sha256"] = sha256_of(p)
    else:
        out["hash_skipped"] = f"size > {RESOURCE_HASH_MAX_BYTES} bytes; sidecars recorded instead"
    for suffix in (".fai", ".md5", ".sha256"):
        side = Path(str(p) + suffix)
        if side.is_file():
            out["sidecars"][suffix] = {
                "sha256": sha256_of(side),
                "content": side.read_text(errors="replace").strip()[:200] if suffix != ".fai" else None,
            }
    return out


def code_repo_state(path_str: str) -> dict:
    """Commit + tracked-file modifications of a code repo outside this one (e.g. shared/utils).
    `--untracked-files=no` on purpose: that repo also contains the conda envs, and scanning them
    for untracked files would make every launch slow. Untracked code there is therefore NOT seen."""
    root = Path(path_str)
    out = {"path": str(root), "commit": None, "dirty_tracked": None, "dirty_files": [],
           "untracked_not_checked": True}
    top = sh(["git", "rev-parse", "--show-toplevel"], cwd=root) if root.is_dir() else None
    if not top:
        return out
    out["git_toplevel"] = top
    out["commit"] = sh(["git", "rev-parse", "HEAD"], cwd=root)
    try:
        raw = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=root,
                             capture_output=True, text=True, check=True, timeout=60).stdout.rstrip("\n")
    except Exception:
        return out
    files = [l[3:] for l in raw.split("\n") if l]
    out["dirty_tracked"] = bool(files)
    out["dirty_files"] = files[:50]
    out["dirty_files_total"] = len(files)
    return out


def env_identity(bin_path: str) -> dict:
    """Fingerprint of the conda env behind a binary, read from conda-meta without running anything.
    The fingerprint is the sha256 of the sorted package-build file names, so ANY package change
    (add/remove/upgrade/rebuild) alters it. None when the env is not readable from here."""
    out = {"binary": bin_path, "env_dir": None, "n_packages": None, "fingerprint_sha256": None,
           "key_packages": {}}
    b = Path(bin_path)
    env_dir = b.parent.parent
    meta = env_dir / "conda-meta"
    if b.parent.name != "bin" or not meta.is_dir():
        return out
    names = sorted(f.name[:-5] for f in meta.glob("*.json"))
    # pip-installed packages (e.g. torch, pangolin, spliceai) are invisible to conda-meta; their
    # *.dist-info directories are the record. Prefixed so they cannot collide with conda names.
    pip_names = sorted("pip:" + d.name[:-len(".dist-info")]
                       for d in env_dir.glob("lib/python*/site-packages/*.dist-info"))
    out["env_dir"] = str(env_dir)
    out["n_packages"] = len(names)
    out["n_pip_packages"] = len(pip_names)
    out["fingerprint_sha256"] = hashlib.sha256("\n".join(names + pip_names).encode()).hexdigest()
    for n in names:
        pkg = n.rsplit("-", 2)[0] if n.count("-") >= 2 else n
        if pkg in KEY_ENV_PACKAGES:
            out["key_packages"][pkg] = n
    for n in pip_names:
        pkg = n[len("pip:"):].split("-", 1)[0].lower().replace("_", "-")
        if pkg in KEY_ENV_PACKAGES:
            out["key_packages"].setdefault(pkg, n)
    return out


def cmd_launch(args):
    outdir = Path(args.outdir)
    info_dir = outdir / INFO_SUBDIR
    info_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = info_dir / "RUN_MANIFEST.json"
    previous = [p for p in (manifest_path, outdir / "RUN_MANIFEST.json") if p.exists()]
    if previous:
        print(
            f"warning: {previous[0]} already exists — this outdir was used before "
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
                REPO_ROOT, info_dir, launched_at.strftime("%Y%m%dT%H%M%SZ"))
            print(f"warning: dirty working tree; launch-time changes saved to "
                  f"{info_dir / repo_state['dirty_patch']['path']}", file=sys.stderr)
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
            # Everything outside this repo that changes the output without changing this repo's
            # commit (TODO.md "Nextflow pipeline: make the run reproducible from the repo commit").
            "external_provenance": {
                "resources": {n: resource_identity(v) for n, v in (getattr(args, "resource", None) or [])},
                "code_repos": {n: code_repo_state(v) for n, v in (getattr(args, "code_repo", None) or [])},
                "envs": {n: env_identity(v) for n, v in (getattr(args, "env", None) or [])},
            },
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


def read_input_check(outdir: Path) -> dict | None:
    """Genotype mode, samples and normalisation counts from the published *.assembly_check.tsv."""
    reports = run_files(outdir, REPORTS_SUBDIR, "*.assembly_check.tsv")
    if not reports:
        return None
    keep = ("status", "genotype_mode", "n_samples", "sample_names", "format_fields_declared", "records_scanned",
            "multiallelic_records", "not_normalised_records", "duplicate_site_records")
    values = {}
    for line in reports[0].read_text().splitlines()[1:]:
        key, _, value = line.partition("\t")
        if key in keep:
            values[key] = value
    values["report"] = reports[0].name
    return values


def cmd_completed(args):
    outdir = Path(args.outdir)
    manifest_path = next((p for p in (outdir / INFO_SUBDIR / "RUN_MANIFEST.json", outdir / "RUN_MANIFEST.json")
                          if p.exists()), outdir / INFO_SUBDIR / "RUN_MANIFEST.json")
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
        "input_check": read_input_check(outdir),
        "genotype_outputs": {p.name: resource_identity(str(p)) for p in run_files(outdir, GENOTYPES_SUBDIR, "*.genotypes.pq")},
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
    p_launch.add_argument("--resource", action="append", type=split_named, metavar="NAME=PATH",
                          help="external file/dir to identify (size, mtime, sha256 if small, sidecars)")
    p_launch.add_argument("--code-repo", action="append", type=split_named, metavar="NAME=PATH",
                          help="git repo outside this one whose code the run executes")
    p_launch.add_argument("--env", action="append", type=split_named, metavar="NAME=BINARY",
                          help="binary inside a conda env; the env is fingerprinted from conda-meta")
    p_launch.set_defaults(func=cmd_launch)

    p_completed = sub.add_parser("completed")
    p_completed.add_argument("--outdir", required=True)
    p_completed.add_argument("--exit-code", required=True, type=int)
    p_completed.add_argument("--nextflow-log")
    p_completed.set_defaults(func=cmd_completed)

    ns = ap.parse_args()
    ns.func(ns)
