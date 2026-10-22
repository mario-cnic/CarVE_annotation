"""
Predictor score cache for SpliceAI and Pangolin. Each (predictor, version) has an append-only
directory of sites-only VCF segments:

    <cache_dir>/<predictor>/<version_id>/VERSION.json
    <cache_dir>/<predictor>/<version_id>/segments/<name>.vcf.gz(.tbi) + <name>.json

A segment holds one record per site the predictor evaluated (biallelic, normalised) with only
that predictor's INFO fields. A record without the field means "computed, no score"; a site with
no record is not in the cache. version_id is derived from the run manifest (FASTA, predictor
database, conda env fingerprint, distance) and the scoring code at the run's commit.

Usage:
    annotation_cache.py version-id --manifest RUN_MANIFEST.json --predictor pangolin --cache-dir DIR
    annotation_cache.py split --predictor pangolin --input subset.vcf.gz --segments-dir DIR \
        --hits hits.vcf.gz --misses misses.vcf.gz --report report.json
    annotation_cache.py make-segment --predictor pangolin --annotated out.vcf.gz \
        --misses misses.vcf.gz --output-prefix segments/<name> --run-id RUN --version-id ID
    annotation_cache.py seed --run-dir RUN_DIR --predictor pangolin --cache-dir DIR \
        [--into-version ID --override-reason TEXT]
    annotation_cache.py verify --version-dir DIR --predictor pangolin
    annotation_cache.py list-segments --version-dir DIR --output segments.txt
"""
from __future__ import annotations  # cluster's system python3 runs the launcher-side commands

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from itertools import groupby
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Per predictor: cached INFO fields and their header lines (as the predictor writes them), manifest
# keys of its identity, scoring code paths, and the Nextflow process whose `script:` block (the
# command line) is part of the identity.
PREDICTORS = {
    "spliceai": {
        "fields": ["SpliceAI"],
        "info_header": [
            '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI predicted splice scores '
            '(SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL)">'],
        "resources": [],
        "env": "python_spliceai",
        "distance": "spliceai_distance",
        "code": ["src/python/annotate_spliceai.py"],
        "module": "modules/local/spliceai.nf",
        "process": "SPLICEAI_ANNOTATE",
        "output_suffix": "annSpliceAI",
    },
    "pangolin": {
        "fields": ["Pangolin", "PangolinTissue"],
        "info_header": [
            '##INFO=<ID=Pangolin,Number=.,Type=String,Description="Pangolin splice score predictions">',
            '##INFO=<ID=PangolinTissue,Number=.,Type=String,Description="Which of Pangolin\'s 4 trained GTEx '
            'tissues (Heart, Liver, Brain, Testis) produced the reported gain/loss score in the Pangolin field, '
            'per gene: gene|gain_tissue|loss_tissue">'],
        "resources": ["pangolin_db"],
        "env": "pangolin_python",
        "distance": "pangolin_distance",
        "code": ["src/external/Pangolin-main/pangolin"],
        "module": "modules/local/pangolin.nf",
        "process": "PANGOLIN_ANNOTATE",
        "output_suffix": "annPangolin",
    },
}


class CacheError(Exception):
    pass


def tool(name: str, override: str | None) -> str:
    path = override or shutil.which(name)
    if not path:
        raise CacheError(f"{name} not found; pass --{name}")
    return path


def run(cmd: list[str], **kw) -> str:
    res = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if res.returncode != 0:
        raise CacheError(f"command failed ({res.returncode}): {' '.join(cmd)}\n{res.stderr.strip()}")
    return res.stdout


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json_atomic(path: Path, data: dict) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ---------------------------------------------------------------- version identity

def script_block(nf_text: str, process: str) -> str:
    """The `script:` section of `process`, up to its `stub:` or the end of the process."""
    body = re.search(r"^process\s+" + re.escape(process) + r"\s*\{(.*?)^}\s*$", nf_text, re.S | re.M)
    m = body and re.search(r"^\s*script:\s*$(.*?)(?=^\s*stub:\s*$|\Z)", body.group(1), re.S | re.M)
    if not m:
        raise CacheError(f"no script: block found for process {process}")
    return m.group(1)


def code_identity(repo_root: Path, commit: str, predictor: str) -> dict:
    """Git blob ids of the scoring code at `commit` (hidden files and caches skipped) plus the
    sha256 of the predictor module's script block."""
    spec = PREDICTORS[predictor]
    out = {}
    for path in spec["code"]:
        listing = run(["git", "ls-tree", "-r", "--full-tree", commit, "--", path], cwd=repo_root)
        n = 0
        for line in listing.splitlines():
            meta, fpath = line.split("\t", 1)
            parts = Path(fpath).parts
            if any(p.startswith(".") or p == "__pycache__" for p in parts) or fpath.endswith(".pyc"):
                continue
            out[fpath] = meta.split()[2]
            n += 1
        if n == 0:
            raise CacheError(f"{path} not tracked at commit {commit}")
    nf = run(["git", "show", f"{commit}:{spec['module']}"], cwd=repo_root)
    out[spec["module"] + "#script"] = hashlib.sha256(script_block(nf, spec["process"]).encode()).hexdigest()
    return out


def version_identity(manifest: dict, predictor: str, repo_root: Path | None = None) -> dict:
    """Identity fields of `predictor` for the run described by `manifest`. Raises when any field
    is unknown: a score with unknown provenance is not cached."""
    spec = PREDICTORS[predictor]
    launch = manifest.get("launch") or {}
    repo = launch.get("annotation_pipeline_repo") or {}
    commit = repo.get("commit")
    if not commit:
        raise CacheError("manifest has no annotation_pipeline_repo.commit")
    code_paths = spec["code"] + [spec["module"]]
    dirty = [f for f in repo.get("dirty_files") or []
             if any(f.split(maxsplit=1)[-1].startswith(p) for p in code_paths)]
    if dirty:
        raise CacheError(f"uncommitted changes in {predictor} scoring code: {dirty}; commit them first")

    ext = launch.get("external_provenance") or {}
    resources = ext.get("resources") or {}
    fasta = resources.get("fasta") or {}
    fai = ((fasta.get("sidecars") or {}).get(".fai") or {}).get("sha256")
    if not fasta.get("size_bytes") or not (fasta.get("sha256") or fai):
        raise CacheError("manifest has no FASTA identity (size and sha256 or .fai sha256)")
    ident = {
        "predictor": predictor,
        "fasta": {"size_bytes": fasta["size_bytes"], "sha256": fasta.get("sha256"), "fai_sha256": fai},
        "resources": {},
        "env_fingerprint_sha256": ((ext.get("envs") or {}).get(spec["env"]) or {}).get("fingerprint_sha256"),
        "distance": launch.get(spec["distance"]),
        "code": code_identity(repo_root or REPO_ROOT, commit, predictor),
    }
    for r in spec["resources"]:
        sha = (resources.get(r) or {}).get("sha256")
        if not sha:
            raise CacheError(f"manifest has no sha256 for resource {r}")
        ident["resources"][r] = sha
    if not ident["env_fingerprint_sha256"]:
        raise CacheError(f"manifest has no env fingerprint for {spec['env']}")
    if ident["distance"] in (None, ""):
        raise CacheError(f"manifest has no {spec['distance']}")
    ident["distance"] = str(ident["distance"])
    return ident


def version_id(identity: dict) -> str:
    canon = json.dumps(identity, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canon.encode()).hexdigest()[:16]


def ensure_version_dir(cache_dir: Path, identity: dict) -> Path:
    """Creates <cache_dir>/<predictor>/<version_id>/ with VERSION.json, or checks the existing one."""
    vid = version_id(identity)
    vdir = cache_dir / identity["predictor"] / vid
    vfile = vdir / "VERSION.json"
    if vfile.exists():
        if json.loads(vfile.read_text()).get("identity") != identity:
            raise CacheError(f"{vfile} identity differs from this run's")
        return vdir
    (vdir / "segments").mkdir(parents=True, exist_ok=True)
    write_json_atomic(vfile, {"version_id": vid, "created_utc": now_utc(), "identity": identity})
    return vdir


def load_manifest(path: Path) -> dict:
    if not path.is_file():
        raise CacheError(f"manifest not found: {path}")
    return json.loads(path.read_text())


# ---------------------------------------------------------------- segments

def list_segments(segments_dir: Path) -> tuple[list[Path], list[str]]:
    """Segments whose .vcf.gz matches its sidecar (size, sha256) and has an index. A segment still
    being copied or otherwise damaged is skipped, so its sites are recomputed, never misread."""
    good, skipped = [], []
    for side in sorted(segments_dir.glob("*.json")):
        vcf = side.with_suffix(".vcf.gz")
        try:
            meta = json.loads(side.read_text())
        except ValueError:
            skipped.append(f"{side.name}: unreadable sidecar")
            continue
        if not vcf.is_file() or not Path(str(vcf) + ".tbi").is_file():
            skipped.append(f"{vcf.name}: missing file or index")
        elif vcf.stat().st_size != meta.get("size_bytes") or sha256_file(vcf) != meta.get("sha256"):
            skipped.append(f"{vcf.name}: size or sha256 differs from sidecar")
        else:
            good.append(vcf)
    return good, skipped


def write_segment_list(version_dir: Path, output: Path) -> tuple[int, list[str]]:
    """Absolute paths of every file of the valid segments, one per line: the frozen set a run reads."""
    segments, skipped = list_segments(version_dir / "segments")
    paths = []
    for vcf in segments:
        paths += [str(vcf.resolve()), str(vcf.resolve()) + ".tbi", str(vcf.resolve().with_name(vcf.name[:-len(".vcf.gz")] + ".json"))]
    output.write_text("".join(p + "\n" for p in paths))
    return len(segments), skipped


def build_lookup(bcftools: str, tabix: str, segments: list[Path], out: Path) -> Path:
    """One indexed file over all segments; duplicate sites (same CHROM/POS/REF/ALT) kept once."""
    if len(segments) == 1:
        return segments[0]
    run([bcftools, "concat", "-a", "-D", *map(str, segments), "-Oz", "-o", str(out)])
    run([tabix, "-f", "-p", "vcf", str(out)])
    return out


def count_records(bcftools: str, vcf: Path) -> int:
    out = subprocess.run([bcftools, "index", "-n", str(vcf)], capture_output=True, text=True)
    if out.returncode == 0 and out.stdout.strip():
        return int(out.stdout.strip())
    return sum(1 for _ in run([bcftools, "view", "-H", str(vcf)]).splitlines())


def split(predictor: str, input_vcf: Path, segments_dir: Path, hits: Path, misses: Path,
          report: Path, bcftools: str, tabix: str) -> dict:
    """Splits `input_vcf` into hits (records in the cache, predictor INFO copied from it) and
    misses (records not in the cache). Both outputs are always written and indexed."""
    fields = PREDICTORS[predictor]["fields"]
    segments, skipped = list_segments(segments_dir) if segments_dir.is_dir() else ([], [])
    with tempfile.TemporaryDirectory(dir=hits.parent) as tmp:
        tmp = Path(tmp)
        if segments:
            lookup = build_lookup(bcftools, tabix, segments, tmp / "lookup.vcf.gz")
            shared = tmp / "shared.vcf.gz"
            run([bcftools, "isec", "-c", "none", "-n=2", "-w1", str(input_vcf), str(lookup), "-Oz", "-o", str(shared)])
            run([tabix, "-f", "-p", "vcf", str(shared)])
            # Records are biallelic, so annotate's default allele matching is exact.
            run([bcftools, "annotate", "-a", str(lookup),
                 "-c", ",".join("INFO/" + f for f in fields), str(shared), "-Oz", "-o", str(hits)])
            run([bcftools, "isec", "-c", "none", "-C", "-w1", str(input_vcf), str(lookup), "-Oz", "-o", str(misses)])
        else:
            run([bcftools, "view", "-h", str(input_vcf), "-Oz", "-o", str(hits)])
            run([bcftools, "view", str(input_vcf), "-Oz", "-o", str(misses)])
    for f in (hits, misses):
        run([tabix, "-f", "-p", "vcf", str(f)])
    n_in, n_hit, n_miss = (count_records(bcftools, f) for f in (input_vcf, hits, misses))
    if n_hit + n_miss != n_in:
        raise CacheError(f"split parity failed: input={n_in} hits={n_hit} misses={n_miss}")
    result = {"predictor": predictor, "input_records": n_in, "hits": n_hit, "misses": n_miss,
              "segments_used": [s.name for s in segments], "segments_skipped": skipped}
    write_json_atomic(report, result)
    return result


def segment_header(bcftools: str, src: Path, predictor: str, extra: dict) -> str:
    """fileformat, contigs, the predictor's INFO lines (the source's, else the canonical ones) and
    the cache metadata; nothing sample- or run-specific from the source header."""
    spec = PREDICTORS[predictor]
    lines, info = [], {}
    for line in run([bcftools, "view", "-h", str(src)]).splitlines():
        if line.startswith(("##fileformat=", "##contig=")):
            lines.append(line)
        for f in spec["fields"]:
            if line.startswith(f"##INFO=<ID={f},"):
                info[f] = line
    for f, canonical in zip(spec["fields"], spec["info_header"]):
        lines.append(info.get(f, canonical))
    lines += [f"##annotation_cache_{k}={v}" for k, v in extra.items() if v is not None]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return "\n".join(lines) + "\n"


def write_segment(predictor: str, src: Path, out_prefix: Path, meta: dict,
                  bcftools: str, tabix: str, bgzip: str) -> dict:
    """Writes <out_prefix>.vcf.gz(.tbi) with the sites of `src` stripped to the predictor fields,
    then the <out_prefix>.json sidecar last (its presence marks the segment complete)."""
    fields = PREDICTORS[predictor]["fields"]
    vcf = out_prefix.with_name(out_prefix.name + ".vcf.gz")
    header = segment_header(bcftools, src, predictor,
                            {"predictor": predictor, "version_id": meta.get("version_id"),
                             "run_id": meta.get("run_id")})
    strip = subprocess.Popen([bcftools, "view", "-G", str(src), "-Ou"], stdout=subprocess.PIPE)
    body = subprocess.Popen([bcftools, "annotate", "-x", "ID,QUAL,FILTER,^" + ",".join("INFO/" + f for f in fields), "-Ov"],
                            stdin=strip.stdout, stdout=subprocess.PIPE, text=True)
    strip.stdout.close()
    with open(vcf, "wb") as out:
        gz = subprocess.Popen([bgzip, "-c"], stdin=subprocess.PIPE, stdout=out)
        gz.stdin.write(header.encode())
        for line in body.stdout:
            if not line.startswith("#"):
                gz.stdin.write(line.encode())
        gz.stdin.close()
        rc = (body.wait(), strip.wait(), gz.wait())
    if any(rc):
        raise CacheError(f"segment write failed for {src} (exit codes {rc})")
    run([tabix, "-f", "-p", "vcf", str(vcf)])
    info = dict(meta, predictor=predictor, created_utc=now_utc(), n_records=count_records(bcftools, vcf),
                size_bytes=vcf.stat().st_size, sha256=sha256_file(vcf))
    write_json_atomic(out_prefix.with_name(out_prefix.name + ".json"), info)
    return info


def make_segment(predictor: str, annotated: Path, misses: Path, out_prefix: Path, meta: dict,
                 bcftools: str, tabix: str, bgzip: str) -> dict:
    """Segment of the records of `annotated` (model output) whose sites are in `misses`."""
    with tempfile.TemporaryDirectory(dir=out_prefix.parent) as tmp:
        src = Path(tmp) / "computed.vcf.gz"
        run([bcftools, "isec", "-c", "none", "-n=2", "-w1", str(annotated), str(misses), "-Oz", "-o", str(src)])
        return write_segment(predictor, src, out_prefix, meta, bcftools, tabix, bgzip)


# ---------------------------------------------------------------- verify / seed

def query_values(bcftools: str, vcf: Path, predictor: str):
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t" + "\t".join(f"%INFO/{f}" for f in PREDICTORS[predictor]["fields"]) + "\n"
    proc = subprocess.Popen([bcftools, "query", "-f", fmt, str(vcf)], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        yield line.rstrip("\n").split("\t")
    if proc.wait():
        raise CacheError(f"bcftools query failed on {vcf}")


def find_conflicts(rows, limit: int = 20) -> tuple[int, int, list]:
    """Rows sorted by (CHROM, POS): counts duplicated sites and those whose values disagree."""
    dups, conflicts, examples = 0, 0, []
    for _, group in groupby(rows, key=lambda r: (r[0], r[1])):
        seen = {}
        for r in group:
            key = (r[2], r[3])
            if key in seen:
                dups += 1
                if seen[key] != r[4:]:
                    conflicts += 1
                    if len(examples) < limit:
                        examples.append({"site": ":".join(r[:4]), "values": [seen[key], r[4:]]})
            else:
                seen[key] = r[4:]
    return dups, conflicts, examples


def verify(version_dir: Path, predictor: str, bcftools: str) -> dict:
    """Checks that a site present in several segments has the same values in all of them."""
    segments, skipped = list_segments(version_dir / "segments")
    if not segments:
        return {"segments": 0, "skipped": skipped, "duplicates": 0, "conflicts": 0, "examples": []}
    with tempfile.TemporaryDirectory() as tmp:
        merged = Path(tmp) / "all.vcf.gz"
        if len(segments) > 1:
            run([bcftools, "concat", "-a", *map(str, segments), "-Oz", "-o", str(merged)])
        else:
            merged = segments[0]
        dups, conflicts, examples = find_conflicts(query_values(bcftools, merged, predictor))
    return {"segments": len(segments), "skipped": skipped, "duplicates": dups,
            "conflicts": conflicts, "examples": examples}


def run_predictor_vcf(run_dir: Path, predictor: str) -> Path:
    suffix = PREDICTORS[predictor]["output_suffix"]
    hits = sorted((run_dir / "predictors").glob(f"*.{suffix}.vcf.gz")) or sorted(run_dir.glob(f"*.{suffix}.vcf.gz"))
    if len(hits) != 1:
        raise CacheError(f"expected one *.{suffix}.vcf.gz in {run_dir}, found {len(hits)}")
    return hits[0]


def identity_diff(a: dict, b: dict, prefix: str = "") -> list[str]:
    out = []
    for k in sorted(set(a) | set(b)):
        va, vb = a.get(k), b.get(k)
        if isinstance(va, dict) and isinstance(vb, dict):
            out += identity_diff(va, vb, f"{prefix}{k}.")
        elif va != vb:
            out.append(f"{prefix}{k}: {va} -> {vb}")
    return out


def seed(run_dir: Path, predictor: str, cache_dir: Path, into_version: str | None,
         override_reason: str | None, bcftools: str, tabix: str, bgzip: str,
         repo_root: Path | None = None) -> dict:
    """Adds a finished run's predictor output to the cache. Sites already cached must carry
    identical values (else nothing is written); only new sites go into the segment."""
    manifest_path = run_dir / "pipeline_info" / "RUN_MANIFEST.json"
    if not manifest_path.is_file():
        manifest_path = run_dir / "RUN_MANIFEST.json"
    manifest = load_manifest(manifest_path)
    identity = version_identity(manifest, predictor, repo_root)
    run_vid = version_id(identity)
    if into_version and into_version != run_vid:
        vfile = cache_dir / predictor / into_version / "VERSION.json"
        if not vfile.is_file():
            raise CacheError(f"--into-version {into_version}: {vfile} not found")
        if not override_reason:
            diff = identity_diff(json.loads(vfile.read_text())["identity"], identity)
            raise CacheError("run identity differs from --into-version; pass --override-reason. Differences:\n  "
                             + "\n  ".join(diff))
        vdir = vfile.parent
    else:
        vdir = ensure_version_dir(cache_dir, identity)
        into_version = run_vid
    src = run_predictor_vcf(run_dir, predictor)
    run_id = run_dir.name
    segdir = vdir / "segments"
    segments, _ = list_segments(segdir)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        new_sites = src
        overlap = {"shared_sites": 0, "conflicts": 0}
        if segments:
            lookup = build_lookup(bcftools, tabix, segments, tmp / "lookup.vcf.gz")
            run([bcftools, "isec", "-c", "none", "-n=2", "-p", str(tmp / "isec"), str(src), str(lookup)])
            a = list(query_values(bcftools, tmp / "isec" / "0000.vcf", predictor))
            b = list(query_values(bcftools, tmp / "isec" / "0001.vcf", predictor))
            ka, kb = sorted(a), sorted(b)
            bad = [(x, y) for x, y in zip(ka, kb) if x != y]
            overlap = {"shared_sites": len(ka), "conflicts": len(bad) + abs(len(ka) - len(kb))}
            if overlap["conflicts"]:
                raise CacheError(f"{overlap['conflicts']} shared sites differ from the cache, e.g. {bad[:3]}")
            new_sites = tmp / "new.vcf.gz"
            run([bcftools, "isec", "-c", "none", "-C", "-w1", str(src), str(lookup), "-Oz", "-o", str(new_sites)])
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        meta = {"version_id": into_version, "run_id": run_id, "source": str(src),
                "source_manifest": str(manifest_path),
                "source_commit": manifest["launch"]["annotation_pipeline_repo"]["commit"],
                "seeded": True, "run_version_id": run_vid, "override_reason": override_reason,
                "overlap_with_cache": overlap}
        info = write_segment(predictor, new_sites, segdir / f"{run_id}.seed.{stamp}", meta, bcftools, tabix, bgzip)
    return info


# ---------------------------------------------------------------- CLI

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--predictor", required=True, choices=sorted(PREDICTORS))
        sp.add_argument("--bcftools")
        sp.add_argument("--tabix")
        sp.add_argument("--bgzip")

    sp = sub.add_parser("version-id", help="print the version_id; with --cache-dir also create/check VERSION.json")
    common(sp)
    sp.add_argument("--manifest", required=True, type=Path)
    sp.add_argument("--cache-dir", type=Path)

    sp = sub.add_parser("list-segments", help="write the files of the version's valid segments to a list")
    sp.add_argument("--predictor", choices=sorted(PREDICTORS))
    sp.add_argument("--version-dir", required=True, type=Path)
    sp.add_argument("--output", required=True, type=Path)

    sp = sub.add_parser("split", help="split a VCF into cache hits and misses")
    common(sp)
    for a in ("--input", "--segments-dir", "--hits", "--misses", "--report"):
        sp.add_argument(a, required=True, type=Path)

    sp = sub.add_parser("make-segment", help="write the computed misses of a run as a new segment")
    common(sp)
    for a in ("--annotated", "--misses", "--output-prefix"):
        sp.add_argument(a, required=True, type=Path)
    sp.add_argument("--run-id", required=True)
    sp.add_argument("--version-id", required=True)

    sp = sub.add_parser("seed", help="add a finished run's predictor output to the cache")
    common(sp)
    sp.add_argument("--run-dir", required=True, type=Path)
    sp.add_argument("--cache-dir", required=True, type=Path)
    sp.add_argument("--into-version")
    sp.add_argument("--override-reason")

    sp = sub.add_parser("verify", help="check that sites cached more than once agree")
    common(sp)
    sp.add_argument("--version-dir", required=True, type=Path)

    args = p.parse_args(argv)
    try:
        if args.cmd == "version-id":
            identity = version_identity(load_manifest(args.manifest), args.predictor)
            if args.cache_dir:
                ensure_version_dir(args.cache_dir, identity)
            print(version_id(identity))
            return 0
        if args.cmd == "list-segments":
            n, skipped = write_segment_list(args.version_dir, args.output)
            for reason in skipped:
                print(f"annotation_cache: skipped segment {reason}", file=sys.stderr)
            print(n)
            return 0
        bcftools = tool("bcftools", args.bcftools)
        if args.cmd == "verify":
            res = verify(args.version_dir, args.predictor, bcftools)
            print(json.dumps(res, indent=2))
            return 1 if res["conflicts"] else 0
        tabix = tool("tabix", args.tabix)
        if args.cmd == "split":
            res = split(args.predictor, args.input, args.segments_dir, args.hits, args.misses,
                        args.report, bcftools, tabix)
        else:
            bgzip = tool("bgzip", args.bgzip)
            if args.cmd == "make-segment":
                res = make_segment(args.predictor, args.annotated, args.misses, args.output_prefix,
                                   {"version_id": args.version_id, "run_id": args.run_id},
                                   bcftools, tabix, bgzip)
            else:
                res = seed(args.run_dir, args.predictor, args.cache_dir, args.into_version,
                           args.override_reason, bcftools, tabix, bgzip)
        print(json.dumps(res, indent=2))
        return 0
    except CacheError as e:
        print(f"annotation_cache: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
