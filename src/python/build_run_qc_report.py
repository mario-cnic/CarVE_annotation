#!/usr/bin/env python3
"""Run QC report for one annotate_vcf.nf run folder: provenance, record counts per stage, time and memory per process, failed attempts.

Reads the run folder (sub-folder layout or the older flat layout), its manifest, the Nextflow trace written
to pipeline_info/ (or, when there is none, a .nextflow.log) and the published VCFs / Parquet tables (index and
column-limited reads only). Writes a timestamped HTML, JSON and per-process TSV into <run>/reports/; never
overwrites an existing report. A value that cannot be determined is reported as unavailable, never as 0.

Usage:
    build_run_qc_report.py --run-dir nf_work/annotation_out/<run_id> [--nextflow-log .nextflow.log]
        [--bcftools PATH] [--out-dir DIR] [--skip-distinct-locus]
"""
from __future__ import annotations

import argparse
import html
import json
import re
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from write_run_manifest import INFO_SUBDIR, REPORTS_SUBDIR, run_files  # noqa: E402

PREDICTOR_VCFS = ["annVEP", "annBranchpoint", "annSpliceAI", "annPangolin", "annSPiP"]
ERROR_TAIL_BYTES = 4000
UNAVAILABLE = None  # a count or time that could not be determined (distinct from 0)


# ----------------------------------------------------------------------------- inputs

def load_manifest(run_dir: Path) -> dict | None:
    found = run_files(run_dir, INFO_SUBDIR, "RUN_MANIFEST.json")
    return json.loads(found[0].read_text()) if found else None


def to_ms(text: str) -> int | None:
    """Trace durations as milliseconds: raw integers or Nextflow's '1d 2h 3m 4.5s' / '343ms' strings."""
    text = text.strip()
    if text in ("", "-"):
        return UNAVAILABLE
    if re.fullmatch(r"\d+(\.\d+)?", text):
        return int(float(text))
    units = {"d": 86400000, "h": 3600000, "m": 60000, "s": 1000, "ms": 1}
    parts = re.findall(r"([\d.]+)\s*(ms|d|h|m|s)", text)
    return int(sum(float(v) * units[u] for v, u in parts)) if parts else UNAVAILABLE


def to_bytes(text: str) -> int | None:
    """Trace memory as bytes: raw integers or '1.5 GB' / '14.3 MB' strings."""
    text = text.strip()
    if text in ("", "-"):
        return UNAVAILABLE
    if re.fullmatch(r"\d+(\.\d+)?", text):
        return int(float(text))
    m = re.fullmatch(r"([\d.]+)\s*(B|KB|MB|GB|TB)", text)
    return int(float(m.group(1)) * 1024 ** ["B", "KB", "MB", "GB", "TB"].index(m.group(2))) if m else UNAVAILABLE


def process_of(name: str) -> str:
    """'WF:SUB:PROCESS (tag)' -> 'PROCESS'."""
    return re.sub(r"\s*\(.*", "", name).split(":")[-1]


def read_trace(path: Path) -> list[dict]:
    lines = path.read_text(errors="replace").splitlines()
    header = lines[0].split("\t")
    rows = []
    for ln in lines[1:]:
        r = dict(zip(header, ln.split("\t")))
        rows.append({
            "name": r.get("name", ""), "process": process_of(r.get("name", "")), "status": r.get("status", ""),
            "exit": r.get("exit", ""), "attempt": int(r["attempt"]) if r.get("attempt", "").isdigit() else UNAVAILABLE,
            "duration_ms": to_ms(r.get("realtime") or r.get("duration", "")),
            "requested_bytes": to_bytes(r.get("memory", "")), "peak_rss_bytes": to_bytes(r.get("peak_rss", "")),
            "host": (r.get("hostname") or "").strip() if (r.get("hostname") or "-").strip() != "-" else UNAVAILABLE,
            "workdir": r.get("workdir", ""),
        })
    return rows


def read_log(path: Path) -> list[dict]:
    """Task rows from the 'Task completed' lines of a .nextflow.log (no memory, no hostname field)."""
    pattern = re.compile(r"name: (.*?); status: (\w+); exit: (\S+);.*?workDir: (\S+?)\s")
    seen, rows = Counter(), []
    for ln in path.read_text(errors="replace").splitlines():
        if "Task completed >" not in ln:
            continue
        m = pattern.search(ln)
        if not m:
            continue
        name, status, exit_code, workdir = m.groups()
        seen[name] += 1
        duration = workdir_runtime_ms(Path(workdir))
        rows.append({
            "name": name, "process": process_of(name), "status": status, "exit": exit_code, "attempt": seen[name],
            "duration_ms": duration, "requested_bytes": UNAVAILABLE, "peak_rss_bytes": UNAVAILABLE,
            "host": log_node(Path(workdir)), "workdir": workdir,
        })
    return rows


def workdir_runtime_ms(workdir: Path) -> int | None:
    """Runtime from the work dir: .command.begin to .exitcode (or to .command.log when SGE killed the task before one was written)."""
    try:
        begin = workdir.joinpath(".command.begin").stat().st_mtime
        end_file = workdir / ".exitcode" if workdir.joinpath(".exitcode").exists() else workdir / ".command.log"
        return int((end_file.stat().st_mtime - begin) * 1000)
    except OSError:
        return UNAVAILABLE


def log_node(workdir: Path) -> str | None:
    """Node name: second line of the task's .command.log (SGE wrapper), when the work dir still exists."""
    try:
        return workdir.joinpath(".command.log").read_text(errors="replace").split("\n")[1].strip() or UNAVAILABLE
    except (OSError, IndexError):
        return UNAVAILABLE


# ----------------------------------------------------------------------------- summaries

def process_table(rows: list[dict]) -> list[dict]:
    by = defaultdict(list)
    for r in rows:
        by[r["process"]].append(r)
    out = []
    for proc, rs in sorted(by.items()):
        ok = [r for r in rs if r["status"] == "COMPLETED" and r["exit"] == "0"]
        failed = [r for r in rs if r["status"] == "FAILED" or (r["status"] == "COMPLETED" and r["exit"] != "0")]
        cached = [r for r in rs if r["status"] == "CACHED"]
        secs = [r["duration_ms"] / 1000 for r in ok if r["duration_ms"] is not None]
        failed_secs = [r["duration_ms"] / 1000 for r in failed if r["duration_ms"] is not None]
        peaks = [r["peak_rss_bytes"] for r in ok if r["peak_rss_bytes"] is not None]
        req = [r["requested_bytes"] for r in ok if r["requested_bytes"] is not None]
        hosts = sorted({r["host"] for r in rs if r["host"]})
        out.append({
            "process": proc, "tasks_ok": len(ok), "tasks_cached": len(cached), "failed_attempts": len(failed),
            "total_h": round(sum(secs) / 3600, 3) if secs else UNAVAILABLE,
            "median_h": round(statistics.median(secs) / 3600, 3) if secs else UNAVAILABLE,
            "max_h": round(max(secs) / 3600, 3) if secs else UNAVAILABLE,
            "failed_attempts_h": round(sum(failed_secs) / 3600, 3) if failed_secs else UNAVAILABLE,
            "requested_mem_gb_median": round(statistics.median(req) / 1024 ** 3, 2) if req else UNAVAILABLE,
            "peak_rss_gb_max": round(max(peaks) / 1024 ** 3, 2) if peaks else UNAVAILABLE,
            "nodes": ",".join(hosts) if hosts else UNAVAILABLE,
        })
    return out


def last_lines(path: Path, n: int = 1) -> str | None:
    try:
        with path.open("rb") as f:
            f.seek(0, 2)
            f.seek(max(0, f.tell() - ERROR_TAIL_BYTES))
            text = f.read().decode(errors="replace")
    except OSError:
        return UNAVAILABLE
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    return " | ".join(lines[-n:]) if lines else None


def failed_attempts(rows: list[dict]) -> list[dict]:
    """Failed or retried attempts with exit code, node, and the last stderr line (else the last .command.log line)."""
    out = []
    for r in rows:
        if r["status"] == "CACHED" or (r["status"] == "COMPLETED" and r["exit"] == "0"):
            continue
        wd = Path(r["workdir"]) if r["workdir"] else None
        if wd and wd.is_dir():
            err = last_lines(wd / ".command.err") or last_lines(wd / ".command.log", 2) or "no output in .command.err / .command.log"
        else:
            err = UNAVAILABLE
        out.append({"name": r["name"], "attempt": r["attempt"], "exit": r["exit"], "host": r["host"],
                    "hours": round(r["duration_ms"] / 3.6e6, 3) if r["duration_ms"] is not None else UNAVAILABLE,
                    "last_error": err})
    return out


def vcf_records(path: Path, bcftools: str) -> int | None:
    """Record count from the tabix index (no VCF scan)."""
    try:
        res = subprocess.run([bcftools, "index", "-n", str(path)], capture_output=True, text=True, timeout=300)
        return int(res.stdout.strip()) if res.returncode == 0 else UNAVAILABLE
    except (OSError, ValueError, subprocess.TimeoutExpired):
        return UNAVAILABLE


def table_counts(table: Path, distinct_locus: bool) -> dict:
    """Rows from Parquet metadata; status / tier / Locus from column-limited reads."""
    try:
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
    except ImportError:
        return {"note": "pyarrow not available in this Python"}
    try:
        pf = pq.ParquetFile(table)
    except Exception as e:  # unreadable or empty file: report it, keep the rest of the report
        return {"note": f"could not read {table.name}: {str(e).splitlines()[0]}"}
    names = pf.schema_arrow.names
    out = {"rows": pf.metadata.num_rows, "columns": len(names), "status": {}, "tier": None, "distinct_locus": None}
    for col in [c for c in names if c.endswith("_status")]:
        counts = pc.value_counts(pf.read(columns=[col]).column(col).combine_chunks())
        out["status"][col] = {str(s["values"]): int(s["counts"]) for s in counts.to_pylist()}
    if "TRANSCRIPT_PRIORITY_TIER" in names:
        counts = pc.value_counts(pf.read(columns=["TRANSCRIPT_PRIORITY_TIER"]).column(0).combine_chunks())
        out["tier"] = {str(s["values"]): int(s["counts"]) for s in counts.to_pylist()}
    if distinct_locus and "Locus" in names:
        out["distinct_locus"] = int(pc.count_distinct(pf.read(columns=["Locus"]).column(0)).as_py())
    return out


def read_check(run_dir: Path) -> dict:
    found = run_files(run_dir, REPORTS_SUBDIR, "*.assembly_check.tsv")
    if not found:
        return {}
    rows = (ln.partition("\t") for ln in found[0].read_text().splitlines()[1:])
    return {k: v for k, _, v in rows}


def find_one(run_dir: Path, subdir: str, suffix: str) -> Path | None:
    found = run_files(run_dir, subdir, f"*{suffix}")
    return found[0] if found else None


# ----------------------------------------------------------------------------- report

def build(run_dir: Path, nextflow_log: Path | None, bcftools: str, distinct_locus: bool) -> dict:
    manifest = load_manifest(run_dir)
    launch, done = (manifest or {}).get("launch", {}), ((manifest or {}).get("completed") or {})
    check = read_check(run_dir)
    notes = []

    traces = sorted((run_dir / INFO_SUBDIR).glob("execution_trace_*.txt"))
    if traces:
        task_rows, task_source = read_trace(traces[-1]), f"trace {traces[-1].name}" + (
            f" (latest of {len(traces)}; cached tasks of a resumed run are listed as CACHED)" if len(traces) > 1 else "")
    elif nextflow_log and nextflow_log.is_file():
        task_rows, task_source = read_log(nextflow_log), f"log {nextflow_log.name} (no trace: runtime and node from the work dirs, no memory figures)"
        notes.append("task data parsed from a Nextflow log; peak memory is not available")
    else:
        task_rows, task_source = [], None
        notes.append("no execution trace and no Nextflow log: per-process times and failed attempts are unavailable")

    run_id = launch.get("cli_flags", {}).get("run_id") or run_dir.name
    wall_h = None
    if launch.get("generated_at") and done.get("generated_at"):
        wall_h = round((datetime.fromisoformat(done["generated_at"]) - datetime.fromisoformat(launch["generated_at"])).total_seconds() / 3600, 2)

    # Record counts per stage, each from the published file when present.
    stages = {}
    in_count = check.get("records_scanned")
    stages["input records"] = int(in_count) if in_count else UNAVAILABLE
    ann = find_one(run_dir, "annotation", ".annotated.vcf.gz")
    pred = {s: find_one(run_dir, "predictors", f".{s}.vcf.gz") for s in PREDICTOR_VCFS}
    for label, path in [(f"{s} records", p) for s, p in pred.items()] + [("annotated (merged) records", ann)]:
        stages[label] = vcf_records(path, bcftools) if path else UNAVAILABLE
    table_path = next((p for p in run_files(run_dir, "tables", "*.parsed.clean.pq") if ".withGT." not in p.name), None)
    table = table_counts(table_path, distinct_locus) if table_path else {"note": "no parsed.clean.pq found (TSV/XLSX output is not summarised)"}
    gt_path = find_one(run_dir, "genotypes", ".genotypes.pq")
    if gt_path:
        try:
            import pyarrow.parquet as pq
            gt_rows, gt_note = pq.ParquetFile(gt_path).metadata.num_rows, None
        except Exception as e:  # pyarrow missing, or unreadable file
            gt_rows, gt_note = UNAVAILABLE, f"could not read {gt_path.name}: {str(e).splitlines()[0]}"
    else:
        gt_rows = UNAVAILABLE
        gt_note = ("not applicable: input has no genotypes" if check.get("genotype_mode") in ("no_samples", "samples_without_GT")
                   else "genotype table not found")

    return {
        "run_id": run_id, "generated_at": datetime.now(timezone.utc).isoformat(), "run_dir": str(run_dir),
        "summary": {
            "commit": launch.get("annotation_pipeline_repo", {}).get("commit"), "dirty": launch.get("annotation_pipeline_repo", {}).get("dirty"),
            "nextflow_version": launch.get("nextflow_version"), "input": launch.get("resolved_paths", {}).get("input"),
            "input_sha256": launch.get("file_hashes_sha256", {}).get("input"), "genotype_mode": check.get("genotype_mode"),
            "samples": check.get("sample_names"), "launched_at": launch.get("generated_at"), "completed_at": done.get("generated_at"),
            "wall_hours": wall_h, "exit_code": done.get("exit_code"), "workflow_stats": done.get("workflow_stats"),
        },
        "stage_counts": stages, "table": table, "genotype_rows": gt_rows, "genotype_note": gt_note,
        "task_source": task_source, "processes": process_table(task_rows), "failed_attempts": failed_attempts(task_rows),
        "provenance": {k: launch.get(k) for k in ("spliceai_distance", "pangolin_distance", "tool_versions", "launched_by", "host")},
        "notes": notes if manifest else notes + ["no RUN_MANIFEST.json: summary and provenance are unavailable"],
    }


def fmt(v) -> str:
    return "unavailable" if v is None else (f"{v:,}" if isinstance(v, int) and not isinstance(v, bool) else str(v))


def table_html(headers: list[str], rows: list[list]) -> str:
    head = "".join(f"<th>{html.escape(h)}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(fmt(c))}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def render_html(rep: dict) -> str:
    s, t = rep["summary"], rep["table"]
    parts = [f"<h1>Run QC: {html.escape(rep['run_id'])}</h1><p class=meta>generated {html.escape(rep['generated_at'])}</p>"]
    if rep["notes"]:
        parts.append("<ul class=notes>" + "".join(f"<li>{html.escape(n)}</li>" for n in rep["notes"]) + "</ul>")
    parts.append("<h2>1. Run summary</h2>" + table_html(["field", "value"], [[k, v] for k, v in s.items()]))
    stage_rows = [[k, v] for k, v in rep["stage_counts"].items()]
    stage_rows.append(["table rows (variant x transcript)", t.get("rows")])
    stage_rows.append(["table distinct Locus", t.get("distinct_locus")])
    stage_rows.append(["genotype rows" + (f" ({rep['genotype_note']})" if rep["genotype_note"] else ""), rep["genotype_rows"]])
    parts.append("<h2>2. Records per stage</h2>" + table_html(["stage", "count"], stage_rows))
    if t.get("note"):
        parts.append(f"<p>{html.escape(t['note'])}</p>")
    if t.get("status"):
        keys = sorted({k for v in t["status"].values() for k in v})
        parts.append("<h3>Predictor status (table rows)</h3>" + table_html(
            ["column"] + keys, [[c] + [v.get(k, 0) for k in keys] for c, v in t["status"].items()]))
    if t.get("tier"):
        parts.append("<h3>TRANSCRIPT_PRIORITY_TIER (table rows)</h3>" + table_html(["tier", "rows"], sorted(t["tier"].items())))
    cols = ["process", "tasks_ok", "tasks_cached", "failed_attempts", "total_h", "median_h", "max_h", "failed_attempts_h",
            "requested_mem_gb_median", "peak_rss_gb_max", "nodes"]
    parts.append(f"<h2>3. Time and memory per process</h2><p>source: {html.escape(fmt(rep['task_source']))}; "
                 "durations are successful attempts, hours; failed attempts counted separately.</p>"
                 + table_html(cols, [[p[c] for c in cols] for p in rep["processes"]]))
    fa = rep["failed_attempts"]
    parts.append(f"<h2>4. Failed or retried attempts ({len(fa)})</h2>" + (table_html(
        ["name", "attempt", "exit", "host", "hours", "last_error"], [[f[c] for c in ("name", "attempt", "exit", "host", "hours", "last_error")] for f in fa])
        if fa else "<p>none recorded</p>"))
    parts.append("<h2>5. Provenance</h2><pre>" + html.escape(json.dumps(rep["provenance"], indent=2)) + "</pre>")
    css = ("body{font:14px system-ui;margin:2rem auto;max-width:1100px;padding:0 1rem;background:#fff;color:#222}"
           "table{border-collapse:collapse;margin:.5rem 0}td,th{border:1px solid #bbb;padding:2px 8px;text-align:left;vertical-align:top}"
           "th{background:#eee}.meta,.notes{color:#666}@media(prefers-color-scheme:dark){body{background:#181818;color:#ddd}"
           "td,th{border-color:#555}th{background:#2a2a2a}}")
    return f"<!doctype html><meta charset=utf-8><title>Run QC {html.escape(rep['run_id'])}</title><style>{css}</style>" + "".join(parts)


def write_reports(rep: dict, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    base = out_dir / f"{rep['run_id']}.run_qc.{stamp}"
    cols = ["process", "tasks_ok", "tasks_cached", "failed_attempts", "total_h", "median_h", "max_h", "failed_attempts_h",
            "requested_mem_gb_median", "peak_rss_gb_max", "nodes"]
    tsv = "\t".join(cols) + "\n" + "".join("\t".join("NA" if p[c] is None else str(p[c]) for c in cols) + "\n" for p in rep["processes"])
    files = {base.with_suffix(".html"): render_html(rep), base.with_suffix(".json"): json.dumps(rep, indent=2) + "\n",
             base.with_suffix(".processes.tsv"): tsv}
    for path, text in files.items():
        with path.open("x") as f:  # 'x': never overwrite an existing report
            f.write(text)
    return list(files)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--nextflow-log", type=Path, help="used only when the run has no execution trace")
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--out-dir", type=Path, help="default: <run-dir>/reports")
    ap.add_argument("--skip-distinct-locus", action="store_true", help="skip the distinct-Locus count (reads the Locus column)")
    a = ap.parse_args()
    if not a.run_dir.is_dir():
        sys.exit(f"error: run folder not found: {a.run_dir}")
    rep = build(a.run_dir, a.nextflow_log, a.bcftools, not a.skip_distinct_locus)
    for p in write_reports(rep, a.out_dir or a.run_dir / REPORTS_SUBDIR):
        print(f"wrote {p}")


if __name__ == "__main__":
    main()
