#!/usr/bin/env python3
"""Per-task, per-process and critical-path tables for annotate_vcf.nf runs (read-only).

Sources per run (repeatable, NAME is a label used in the output):
  --trace NAME PATH   Nextflow execution_trace_*.txt (all fields come from the trace)
  --log   NAME PATH   .nextflow.log of a run without a trace. Submit time comes from the
                      "[hash] (Re-)submitted process" lines, start/end from the "Task completed"
                      lines; requested cpus/h_rt/h_rss/h_vmem from each task's .command.run,
                      node from line 2 of .command.log. %cpu and peak memory are not in the log:
                      pass --qacct to join SGE accounting (cpu, maxvmem) by job id.
  --qacct PATH        TSV written by resources/pull_qacct_snapshot.sh (optional, repeatable)

Outputs in --outdir (TSV): tasks.tsv (one row per attempt), per_process.tsv, critical_path.tsv,
run_summary.tsv. Times are local wall-clock seconds. Only run labels, process names and chunk
tags are written; no sample identifiers beyond the run label.

The critical path is walked backwards from the last task to finish: the predecessor of a task is
the latest-finishing attempt, among the processes it depends on in annotate_vcf.nf (DEPS below) or
a failed attempt of its own process, that finished before the task was submitted. Each step splits
into: gap (predecessor end -> submit; Nextflow polling and queueSize throttling), queue wait
(submit -> start) and runtime (start -> end).
"""
import argparse, collections, csv, datetime as dt, json, math, os, re, sys

# Upstream processes of each process in annotate_vcf.nf (short names, see short_name()).
DEPS = {
    "GENE_SUBSET": ["CHECK_INPUT_ASSEMBLY"],
    "BROAD:CHUNK_VCF": ["CHECK_INPUT_ASSEMBLY"],
    "RESTR:CHUNK_VCF": ["GENE_SUBSET"],
    "VEP_ANNOTATE": ["BROAD:CHUNK_VCF"],
    "BRANCHPOINT_ANNOTATE": ["BROAD:CHUNK_VCF"],
    "SPLICEAI_ANNOTATE": ["RESTR:CHUNK_VCF"],
    "PANGOLIN_ANNOTATE": ["RESTR:CHUNK_VCF"],
    "SPIP_ANNOTATE": ["RESTR:CHUNK_VCF", "BROAD:CHUNK_VCF"],
    "CONCAT_VEP": ["VEP_ANNOTATE"],
    "CONCAT_BRANCHPOINT": ["BRANCHPOINT_ANNOTATE"],
    "CONCAT_SPLICEAI": ["SPLICEAI_ANNOTATE"],
    "CONCAT_PANGOLIN": ["PANGOLIN_ANNOTATE"],
    "CONCAT_SPIP_RESTRICTED": ["SPIP_ANNOTATE"],
    "CONCAT_SPIP_BROAD": ["SPIP_ANNOTATE"],
    "MERGE_ANNOTATIONS": ["CONCAT_VEP", "CONCAT_BRANCHPOINT", "CONCAT_SPLICEAI", "CONCAT_PANGOLIN",
                          "CONCAT_SPIP_RESTRICTED", "CONCAT_SPIP_BROAD"],
    "VCF_TO_TABLE": ["MERGE_ANNOTATIONS"],
    "VCF_TO_TSV": ["MERGE_ANNOTATIONS"],
    "TSV_TO_TABLE": ["VCF_TO_TSV"],
    "EXTRACT_GENOTYPES": ["MERGE_ANNOTATIONS"],
    "TAG_TRANSCRIPT_PRIORITY": ["VCF_TO_TABLE", "TSV_TO_TABLE"],
    "ADD_WIDE_GENOTYPES": ["TAG_TRANSCRIPT_PRIORITY", "EXTRACT_GENOTYPES"],
}
TOL = 120  # s: a task is submitted at most this long after its predecessor is seen complete
MIN_CPU_S = 60  # s: %cpu of shorter tasks is not summarised (sampling noise)


def short_name(proc):
    """CHUNK_VCF keeps its tier; every other process drops the subworkflow prefix."""
    last = proc.split(":")[-1]
    if last == "CHUNK_VCF":
        return ("BROAD:" if proc.startswith("BROAD") else "RESTR:") + last
    return last


def dur_s(txt):
    """Trace duration ('1h 2m 3s', '70ms', '-') -> seconds."""
    if txt in ("", "-"):
        return None
    tot = 0.0
    for num, unit in re.findall(r"([\d.]+)\s*(ms|d|h|m|s)", txt):
        tot += float(num) * {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}[unit]
    return tot


def mem_gb(txt):
    """Trace memory ('12 GB', '443.7 MB', '0', '-') -> GiB."""
    if txt in ("", "-"):
        return None
    m = re.match(r"([\d.]+)\s*([KMGT]?B)?", txt.strip())
    if not m:
        return None
    return float(m.group(1)) * {None: 1 / 2**30, "B": 1 / 2**30, "KB": 1 / 2**20, "MB": 1 / 1024,
                                "GB": 1, "TB": 1024}[m.group(2)]


def ts(txt):
    """Trace timestamp 'YYYY-mm-dd HH:MM:SS.fff' (local) -> epoch-like seconds."""
    if txt in ("", "-"):
        return None
    return dt.datetime.strptime(txt, "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=dt.timezone.utc).timestamp()


def read_trace(run, path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            submit, start, end = ts(r["submit"]), ts(r["start"]), ts(r["complete"])
            pc = r.get("%cpu", "-")
            rows.append(dict(
                run=run, source="trace", job_id=r.get("native_id", ""), process=short_name(r["process"]),
                tag=r.get("tag", ""), attempt=int(r["attempt"]), status=r["status"], exit=r["exit"],
                submit=submit, start=start, end=end, realtime_s=dur_s(r.get("realtime", "-")),
                cpus_req=int(r["cpus"]) if r.get("cpus", "-") not in ("-", "") else None,
                mem_req_gb=mem_gb(r.get("memory", "-")), h_rt_req_s=None, h_vmem_req_gb=None,
                cores_used=float(pc.rstrip("%")) / 100 if pc not in ("-", "") else None,
                peak_rss_gb=mem_gb(r.get("peak_rss", "-")), host=r.get("hostname", "-")))
    return rows


LOG_TIME = re.compile(r"^(\w{3}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3})")
SUBMIT = re.compile(r"\[([0-9a-f]{2}/[0-9a-f]{6})\] (Re-)?[Ss]ubmitted process > (.*)$")
DONE = re.compile(r"Task completed > TaskHandler\[jobId: (\S+); id: (\d+); name: (.*?); status: (\w+); "
                  r"exit: (\S+); error: .*?; workDir: (\S+) started: (\S+); exited: (\S+); \]")


def iso_utc(txt):
    """'2026-09-25T08:39:23.64888Z' (any fraction length) -> naive UTC datetime."""
    m = re.match(r"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})(?:\.(\d+))?", txt)
    frac = (m.group(2) or "0")[:6].ljust(6, "0")
    return dt.datetime.strptime(m.group(1), "%Y-%m-%dT%H:%M:%S").replace(microsecond=int(frac))


def hms(txt):
    h, m, s = (int(x) for x in txt.split(":"))
    return h * 3600 + m * 60 + s


CACHE = {}


def read_command_run(wd):
    if wd in CACHE:
        return dict(CACHE[wd])
    out = dict(cpus_req=None, h_rt_req_s=None, mem_req_gb=None, h_vmem_req_gb=None, host="?")
    try:
        with open(os.path.join(wd, ".command.run"), errors="replace") as fh:
            head = "".join(next(fh) for _ in range(12))
    except (OSError, StopIteration):
        head = ""
    m = re.search(r"-pe smp (\d+)", head)
    out["cpus_req"] = int(m.group(1)) if m else None
    m = re.search(r"h_rt=(\d+:\d+:\d+)", head)
    out["h_rt_req_s"] = hms(m.group(1)) if m else None
    m = re.search(r"h_rss=(\d+)M", head)
    out["mem_req_gb"] = int(m.group(1)) / 1024 if m else None
    m = re.search(r"h_vmem=(\d+)G", head)
    out["h_vmem_req_gb"] = int(m.group(1)) if m else None
    try:
        with open(os.path.join(wd, ".command.log"), errors="replace") as fh:
            out["host"] = fh.read().split("\n")[1].strip() or "?"
    except (OSError, IndexError):
        pass
    CACHE[wd] = dict(out)
    return out


def read_log(run, path):
    year = None
    lines = open(path, errors="replace").read().splitlines()
    submits, dones = {}, []
    for ln in lines:
        mt = LOG_TIME.match(ln)
        if not mt:
            continue
        m = SUBMIT.search(ln)
        if m:
            submits.setdefault(m.group(1), (mt.group(1), bool(m.group(2))))
            continue
        m = DONE.search(ln)
        if m:
            dones.append((mt.group(1), m.groups()))
    # The log has no year and is local time; derive both from the first UTC 'exited' stamp.
    offset = None
    for logt, g in dones:
        if g[7] != "-":
            exited = iso_utc(g[7])
            year = exited.year
            local = dt.datetime.strptime(f"{year}-{logt}", "%Y-%b-%d %H:%M:%S.%f")
            offset = round((local - exited).total_seconds() / 1800) * 1800
            break
    if offset is None:
        sys.exit(f"{path}: no completed task with an exit timestamp")

    def logsec(t):
        return dt.datetime.strptime(f"{year}-{t}", "%Y-%b-%d %H:%M:%S.%f").replace(
            tzinfo=dt.timezone.utc).timestamp()

    rows = []
    for logt, (job, tid, name, status, ex, wd, started, exited) in dones:
        wd = wd.rstrip("]")
        key = "/".join(wd.split("/")[-2:])[:9]
        sub = submits.get(key)
        proc = re.sub(r"\s*\(.*", "", name)
        tag = re.search(r"\((.*)\)", name).group(1) if "(" in name else ""
        start = int(started) / 1000 + offset if started not in ("-", "") else None
        end = (iso_utc(exited).replace(tzinfo=dt.timezone.utc).timestamp()
               + offset) if exited != "-" else logsec(logt)
        r = dict(run=run, source="log", job_id=job, process=short_name(proc), tag=tag, attempt=None,
                 resubmitted=sub[1] if sub else None, status=status, exit=ex,
                 submit=logsec(sub[0]) if sub else None, start=start, end=end,
                 realtime_s=(end - start) if start else None, cores_used=None, peak_rss_gb=None)
        r.update(read_command_run(wd))
        rows.append(r)
    # Attempt from the time limit: every label scales h_rt linearly with task.attempt.
    min_rt = collections.defaultdict(lambda: math.inf)
    for r in rows:
        if r["h_rt_req_s"]:
            min_rt[r["process"]] = min(min_rt[r["process"]], r["h_rt_req_s"])
    for r in rows:
        if r["h_rt_req_s"] and min_rt[r["process"]] < math.inf:
            r["attempt"] = round(r["h_rt_req_s"] / min_rt[r["process"]])
        elif r["resubmitted"] is not None:
            r["attempt"] = 2 if r["resubmitted"] else 1
    return rows


def read_qacct(paths):
    acc = {}
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                acc[r["jobnumber"]] = r
    return acc


def to_gb(txt):
    m = re.match(r"([\d.]+)([KMGT]?)B?", txt or "")
    return float(m.group(1)) * {"": 1 / 2**30, "K": 1 / 2**20, "M": 1 / 1024, "G": 1, "T": 1024}[m.group(2)] if m else None


def pct(v, q):
    if not v:
        return None
    v = sorted(v)
    k = (len(v) - 1) * q
    lo, hi = math.floor(k), math.ceil(k)
    return v[lo] + (v[hi] - v[lo]) * (k - lo)


def ok(r):
    return r["exit"] == "0"


def per_process(rows):
    out = []
    by = collections.defaultdict(list)
    for r in rows:
        by[(r["run"], r["process"])].append(r)
    for (run, proc), rs in sorted(by.items(), key=lambda x: (x[0][0], min(r["submit"] or 0 for r in x[1]))):
        good = [r for r in rs if ok(r)]
        rt = [r["realtime_s"] for r in good if r["realtime_s"] is not None]
        qw = [r["start"] - r["submit"] for r in rs if r["start"] and r["submit"]]
        cu = [r["cores_used"] for r in good if r["cores_used"] is not None and (r["realtime_s"] or 0) >= MIN_CPU_S]
        rss = [r["peak_rss_gb"] for r in good if r["peak_rss_gb"] is not None]
        first = [r for r in good if r["attempt"] == 1]
        hosts = collections.Counter(r["host"] for r in rs)
        task_h = sum((r["end"] - r["start"]) for r in rs if r["start"] and r["end"]) / 3600
        lost_h = sum((r["end"] - r["start"]) for r in rs if r["start"] and r["end"] and not ok(r)) / 3600
        cpus = collections.Counter(r["cpus_req"] for r in rs)
        mem = collections.Counter(r["mem_req_gb"] for r in rs)
        out.append(dict(
            run=run, process=proc, tasks_ok=len(good), attempts=len(rs), failed_attempts=len(rs) - len(good),
            ok_on_attempt1=len(first),
            exit_codes=";".join(f"{k}:{v}" for k, v in sorted(collections.Counter(r["exit"] for r in rs).items())),
            realtime_median_h=pct(rt, .5) and pct(rt, .5) / 3600, realtime_p95_h=pct(rt, .95) and pct(rt, .95) / 3600,
            realtime_max_h=max(rt) / 3600 if rt else None,
            cpus_req=";".join(f"{k}:{v}" for k, v in sorted(cpus.items(), key=lambda x: str(x[0]))),
            cores_used_median=pct(cu, .5), cores_used_max=max(cu) if cu else None,
            mem_req_gb=";".join(f"{k:g}:{v}" if k else f"-:{v}" for k, v in sorted(mem.items(), key=lambda x: x[0] or 0)),
            peak_rss_median_gb=pct(rss, .5), peak_rss_max_gb=max(rss) if rss else None,
            queue_wait_median_min=pct(qw, .5) and pct(qw, .5) / 60, queue_wait_p95_min=pct(qw, .95) and pct(qw, .95) / 60,
            queue_wait_max_min=max(qw) / 60 if qw else None,
            first_submit=min((r["submit"] for r in rs if r["submit"]), default=None),
            last_end=max((r["end"] for r in rs if r["end"]), default=None),
            task_hours=task_h, failed_task_hours=lost_h,
            hosts=";".join(f"{h}:{n}" for h, n in hosts.most_common())))
    return out


def critical_path(rows):
    rs = [r for r in rows if r["end"] is not None and r["submit"] is not None]
    cur = max(rs, key=lambda r: r["end"])
    path = []
    while cur is not None:
        cands = [r for r in rs if r["end"] < cur["end"] and r["end"] <= cur["submit"] + TOL and (
            (r["process"] in DEPS.get(cur["process"], []) and ok(r))
            or (r["process"] == cur["process"] and not ok(r)))]
        pred = max(cands, key=lambda r: r["end"]) if cands else None
        path.append((cur, pred))
        cur = pred
    path.reverse()
    out = []
    for cur, pred in path:
        out.append(dict(
            run=cur["run"], process=cur["process"], tag=cur["tag"], attempt=cur["attempt"], exit=cur["exit"],
            host=cur["host"], gap_min=(max(0, cur["submit"] - pred["end"]) / 60) if pred else 0,
            queue_wait_min=(cur["start"] - cur["submit"]) / 60 if cur["start"] else None,
            runtime_h=(cur["end"] - cur["start"]) / 3600 if cur["start"] else None,
            end_h_since_launch=None))
    t0 = min(r["submit"] for r in rs)
    for o, (cur, _) in zip(out, path):
        o["end_h_since_launch"] = (cur["end"] - t0) / 3600
    return out


def run_summary(rows, cp, procs):
    out = []
    for run in dict.fromkeys(r["run"] for r in rows):
        rs = [r for r in rows if r["run"] == run]
        t0 = min(r["submit"] for r in rs if r["submit"])
        t1 = max(r["end"] for r in rs if r["end"])
        c = [x for x in cp if x["run"] == run]
        pp = {p["process"]: p for p in procs if p["run"] == run}
        # Tasks inside SGE (submitted, not yet finished) over time: Nextflow caps this at queueSize.
        ev = sorted([(r["submit"], 1) for r in rs if r["submit"]] + [(r["end"], -1) for r in rs if r["end"]])
        n = peak = 0
        at_peak = 0.0
        prev = None
        for t, d in ev:
            if prev is not None and n >= 95:
                at_peak += t - prev
            n += d
            peak = max(peak, n)
            prev = t

        def end_of(p):
            return (pp[p]["last_end"] - t0) / 3600 if p in pp else None
        out.append(dict(
            run=run, wall_h=(t1 - t0) / 3600, attempts=len(rs), failed_attempts=sum(not ok(r) for r in rs),
            task_hours=sum((r["end"] - r["start"]) for r in rs if r["start"] and r["end"]) / 3600,
            failed_task_hours=sum((r["end"] - r["start"]) for r in rs if r["start"] and r["end"] and not ok(r)) / 3600,
            cp_runtime_h=sum(x["runtime_h"] or 0 for x in c if x["exit"] == "0"),
            cp_failed_attempt_h=sum(x["runtime_h"] or 0 for x in c if x["exit"] != "0"),
            cp_queue_wait_h=sum(x["queue_wait_min"] or 0 for x in c) / 60,
            cp_gap_h=sum(x["gap_min"] for x in c) / 60,
            max_tasks_in_sge=peak, hours_with_ge95_in_sge=at_peak / 3600,
            concat_vep_end_h=end_of("CONCAT_VEP"), concat_branchpoint_end_h=end_of("CONCAT_BRANCHPOINT"),
            concat_spliceai_end_h=end_of("CONCAT_SPLICEAI"), concat_pangolin_end_h=end_of("CONCAT_PANGOLIN"),
            concat_spip_end_h=end_of("CONCAT_SPIP_RESTRICTED"), merge_end_h=end_of("MERGE_ANNOTATIONS"),
            tag_end_h=end_of("TAG_TRANSCRIPT_PRIORITY")))
    return out


def write(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.3f}" if isinstance(v, float) else ("-" if v is None else v)) for k, v in r.items()})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--trace", nargs=2, action="append", default=[], metavar=("NAME", "PATH"))
    ap.add_argument("--log", nargs=2, action="append", default=[], metavar=("NAME", "PATH"))
    ap.add_argument("--qacct", action="append", default=[])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cache", help="JSON cache of the per-work-dir metadata (speeds up re-runs on NFS)")
    a = ap.parse_args()
    if a.cache and os.path.exists(a.cache):
        CACHE.update(json.load(open(a.cache)))
    rows = []
    for run, p in a.trace:
        rows += read_trace(run, p)
    for run, p in a.log:
        rows += read_log(run, p)
    acc = read_qacct(a.qacct)
    for r in rows:
        q = acc.get(str(r["job_id"]))
        if q and r["cores_used"] is None:
            wall = float(q["ru_wallclock"].rstrip("s"))
            r["cores_used"] = float(q["cpu"].rstrip("s")) / wall if wall > 0 else None
            r["peak_rss_gb"] = to_gb(q["maxvmem"])  # maxvmem: virtual, an upper bound of RSS
            r["source"] += "+qacct"
    if a.cache:
        json.dump(CACHE, open(a.cache, "w"))
    rows.sort(key=lambda r: (r["run"], r["submit"] or 0))
    for r in rows:
        r["queue_wait_s"] = (r["start"] - r["submit"]) if r["start"] and r["submit"] else None
    os.makedirs(a.outdir, exist_ok=True)
    cols = ["run", "source", "job_id", "process", "tag", "attempt", "status", "exit", "host", "submit", "start",
            "end", "queue_wait_s", "realtime_s", "cpus_req", "cores_used", "mem_req_gb", "peak_rss_gb",
            "h_rt_req_s", "h_vmem_req_gb"]
    write(os.path.join(a.outdir, "tasks.tsv"), [{c: r.get(c) for c in cols} for r in rows])
    procs = per_process(rows)
    write(os.path.join(a.outdir, "per_process.tsv"), procs)
    cp = []
    for run in dict.fromkeys(r["run"] for r in rows):
        cp += critical_path([r for r in rows if r["run"] == run])
    write(os.path.join(a.outdir, "critical_path.tsv"), cp)
    write(os.path.join(a.outdir, "run_summary.tsv"), run_summary(rows, cp, procs))


if __name__ == "__main__":
    main()
