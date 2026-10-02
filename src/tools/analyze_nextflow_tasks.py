#!/usr/bin/env python3
"""Per-task runtime / node / outcome table from a Nextflow log and its work directories (read-only).

For every "Task completed" line of a .nextflow.log it reads the task's work dir: node (2nd line of .command.log),
runtime (.command.begin -> .exitcode, or -> .command.log mtime when no .exitcode exists), and whether an .exitcode
file exists. Tasks that Nextflow reports with exit '-' ("terminated for an unknown reason") are tasks SGE killed at
their h_rt limit: no .exitcode is written. Prints: outcome counts, the list of such tasks (start, process, node, hours,
bytes of stdout/stderr = did the task produce output before being killed), and runtime by process and node for
successful tasks. Aggregates only; no sample identifiers are printed.
"""
import argparse, collections, datetime, os, re, statistics as st

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--log", required=True, help=".nextflow.log (or a rotated .nextflow.log.N)")
ap.add_argument("--processes", default="PANGOLIN_ANNOTATE,SPLICEAI_ANNOTATE,SPIP_ANNOTATE",
                help="processes for the per-node runtime table")
a = ap.parse_args()

rows = []
for ln in open(a.log, errors="replace"):
    if "Task completed >" not in ln:
        continue
    m = re.search(r"name: (.*?); status: COMPLETED; exit: (\S+);.*?workDir: (\S+?)\s", ln)
    if not m:
        continue
    name, ex, d = m.group(1), m.group(2), m.group(3).rstrip("]")
    proc = re.sub(r"\s*\(.*", "", name).split(":")[-1]
    try:
        node = open(d + "/.command.log", errors="replace").read().split("\n")[1].strip()
    except Exception:
        node = "?"
    try:
        b = os.stat(d + "/.command.begin").st_mtime
        e = os.stat(d + ("/.exitcode" if os.path.exists(d + "/.exitcode") else "/.command.log")).st_mtime
        hours = (e - b) / 3600
        start = datetime.datetime.fromtimestamp(b)
    except Exception:
        hours, start = None, None
    sizes = []
    for f in (".command.out", ".command.err"):
        try: sizes.append(os.stat(d + "/" + f).st_size)
        except Exception: sizes.append(None)
    rows.append(dict(proc=proc, exit=ex, node=node, hours=hours, start=start, out_err=tuple(sizes)))

print(f"completed-task lines: {len(rows)}  exit codes: {dict(collections.Counter(r['exit'] for r in rows))}")
print("\nTasks without an exit status (killed by SGE at the h_rt limit; Nextflow: 'terminated for an unknown reason'):")
c = collections.Counter()
for r in rows:
    if r["exit"] == "-":
        c[(r["proc"], r["node"], None if r["hours"] is None else round(r["hours"], 0))] += 1
for k, v in sorted(c.items(), key=lambda x: -x[1]):
    print(f"  {v:3d} x process={k[0]:24s} node={k[1]:10s} ran={k[2]} h")
lost = sum(r["hours"] for r in rows if r["exit"] == "-" and r["hours"])
print(f"  task-hours spent in these tasks: {lost:.0f}")
print("\nSuccessful tasks, runtime in hours by process and node (n, median, max):")
by = collections.defaultdict(list)
for r in rows:
    if r["exit"] == "0" and r["hours"] is not None and r["proc"] in a.processes.split(","):
        by[(r["proc"], r["node"])].append(r["hours"])
for (p, n), v in sorted(by.items()):
    print(f"  {p:18s} {n:10s} n={len(v):3d} median={st.median(v):5.1f} max={max(v):5.1f}")
print("\nTasks per node (all outcomes) and share without exit status:")
tot = collections.Counter(r["node"] for r in rows); bad = collections.Counter(r["node"] for r in rows if r["exit"] == "-")
for n in sorted(tot):
    print(f"  {n:10s} tasks={tot[n]:4d} no-exit-status={bad[n]:3d} ({100*bad[n]/tot[n]:.0f}%)")
