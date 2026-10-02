#!/usr/bin/env python3
"""MISC-17: compare qacct accounting of healthy vs killed SpliceAI/Pangolin chunks.

Inputs (written by resources/pull_qacct_snapshot.sh on the login node):
  --healthy  misc17_healthy_qacct.tsv  (step label: <PROC>_ok_<node>_<hours>h)
  --killed   misc16_17_qacct.tsv       (step label: <PROC>_KILLED_at_h_rt; the
                                        killed rows carry NO node: qacct's
                                        `hostname` field was not pulled)
Units: cpu = SGE summed CPU seconds over all slots; ru_wallclock = seconds;
cores_busy = cpu / ru_wallclock; cpu_eff = cores_busy / slots.
Read-only; prints tables to stdout.
"""
import argparse
import re

import pandas as pd


def to_gb(s):
    m = re.fullmatch(r"([\d.]+)(B|MB|GB)", str(s))
    return float(m.group(1)) / {"B": 1e9, "MB": 1e3, "GB": 1}[m.group(2)]


def load(path):
    d = pd.read_csv(path, sep="\t")
    d["cpu_s"] = d["cpu"].str.rstrip("s").astype(float)
    d["wall_s"] = d["ru_wallclock"].astype(str).str.rstrip("s").astype(float)
    d["wall_h"] = d["wall_s"] / 3600
    d["cpu_h"] = d["cpu_s"] / 3600
    d["cores_busy"] = d["cpu_s"] / d["wall_s"]
    d["cpu_eff"] = d["cores_busy"] / d["slots"]
    d["maxvmem_gb"] = d["maxvmem"].map(to_gb)
    d["proc"] = d["step"].str.extract(r"^(SPLICEAI|PANGOLIN)_ANNOTATE")[0]
    return d


def job_nodes(nxf_log):
    """SGE job id -> (node, submit 'Mon-DD HH:MM') from the Nextflow log and the task work dirs."""
    out = {}
    pat = re.compile(r"^(\w{3}-\d\d \d\d:\d\d):\S+ .*submitted process .* > jobId: (\d+); workDir: (\S+)")
    with open(nxf_log, errors="replace") as fh:
        for line in fh:
            m = pat.match(line)
            if not m:
                continue
            try:
                node = open(m.group(3) + "/.command.log", errors="replace").read().split("\n")[1].strip()
            except (OSError, IndexError):
                node = "?"
            out[int(m.group(2))] = (node, m.group(1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--healthy", required=True)
    ap.add_argument("--killed", required=True)
    ap.add_argument("--nxf-log", help="Nextflow log of the same run (e.g. .nextflow.log.5): maps SGE job id -> work dir "
                    "-> node (2nd line of .command.log, read-only), so killed rows get a node")
    a = ap.parse_args()

    h = load(a.healthy)
    h["node"] = h["step"].str.extract(r"_ok_(c\d+-cn\d+)_")[0]
    tail = h["wall_h"] < 0.1  # final, tiny chunks (96-115 s): not comparable to full 8,000-variant chunks
    print(f"excluded {int(tail.sum())} tail chunks with wall < 0.1 h from all healthy statistics\n")
    h = h[~tail].copy()
    k = load(a.killed)
    k = k[k["step"].str.contains("KILLED_at_h_rt")].copy()
    k["node"] = "unknown"
    jn = job_nodes(a.nxf_log) if a.nxf_log else {}
    if jn:
        k["node"] = k["jobnumber"].map(lambda j: jn.get(int(j), ("unknown", ""))[0])
        k["submit"] = k["jobnumber"].map(lambda j: jn.get(int(j), ("", ""))[1])
        h["log_node"] = h["jobnumber"].map(lambda j: jn.get(int(j), ("unknown", ""))[0])
        h["submit"] = h["jobnumber"].map(lambda j: jn.get(int(j), ("", ""))[1])
        bad = h[h.node != h.log_node]
        print(f"healthy label node vs Nextflow-log node: {len(h) - len(bad)} agree, {len(bad)} disagree/missing")
    k["wall_cls"] = k["wall_h"].round().astype(int).astype(str) + "h"

    pd.set_option("display.width", 200, "display.float_format", "{:.2f}".format)
    cols = ["cpu_h", "wall_h", "cores_busy", "cpu_eff", "maxvmem_gb"]

    print("== healthy: n, slots, per process ==")
    print(h.groupby(["proc", "slots"]).size().unstack(fill_value=0))
    print("\n== healthy: per process, median [min-max] ==")
    for p, g in h.groupby("proc"):
        print(p, len(g))
        print(g[cols].describe().loc[["min", "50%", "max"]])
    print("\n== healthy: per process x node (median; n) ==")
    t = h.groupby(["proc", "node"])[cols].median()
    t["n"] = h.groupby(["proc", "node"]).size()
    t["slots_set"] = h.groupby(["proc", "node"])["slots"].agg(lambda s: ",".join(map(str, sorted(set(s)))))
    print(t)
    print("\n== healthy: same, only 4-slot jobs ==")
    h4 = h[h.slots == 4]
    t = h4.groupby(["proc", "node"])[cols].median()
    t["n"] = h4.groupby(["proc", "node"]).size()
    print(t)
    print("\n== healthy 8-slot jobs ==")
    print(h[h.slots == 8][["step", "slots", "cpu_h", "wall_h", "cores_busy", "cpu_eff"]])
    print("\n== killed (node unknown) ==")
    print(k.groupby(["proc", "slots", "wall_cls"])[cols].median().join(
        k.groupby(["proc", "slots", "wall_cls"]).size().rename("n")))
    if jn:
        print("\n== killed per node (node from Nextflow log) ==")
        print(k.groupby(["proc", "node", "slots", "wall_cls"]).agg(
            n=("cpu_h", "size"), cpu_h=("cpu_h", "median"), cores_busy=("cores_busy", "median")))
        print("\n== submit time (Mon-DD HH:MM) of every S223 SpliceAI/Pangolin chunk, by node and outcome ==")
        a1 = h[["proc", "node", "slots", "wall_h", "cpu_h", "submit"]].assign(outcome="ok")
        a2 = k[["proc", "node", "slots", "wall_h", "cpu_h", "submit"]].assign(outcome="KILLED")
        allr = pd.concat([a1, a2]).sort_values(["node", "submit"])
        print(allr[allr.node.isin(["c0053-cn1", "c0051-cn1", "c0051-cn4"])].to_string(index=False))
    print("\n== killed: cpu_h spread per process ==")
    print(k.groupby("proc")["cpu_h"].describe())

    # CPU-hours needed for a full chunk: healthy vs killed (killed did NOT finish)
    print("\n== CPU-hours per chunk: healthy (finished) vs killed (unfinished) ==")
    for p in ["SPLICEAI", "PANGOLIN"]:
        hh, kk = h[h.proc == p]["cpu_h"], k[k.proc == p]["cpu_h"]
        print(f"{p}: healthy n={len(hh)} median={hh.median():.1f} "
              f"[{hh.min():.1f}-{hh.max():.1f}] CPU-h; killed n={len(kk)} "
              f"median={kk.median():.1f} [{kk.min():.1f}-{kk.max():.1f}] CPU-h")
    # SpliceAI: reference = 4-slot chunks on the five nodes with no slow chunk seen.
    ref = h[(h.proc == "SPLICEAI") & (h.slots == 4) & ~h.node.isin(["c0051-cn1", "c0053-cn1"])]
    base = ref["cpu_h"].median()
    print(f"\n== SpliceAI CPU-h relative to reference ({len(ref)} chunks, 5 nodes, median {base:.1f} CPU-h, "
          f"range {ref.cpu_h.min():.1f}-{ref.cpu_h.max():.1f}; cores_busy {ref.cores_busy.min():.2f}-{ref.cores_busy.max():.2f}) ==")
    s = h[h.proc == "SPLICEAI"].assign(ratio=lambda d: d.cpu_h / base)
    print(s.groupby("node")["ratio"].agg(["count", "min", "median", "max"]))
    sk = k[k.proc == "SPLICEAI"].assign(ratio=lambda d: d.cpu_h / base)
    print("killed SpliceAI (unfinished, so lower bounds):", sk.groupby("wall_cls")["ratio"].agg(["count", "min", "max"]).to_dict("index"))
    print("cores_busy, SpliceAI, all healthy: %.2f-%.2f; killed: %.2f-%.2f" % (
        s.cores_busy.min(), s.cores_busy.max(), sk.cores_busy.min(), sk.cores_busy.max()))
    print("\n== healthy: wall vs cpu correlation (does CPU grow with wall time?) ==")
    for p, g in h.groupby("proc"):
        print(p, "spearman(cpu_h, wall_h) = %.2f" % g["cpu_h"].corr(g["wall_h"], method="spearman"),
              "| spearman(cores_busy, wall_h) = %.2f" % g["cores_busy"].corr(g["wall_h"], method="spearman"))
    print("\n== healthy: slowest 5 per process ==")
    print(h.sort_values("wall_h", ascending=False).groupby("proc").head(5)[
        ["step", "slots", "cpu_h", "wall_h", "cores_busy", "cpu_eff"]])


if __name__ == "__main__":
    main()
