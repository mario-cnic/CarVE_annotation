#!/usr/bin/env python3
"""MISC-17/MISC-12: what explains the run-time variability of S223 Pangolin chunks? (read-only)

For every PANGOLIN_ANNOTATE job in the two qacct TSVs: locate its work dir from the Nextflow log
(job id -> workDir), read raw_input.vcf (chunk input) and .command.err (Pangolin prints one
"[Line N] WARNING, skipping variant ..." per variant it does not score), and relate the share of
variants Pangolin actually scored to wall time / CPU. Also reports which Pangolin db each task used
(from .command.sh). Nothing is written except stdout.
"""
import argparse
import collections
import os
import re

import pandas as pd

from analyze_misc17_qacct import load, job_nodes  # same directory


def work_dirs(nxf_log):
    out = {}
    pat = re.compile(r"submitted process .* > jobId: (\d+); workDir: (\S+)")
    with open(nxf_log, errors="replace") as fh:
        for line in fh:
            m = pat.search(line)
            if m:
                out[int(m.group(1))] = m.group(2)
    return out


def inspect(wd):
    r = {"n_chr1_19": None, "n_input": None, "n_warn": None, "last_line": None, "db": None, "warn_kinds": None}
    try:
        with open(wd + "/raw_input.vcf", errors="replace") as fh:
            ch = [l.split("\t", 1)[0] for l in fh if not l.startswith("#")]
        r["n_input"] = len(ch)
        # the S223 db is MOUSE (MISC-12): contigs chr1-chr19 only, so human chr20-22/X/Y/M cannot be looked up
        r["n_chr1_19"] = sum(1 for c in ch if re.fullmatch(r"chr(\d+)", c) and 1 <= int(c[3:]) <= 19)
    except OSError:
        pass
    try:
        kinds, last = collections.Counter(), 0
        with open(wd + "/.command.err", errors="replace") as fh:
            for l in fh:
                m = re.match(r"\[Line (\d+)\] (.*)", l)
                if m:
                    last = max(last, int(m.group(1)))
                    kinds[re.sub(r"\s+", " ", m.group(2))[:60]] += 1
        r["n_warn"], r["last_line"], r["warn_kinds"] = sum(kinds.values()), last, dict(kinds)
    except OSError:
        pass
    try:
        m = re.search(r"(\S+\.db)\b", open(wd + "/.command.sh").read())
        r["db"] = os.path.basename(m.group(1)) if m else None
    except OSError:
        pass
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--healthy", required=True)
    ap.add_argument("--killed", required=True)
    ap.add_argument("--nxf-log", required=True)
    a = ap.parse_args()
    wd, jn = work_dirs(a.nxf_log), job_nodes(a.nxf_log)
    h = load(a.healthy); h["outcome"] = "ok"
    k = load(a.killed); k = k[k.step.str.contains("KILLED")].copy(); k["outcome"] = "KILLED"
    d = pd.concat([h, k])
    d = d[d.proc == "PANGOLIN"].copy()
    d["node"] = d.jobnumber.map(lambda j: jn.get(int(j), ("?", ""))[0])
    info = d.jobnumber.map(lambda j: inspect(wd[int(j)]) if int(j) in wd else {})
    for c in ["n_chr1_19", "n_input", "n_warn", "last_line", "db", "warn_kinds"]:
        d[c] = info.map(lambda x: x.get(c))
    d["skip_frac"] = d.n_warn / d.last_line
    d["scored_lines"] = d.last_line - d.n_warn
    pd.set_option("display.width", 220, "display.max_rows", 200)
    print("db used by the tasks:", d.db.value_counts(dropna=False).to_dict())
    print("distinct warning kinds:", collections.Counter(k for w in d.warn_kinds.dropna() for k in w))
    cols = ["outcome", "node", "slots", "wall_h", "cpu_h", "cores_busy", "n_input", "last_line", "n_warn", "scored_lines", "n_chr1_19", "skip_frac"]
    print(d.sort_values("wall_h")[cols].to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    f = d[(d.wall_h > 0.1)]
    print("\nSpearman vs wall_h (chunks > 0.1 h):")
    for c in ["n_chr1_19", "scored_lines", "n_warn", "skip_frac", "n_input", "cpu_h", "cores_busy"]:
        print(f"  {c:13s} {f[c].astype(float).corr(f.wall_h, method='spearman'):.2f}")
    ok = f[f.outcome == "ok"]
    print("\nfinished chunks only (n=%d), Spearman wall_h vs n_chr1_19: %.2f; vs scored_lines: %.2f" % (
        len(ok), ok.wall_h.corr(ok.n_chr1_19.astype(float), method="spearman"), ok.wall_h.corr(ok.scored_lines.astype(float), method="spearman")))
    for lo, hi in [(0, 1), (1, 2500), (2500, 5500), (5500, 8001)]:
        g = ok[(ok.n_chr1_19 >= lo) & (ok.n_chr1_19 < hi)]
        if len(g): print(f"  n_chr1_19 in [{lo},{hi}): n={len(g)} wall_h median {g.wall_h.median():.1f} [{g.wall_h.min():.1f}-{g.wall_h.max():.1f}]")
    print("\nSpearman vs cpu_h:", {c: round(f[c].astype(float).corr(f.cpu_h, method='spearman'), 2) for c in ["scored_lines", "skip_frac"]})


if __name__ == "__main__":
    main()
