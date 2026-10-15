#!/usr/bin/env python3
"""What-if replay of the predictor fan-out of one run under Nextflow's queueSize cap (read-only).

Reads tasks.tsv from parse_nf_runs.py. Replays the per-chunk predictor tasks (VEP, Branchpointer,
SpliceAI, Pangolin, SPiP) in their observed submission order: a task starts when its chunking step
has finished, fewer than QUEUE tasks are in SGE and SLOTS has room for its cpus (observed SGE queue
wait was < 1 min in F2005 and S223 at 400 slots in use; total free slots above that are an
ASSUMPTION). Prints, per scenario, when each tool's last chunk ends and the estimated run wall-clock
= fan-out end + the run's measured post-fan-out time (concat, merge, table steps).
Needs a per-chunk tag on SpliceAI/Pangolin tasks (absent in S223's Pangolin tasks: refused).

Scenario knobs:
  --queue N            Nextflow executor.queueSize (100 = current default)
  --slots N            cluster slots available to the run (default: unlimited)
  --slow-hosts H,...   runtimes of successful tasks on these hosts are replaced by the process median on
                       the other hosts and killed attempts on them are dropped (= the host behaves like
                       a healthy one; to model excluding it, also lower --slots by its slot count)
  --restricted-scale F SpliceAI/Pangolin/SPiP chunks F times their size (0.5 = half-size chunks:
                       twice as many tasks, each half the runtime plus --task-overhead-min)
  --tail-ok-only       post-fan-out time with only the successful attempts (no OOM retries)
  --h-rt-hours H       SpliceAI/Pangolin time limit (see --help)
Killed attempts that are kept are replayed as observed: the slot is held until the kill, then the
retry runs for its observed runtime.
"""
import argparse, collections, csv, statistics as st

FANOUT = {"VEP_ANNOTATE": "BROAD:CHUNK_VCF", "BRANCHPOINT_ANNOTATE": "BROAD:CHUNK_VCF",
          "SPLICEAI_ANNOTATE": "RESTR:CHUNK_VCF", "PANGOLIN_ANNOTATE": "RESTR:CHUNK_VCF",
          "SPIP_ANNOTATE": "RESTR:CHUNK_VCF"}
RESTRICTED = {"SPLICEAI_ANNOTATE", "PANGOLIN_ANNOTATE", "SPIP_ANNOTATE"}
TAIL = ["CONCAT_PANGOLIN", "CONCAT_SPLICEAI", "MERGE_ANNOTATIONS", "VCF_TO_TABLE", "VCF_TO_TSV", "TSV_TO_TABLE",
        "TAG_TRANSCRIPT_PRIORITY", "ADD_WIDE_GENOTYPES"]


def f(x):
    return None if x in ("-", "") else float(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tasks", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--queue", type=int, default=100)
    ap.add_argument("--slots", type=int, default=10**9)
    ap.add_argument("--slow-hosts", default="")
    ap.add_argument("--restricted-scale", type=float, default=1.0)
    ap.add_argument("--task-overhead-min", type=float, default=2.0,
                    help="fixed per-task cost (model load, submit/poll) added to each split chunk")
    ap.add_argument("--tail-ok-only", action="store_true")
    ap.add_argument("--h-rt-hours", type=float, default=None,
                    help="time limit for SpliceAI/Pangolin attempts: a longer attempt holds its slot until the limit, "
                         "then the chunk re-runs at the process median (the retry keeps the slot: slightly optimistic)")
    a = ap.parse_args()

    rows = [r for r in csv.DictReader(open(a.tasks), delimiter="\t") if r["run"] == a.run]
    for r in rows:
        for k in ("submit", "start", "end", "realtime_s"):
            r[k] = f(r[k])
    t0 = min(r["submit"] for r in rows if r["submit"])
    ready = {p: max(r["end"] for r in rows if r["process"] == p) for p in set(FANOUT.values())}
    slow = set(filter(None, a.slow_hosts.split(",")))

    healthy = collections.defaultdict(list)
    for r in rows:
        if r["process"] in FANOUT and r["exit"] == "0" and r["realtime_s"] and r["host"] not in slow:
            healthy[r["process"]].append(r["realtime_s"])
    med = {p: st.median(v) for p, v in healthy.items()}

    # One job per chunk: its successful runtime, plus the held slot time of killed attempts that are kept.
    by_key = collections.defaultdict(list)
    for r in rows:
        if r["process"] in FANOUT:
            by_key[(r["process"], r["tag"] if r["process"] in ("SPLICEAI_ANNOTATE", "PANGOLIN_ANNOTATE")
                    else r["job_id"])].append(r)
    for p in ("SPLICEAI_ANNOTATE", "PANGOLIN_ANNOTATE"):
        n_ok = sum(1 for r in rows if r["process"] == p and r["exit"] == "0")
        n_key = sum(1 for k in by_key if k[0] == p)
        if n_ok and n_key < n_ok:
            raise SystemExit(f"{a.run}: {p} tasks have no per-chunk tag ({n_key} keys for {n_ok} chunks); cannot replay")
    jobs = []
    for (p, _), rs in by_key.items():
        rs.sort(key=lambda r: r["submit"])
        okr = [r for r in rs if r["exit"] == "0"]
        if not okr:
            continue
        okr = okr[-1]
        cpus = int(rs[0]["cpus_req"]) if rs[0]["cpus_req"] not in ("-", "") else 1
        run_s = okr["realtime_s"] if okr["host"] not in slow else med[p]
        killed = [(r["end"] - r["start"]) for r in rs if r["exit"] != "0" and r["start"] and r["host"] not in slow]
        if a.h_rt_hours and p in ("SPLICEAI_ANNOTATE", "PANGOLIN_ANNOTATE"):
            limit = a.h_rt_hours * 3600
            killed = [min(k, limit) for k in killed]
            if run_s > limit:
                killed.append(limit)
                run_s = med[p]
        held = sum(killed)
        if p in RESTRICTED and a.restricted_scale != 1.0:
            n = max(1, round(1 / a.restricted_scale))
            for i in range(n):
                jobs.append((rs[0]["submit"], p, run_s * a.restricted_scale + a.task_overhead_min * 60, held / n, cpus))
        else:
            jobs.append((rs[0]["submit"], p, run_s, held, cpus))
    jobs.sort()

    # List scheduling: FIFO in observed submit order, at most QUEUE tasks and SLOTS cpus in use at once.
    import heapq
    running, used, now, last_end = [], 0, 0.0, collections.defaultdict(float)
    for _, p, run_s, held, cpus in jobs:
        now = max(now, ready[FANOUT[p]])
        while running and (len(running) >= a.queue or used + cpus > a.slots):
            end, c = heapq.heappop(running)
            used -= c
            now = max(now, end)
        end = now + held + run_s
        heapq.heappush(running, (end, cpus))
        used += cpus
        last_end[p] = max(last_end[p], end)
    fan_end = max(last_end.values())

    # Post-fan-out: measured time from the last predictor chunk end to the run end.
    obs_fan_end = max(r["end"] for r in rows if r["process"] in FANOUT)
    run_end = max(r["end"] for r in rows if r["end"])
    tail_h = (run_end - obs_fan_end) / 3600
    if a.tail_ok_only:
        tail_h = 0.0
        for p in TAIL:
            okr = [r for r in rows if r["process"] == p and r["exit"] == "0"]
            if okr and p not in ("CONCAT_SPLICEAI",):  # concats run in parallel; count one
                tail_h += (okr[-1]["end"] - okr[-1]["start"] + 60) / 3600
    print("\t".join(["run", "queue", "slots", "slow_hosts", "restricted_scale", "overhead_min", "h_rt_h", "tail_ok_only"] +
                     [f"{p}_end_h" for p in FANOUT] + ["fanout_end_h", "tail_h", "est_wall_h", "observed_wall_h"]))
    print("\t".join([a.run, str(a.queue), str(a.slots) if a.slots < 10**9 else "inf", a.slow_hosts or "-", str(a.restricted_scale), str(a.task_overhead_min), str(a.h_rt_hours or "-"), str(a.tail_ok_only)] +
                    [f"{(last_end[p] - t0) / 3600:.1f}" for p in FANOUT] +
                    [f"{(fan_end - t0) / 3600:.1f}", f"{tail_h:.1f}", f"{(fan_end - t0) / 3600 + tail_h:.1f}",
                     f"{(run_end - t0) / 3600:.1f}"]))


if __name__ == "__main__":
    main()
