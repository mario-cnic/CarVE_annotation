# Slow and hung Nextflow tasks: S223 (2026-09-25..29) and `panel7_v45_retest2` (2026-10-01)

**2026-10-02.** Question: is `MISC-16` (six trivial tasks hung on the `panel7` run) related to the performance problems of the S223 run?
Method: `src/tools/analyze_nextflow_tasks.py` on `.nextflow.log.5` (S223) and `.nextflow.log` (panel7) plus the task work directories (read-only). Outputs of both runs: `RUNS/transcript_mapping/{s223,panel7_retest2}_task_analysis.txt` (not tracked).

## FACTS
**The same symptom, two different causes.** Nextflow reports "terminated for an unknown reason" when SGE kills a task at its `h_rt` limit (no `.exitcode` is written). S223: 907 completed-task lines, **26 without exit status = 540 task-hours** (18 chunks killed at 24 h, 2 at 48 h, 6 trivial tasks at 2 h). panel7: 6 of 6, all trivial.

**1. Long chunks (SpliceAI, Pangolin): ~10x slower on specific nodes, not stuck.**
| Node | tasks | no exit status | SpliceAI / Pangolin successes |
|---|---|---|---|
| c0053-cn1 | 29 | **15 (52%)** | one SpliceAI in 47.9 h, one Pangolin in 34.9 h |
| c0051-cn1 | 95 | 3 | SpliceAI median 8.6 h, max 46.1 h |
| c0051-cn4 | 196 | 3 (two of them Pangolin at 24 h) | normal |
| c0051-cn2, -cn3, c0053-cn2, -cn3 | 136-166 each | 0-2 (trivial tasks only) | SpliceAI median 6.1-6.4 h, Pangolin 14.8-21.2 h, max 23.7 h |
- Healthy baseline per 8,000-variant chunk: SpliceAI ~6.2 h, Pangolin ~15-21 h, SPiP ~0.1 h. Chunk size is constant, so variance on healthy nodes is small.
- Evidence they were slow, not hung: a killed SpliceAI chunk was still printing `3 s/step` (healthy: `263 ms/step`); a killed Pangolin chunk had reached line 3,318 of 8,000 after 24 h (healthy chunk: 8,000 lines in 1.79 h). Both left a 1 MiB (buffered) partial output and no `.exitcode`.
- Two SpliceAI chunks were retried (attempt 2, 48 h limit) and landed on c0053-cn1 again: killed again at 48 h.
- This supports, and sharpens, the open hypothesis in `TODO.md` (slowness tracks the node, not variant count or gene density). Cause of the slowdown (CPU contention from co-located 8-12 slot jobs, throttling, `OMP_NUM_THREADS=task.cpus` oversubscription, a node problem) is NOT established.

**2. Trivial tasks (`CONCAT_*`, `MERGE_ANNOTATIONS`): the FIRST attempt hung, the retry took ~1 s - in both runs.**
- All six `process_single` tasks, in both S223 and panel7 (12 of 12): attempt 1 ran until `h_rt=02:00:00` with zero stdout/stderr; in panel7 the script never created its first output file (`tmp_spip.vcf.gz`); nodes c0051-cn2/-cn3/-cn4, c0053-cn3/-cn4 (not node-specific).
- The retry ran the byte-identical `.command.sh` and finished in 0.8 s, with different SGE resources: `h_rt` 02:00 -> 04:00, `h_rss`/`mem_free` 4096M -> 8192M, `h_vmem` 14G -> 18G (`process_single`: memory `4.GB * attempt`, time `2.h * attempt`).
- Not seen for other labels in the same runs (`CHUNK_VCF`/`process_low`, `process_medium` ran normally). So it looks tied to attempt 1 of `process_single`, but the mechanism is unknown.
- S223 cost of this: 12 task-hours, but each hung attempt also delayed every downstream step by 2 h.

## INFERENCE / not done
- The two problems are related only in how SGE/Nextflow report them. They need separate fixes.
- Not analyzed: the 6 `exit 1` tasks of S223 (2x `VCF_TO_TABLE`, 4x `TAG_TRANSCRIPT_PRIORITY`) - different failure type, S223 was later completed.
- Options (none applied; configuration changes need a decision): (a) SGE host exclusion of c0053-cn1 for `process_long` tasks (and watching c0051-cn1) or ask the cluster admins to check the node; (b) shorter `h_rt` plus smaller chunks so a slow chunk is retried early instead of after 24 h; (c) for `process_single`, give attempt 1 the resources of the successful retry and a short `h_rt` (15-30 min) so a hung attempt is retried quickly; (d) test (c)'s premise first: a `qrsh -l h_rss=4096M,mem_free=4096M,h_vmem=14G` session running a COPY of the hung task's inputs and `timeout 120 bash .command.sh`.

## Update 2026-10-02 (later): MISC-16 mitigation applied; cause still open
**More facts gathered before choosing the fix**
- Legitimate runtimes of the trivial tasks at S223 scale (successful attempts): `CONCAT_*` 0.3-11 min, `MERGE_ANNOTATIONS` 16.1 min, `CHUNK_VCF` <1 min.
- The retry ran on a DIFFERENT node than the hung first attempt in 11 of 12 cases, and once (S223 `CONCAT_BRANCHPOINT`, c0053-cn3) on the SAME node, where it finished in 1.5 min. So it is not a bad node.
- `CHECK_INPUT_ASSEMBLY` also carries the `process_single` label and its first attempt ran normally (it runs first, when nothing else is on the cluster). So the label alone is not the trigger. The hung tasks all started while heavy SpliceAI/Pangolin jobs were running on the cluster.
- Only two runs ever reached these tasks in the available logs (S223, `panel7_v45_retest2`); both: 6 of 6 first attempts hung.

**Hypotheses (not tested):** (1) memory/swap pressure on nodes shared with 12-slot SpliceAI/Pangolin jobs: `mem_free=4096M` is only checked at dispatch, and the retry's `mem_free=8192M` selects less loaded nodes (this would also fit the ~10x slow chunks of `MISC-17`; `nextflow.config` itself notes that a missing h_rss reservation "causes silent disk-thrashing"); (2) a filesystem stall on inputs written a moment earlier by other nodes (does not explain why `VCF_TO_TABLE`, with equally fresh input, ran fine); (3) a scheduler effect of the resource request.

**Mitigation applied (`nextflow.config`, `process_single`, not a verified fix):** memory `4.GB * attempt` -> `8.GB * attempt` (attempt 1 now gets what the successful retries had) and time `2.h * attempt` -> `1.h * attempt` (a hung attempt is retried after 1 h instead of 2 h; 3.7x above the longest legitimate runtime seen). Checked: config parses; a local `-stub-run` of the whole-VCF workflow completes (wiring only; a stub does not run the scripts, so it cannot show whether the hang is gone).

**How to tell whether it worked / find the cause**
- Next full run: look at the `CONCAT_*` / `MERGE` tasks: `python3 src/tools/analyze_nextflow_tasks.py --log .nextflow.log` should list no `process_single` task without exit status.
- Decisive diagnostic (cluster login node; I cannot run `qacct`): `bash resources/pull_qacct_snapshot.sh RUNS/transcript_mapping/misc16_17_job_ids.tsv RUNS/transcript_mapping/misc16_17_qacct.tsv`. The list holds 44 jobs: the 12 hung first attempts, their 12 successful retries and the 20 SpliceAI/Pangolin chunks killed at `h_rt`. Read it as: `cpu` close to 0 with `ru_wallclock` 7200 s = blocked (I/O or lock wait, supports hypotheses 2/3); `cpu` close to wallclock = busy; `maxvmem` near the node's memory = pressure (hypothesis 1).

## Result of the `qacct` diagnostic (2026-10-02; `RUNS/transcript_mapping/misc16_17_qacct.tsv`, 44 jobs)
**MISC-16 - the hung first attempts were busy, not blocked, and tiny.** All 12 (S223 and panel7): exit 137 (SGE kill at `h_rt`), CPU time 7,127-7,183 s over 7,201 s wall = **99-100% of one CPU for the whole 2 h**, `maxvmem` only **0.022-0.027 GB** (about 25 MB). No output file was ever created. The 12 successful retries: 1 s (panel7) and 21-965 s (S223: `CONCAT_VEP` 659 s, `MERGE` 965 s), `maxvmem` up to 1.0 GB.
- Consequences: hypothesis (1) memory/swap pressure is REFUTED (25 MB footprint); a blocked I/O or lock wait is REFUTED (the CPU was fully used). What remains: a CPU-bound spin at the very start of the command (before its first output is created), only on attempt 1, in a process using ~25 MB. Not identified which process spins (the wrapper shell or the first `bcftools`/`tabix` call).
- Consequence for the mitigation: the memory increase (4 -> 8 GB) is NOT supported by this evidence (a 25 MB process is not memory-limited); the retry differed in `h_rt`, `h_rss`, `mem_free` and `h_vmem` together, so which of them matters is unknown. The 1 h time limit still limits the cost (1 h instead of 2 h per hung attempt). Both stay in place until a run shows whether the hang is gone.

**MISC-17 - the killed long chunks were busy and small as well.** 20 jobs (S223, 4 slots each), exit 137 at 86,401 s (24 h; two at 48 h): CPU time ~227,000 s each (63 CPU-hours), i.e. 26-66% of the slot time (median 57% Pangolin, 66% SpliceAI = about 2.3-2.6 of 4 cores busy), `maxvmem` only 1.1-1.7 GB. So memory pressure is refuted here too. They kept consuming CPU for 24 h without finishing a chunk that takes ~6 h (SpliceAI) on a healthy node. What this suggests (not proven): the CPU time was spent without progress (spinning threads, contention or oversubscription) - to be compared with healthy chunks.
- Next diagnostic (cluster login node): `bash resources/pull_qacct_snapshot.sh RUNS/transcript_mapping/misc17_healthy_job_ids.tsv RUNS/transcript_mapping/misc17_healthy_qacct.tsv` (74 healthy chunk jobs from all 7 nodes, labelled with process, node and hours). Compare CPU time per chunk and cores busy between healthy chunks and the killed ones: if healthy chunks need far less CPU time for the same work, the killed ones were wasting CPU.

## MISC-17: healthy vs killed chunks, per node (2026-10-02T09:39+02:00)
Inputs: `RUNS/transcript_mapping/misc17_healthy_qacct.tsv` (74 jobs; 3 tail chunks of 96-115 s excluded -> 71: 37 SpliceAI, 34 Pangolin) and `misc16_17_qacct.tsv` (20 killed chunks). Node of every job taken from the Nextflow log (`.nextflow.log.5`: job id -> work dir -> 2nd line of `.command.log`; 71/71 healthy labels agree with it); qacct itself has no hostname in these TSVs. Script: `src/tools/analyze_misc17_qacct.py` (read-only); full output `RUNS/transcript_mapping/misc17_qacct_comparison.txt` (gitignored). Units: `cpu` = SGE CPU seconds summed over slots; cores busy = cpu / wallclock.

### FACTS - SpliceAI (the clean signal)
- **Cores busy is the same everywhere: 2.56-2.72 for all 37 healthy chunks and 2.58-2.65 for the 9 killed ones**, on every node, 4 or 8 slots. Threads are pinned (TF 2+2, `src/python/annotate_spliceai.py:17-18`), so a chunk is ~2.6 cores whatever the slot request; the 8-slot retries bought nothing.
- CPU per chunk is what varies. Reference = 30 four-slot chunks on c0051-cn2/-cn3/-cn4, c0053-cn2/-cn3: **median 17.0 CPU-h (10.9-19.3), wall median 6.4 h**. Relative CPU-h per node: five reference nodes 0.64-1.14x; **c0051-cn1 1.07-2.18x (median 1.5x; three chunks at 31-37 CPU-h, wall 11.9-14.0 h)**; **c0053-cn1 7.2x**. Across all 37 chunks wall time tracks CPU time (Spearman 0.98) while cores busy does not (-0.15): a slow chunk is a chunk that burned more CPU-seconds for the same 8,000 variants.
- c0053-cn1: the 7 killed 4-slot chunks stopped at 63 CPU-h (>= 3.7x the reference); the 3 attempt-2 chunks (8 slots, 48 h limit) needed 122-124 CPU-h (7.2-7.3x): one finished at 47.9 h, two were killed at 48.0 h with 123.9-124.1 CPU-h, i.e. within ~1.5% of the finishing CPU of the one that made it.

### FACTS - c0053-cn1 as a node
- **It never completed a 4-slot SpliceAI/Pangolin chunk: 13 of 13 first attempts (7 SpliceAI, 6 Pangolin) were killed at 24 h**, submitted Sep-25 10:38-12:24 (+ attempt-2 chunks on Sep-26 10:44). Chunks submitted to other nodes in the same minutes ran normally (e.g. c0051-cn4 SpliceAI Sep-25 10:45, 6.7 h). So not a time-window effect within S223; the slowness lasted >= 26 h.
- Its only finished chunks are the two slowest of all 71 healthy ones: SpliceAI 47.9 h / 122 CPU-h; Pangolin 34.9 h / 112 CPU-h at 3.2 cores busy (next highest Pangolin anywhere: 45 CPU-h, 23.7 h). Both are 8-slot attempt-2 chunks, n = 1 each.
- Killed Pangolin on c0053-cn1: 55-69 CPU-h in 24 h at 2.0-2.9 cores; Pangolin elsewhere is mostly wait-dominated (median 0.9 cores, see below). Same node, both tools CPU-saturated and slow.
- Loss accounting: of the 528 killed task-hours, **408 (77%) are c0053-cn1** (13 x 24 h + 2 x 48 h), 120 (23%) are Pangolin on c0051-cn1 (3) and c0051-cn4 (2).

### FACTS - Pangolin (no clean baseline)
- Healthy Pangolin chunks are not one population: 3 finished at ~3.9 cores in 0.5-5.1 h (2-20 CPU-h), most of the rest took 8-23.7 h at 0.03-2.2 cores (0.35-45 CPU-h; one 8-slot chunk 5.9 h at 2.4 cores). CPU time is not a proxy for work here; no per-node CPU-ratio is computed. Why Pangolin chunks differ so much on healthy nodes is not explained by these data (chunk content not examined).
- **The 24 h `h_rt` sits inside the healthy wall-time range.** Healthy wall per node (median / max): c0051-cn4 22.1 / 23.2 h, c0053-cn2 20.5 / 23.7 h, c0051-cn1 19.1 / 22.5 h, others 14-16 h. All 4 finished 4-slot Pangolin chunks on c0051-cn4 took 21.8-23.3 h; 2 more on that node were killed at 24 h with 25 CPU-h at 1.04 cores - the same CPU profile as the finished ones (16-21 CPU-h), i.e. the tail beyond the limit, not the c0053-cn1 pattern.
- c0051-cn1 Pangolin: 3 killed (43, 93, 46 CPU-h in 24 h) vs 6 finished (12.8-22.5 h): mixed, not classified.

### What the data supports / does not support
- **Supports (FACT):** (1) c0053-cn1 differs from healthy nodes even when its chunks finish: same ~2.6 cores busy but 7x the CPU-seconds per SpliceAI chunk, 13/13 first attempts killed. (2) c0051-cn1 shows the same SpliceAI signature at a milder, graded level (1.1-2.2x). (3) Not memory (<= 1.7 GB), not an allocation of too few cores (cores busy unchanged), not too many threads inside the job for SpliceAI (pinned 2+2). (4) 77% of the lost task-hours trace to one node.
- **Supports as INFERENCE only:** each CPU-second does less work on c0053-cn1 (slower effective cores: frequency/power throttling, a different CPU, SMT or memory-bandwidth sharing, hypervisor effects), or the pinned threads burn CPU without progress (spinning on a stalled resource). Constant cores busy rules out starvation by time-slicing (that would LOWER cores busy), but does not separate these two.
- **Does not support / unknown:** the cause; whether c0053-cn1 is still slow (S223 ran 2026-09-25..29; one run only; panel7 had no long chunks); whether Pangolin's 24 h kills elsewhere are a node effect (they look like the wall-time tail); thread settings as a fix (SpliceAI threads already pinned and not varying; Pangolin sets `OMP_NUM_THREADS=${task.cpus}`, `modules/local/pangolin.nf:26`, so its 8-slot chunks (1.7-3.2 cores) are not comparable to 4-slot ones). n = 1 successful chunk per tool on c0053-cn1, both attempt-2.

### PROPOSAL (for review; `nextflow.config` NOT edited)
1. **Exclude c0053-cn1 for `process_long`** (SPLICEAI_ANNOTATE, PANGOLIN_ANNOTATE: the only two users of the label). Works whatever the cause, reversible, would have removed 408 of 528 lost task-hours and the retry-lands-on-the-same-node problem. Implementation caveat: `process.clusterOptions` is one closure in the `standard` profile (`nextflow.config:169-171`: `-P`, `-A`, `h_vmem`, `TMPDIR`); a `clusterOptions` inside `withLabel: process_long` would REPLACE it and silently drop those flags, so add a conditional inside the existing closure instead. SGE syntax to verify first on the login node: `qsub -w v -l h='!c0053-cn1' -b y /bin/true` (verification only, nothing runs). After approval I check the config renders and a local `-stub-run` completes.
2. **Ask the admins about c0053-cn1** (and c0051-cn1) with evidence: run `resources/misc17_node_benchmark.sh` pinned to c0053-cn1, c0051-cn1 and one healthy node (commands in `TODO.md`). If the fixed pure-Python loop is slower on c0053-cn1, it is the hardware/hypervisor; if the loop is normal but the real SpliceAI chunk is slow, look at memory/threads/IO there. The same run also shows whether the node is still slow today.
3. **Pangolin 24 h limit as a separate item**, not fixed by (1): healthy wall up to 23.7 h. Option: `withName: PANGOLIN_ANNOTATE { time = { 36.h * task.attempt } }` (36 h = 1.5x the longest healthy chunk; a guess, not derived). A SHORTER `h_rt` at the current chunk size would kill healthy Pangolin chunks; shorter `h_rt` + smaller chunks only makes sense after (1) and only if Pangolin time scales with variant count (not shown).
4. **Not proposed:** thread-setting changes (no evidence for SpliceAI), shortening SpliceAI `h_rt` (c0051-cn1's healthy-but-slow chunks reach 14 h), changing slots (SpliceAI uses ~2.6 of 4; the 8-slot retry wastes slots but is not the cause).
5. **Watch c0051-cn1** (SpliceAI 1.1-2.2x, 3 Pangolin kills) after the next run: `analyze_misc17_qacct.py` on the new qacct pull.

### Commands for the login node (I cannot run SGE; nothing below has been run)
```
# (a) syntax check of the proposed exclusion - verification only, nothing is submitted
qsub -w v -l h='!c0053-cn1' -b y /bin/true

# (b) qacct detail for 5 jobs: stime/utime split and involuntary context switches (spinning vs slow cores)
for j in 5012707 5011431 5011808 5011581 5012021; do echo "== $j"; \
  qacct -j $j | egrep '^(hostname|slots|ru_wallclock|ru_utime|ru_stime|cpu|ru_nvcsw|ru_nivcsw|ru_minflt|ru_majflt|io|iow|maxrss|start_time|end_time)'; done
#   5012707 c0053-cn1 SpliceAI ok 47.9 h | 5011431 c0053-cn1 SpliceAI killed 24 h | 5011808 c0051-cn1 12.9 h
#   5011581 c0051-cn1 7.2 h | 5012021 c0053-cn3 6.9 h (reference)

# (c) node benchmark: same SpliceAI chunk (first 100 variants of S223 chunk_0009), same flags as process_long attempt 1
R=/data_lab_PGP/pipelines/annotation_pipeline_new; cd $R; mkdir -p RUNS/transcript_mapping/misc17_bench
CH=$R/work/4a/2bb9ec1a24d841fb820c12ec6ef3c0/chunk_0009.vcf.gz
FA=/references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta
for n in c0053-cn1 c0051-cn1 c0051-cn2; do
  qsub -N misc17_bench_$n -S /bin/bash -P BIGN -A PGP -pe smp 4 -l h=$n -l h_rt=02:00:00 \
       -l h_rss=12G -l mem_free=12G -l h_vmem=22G -j y -o RUNS/transcript_mapping/misc17_bench/$n.out \
       resources/misc17_node_benchmark.sh $CH 100 $FA
done
# afterwards: grep BENCH RUNS/transcript_mapping/misc17_bench/*.out
# (if `-l h=<node>` needs the FQDN on this cluster, use it; healthy SpliceAI: ~6.4 h / 8,000 variants = ~5 min for 100)
```
Reading (c): the fixed loops (`python_loop`, `numpy_matmul_1thread`, `numpy_memcopy`) report CPU seconds for identical work. Several times higher on c0053-cn1 than on c0051-cn2 -> slow hardware/hypervisor (evidence for the admins). Loops equal but `spliceai_real` slower -> something specific to the TF workload there (memory bandwidth, IO, threads). `BENCH rusage` and (b): high `ru_stime` / `ru_nivcsw` on the slow node points to spinning/contention rather than slower cores.

## Update 2026-10-02: host exclusion applied (user decision)
`nextflow.config`: new param `long_task_excluded_hosts = ['c0053-cn1']`; the `standard` profile's `clusterOptions` closure appends ` -l h=!c0053-cn1` when `task.process` ends in `SPLICEAI_ANNOTATE` or `PANGOLIN_ANNOTATE` (the only `process_long` users). Checks: SGE verification `qsub -w v -P BIGN -A PGP -l h='!c0053-cn1' -b y /bin/true` -> "found suitable queue(s)" (run by the user; the first attempt failed only because `-P` came after the command); the closure, extracted verbatim from the config, evaluated in a throwaway pipeline: SpliceAI/Pangolin get `... -l h_vmem=22G -l h=!c0053-cn1`, `CONCAT` keeps `... -l h_vmem=18G`, `--long_task_excluded_hosts false` removes it; `nextflow config -profile standard` parses; local `-stub-run` (`local_dev`, log/work/outdir in the job tmp dir) completes. Not verified on the cluster (the local profile never evaluates this closure for SGE). Next full run: check `#$ -l h=` in a task's `.command.run` and `python3 src/tools/analyze_nextflow_tasks.py --log .nextflow.log` for tasks on c0053-cn1. Still open: Pangolin `h_rt`, admin ticket, c0051-cn1.

## Benchmark and qacct detail results (2026-10-02, ~11:00+02:00): c0053-cn1 is slow hardware
Jobs 5014281 (c0053-cn1), 5014282 (c0051-cn1), 5014283 (c0051-cn2), `RUNS/transcript_mapping/misc17_bench/*.out`; qacct detail of 5 S223 jobs (pasted by the user).

### FACTS
- **Same fixed single-thread work, CPU seconds (= wall seconds in every case):** pure-Python loop 17.5 s on c0053-cn1 vs 1.9 s on c0051-cn1 and c0051-cn2 (**9.2x**); numpy matmul 3.3 vs 0.36 s (9.2x); 1.6 GB memcopy x5 5.6 vs 1.24 s (4.5x). c0051-cn1 equals c0051-cn2 today (one 2.26 s blip in the loop). CPU = wall means the extra time is real CPU, not waiting, contention or spinning.
- **Node facts** (same CPU model on both: AMD EPYC 9654, 384 logical CPUs, min/max 1500/3709 MHz, governor `performance` on both): c0053-cn1 had **18 logical CPUs at ~400 MHz (below the 1500 MHz minimum)** and the rest at 2400 MHz, loadavg 3.3/3.1/3.1 with no visible CPU-heavy process; c0051-cn2 idle (loadavg 0.0), the active cores at ~3,700 MHz.
- **qacct, finished S223 chunks (SpliceAI):** system-time share of CPU is 17-20% on every node, including c0053-cn1 (16.8%) -> no excess kernel/spin time. Involuntary context switches per CPU-second: 8.6 on c0053-cn1 vs 21.7 on the reference node -> no contention. Total voluntary context switches per chunk are about the same on all (c0053-cn1 2.9e9; c0051-cn1 slow-mode 2.5e9, normal 3.0e9; reference 3.1e9) while CPU time differs 1.7-6.7x. Killed jobs' `ru_utime`/`ru_stime` (0.6 s / 1.4 s for 5011431) are just the wrapper, unusable.
- c0051-cn1 benchmarks normal today, so its earlier graded slowness (SpliceAI 1.1-2.2x) is not reproduced.

### INFERENCE
- c0053-cn1's cores deliver ~1/9 of the work per CPU-second; the 9.2x equals 3,700 MHz / 400 MHz (9.25), and 400 MHz is below the 1500 MHz P-state floor with the `performance` governor, so the limit is probably firmware/hardware (power cap, thermal throttle, stuck low-power state), not the OS. The match could be a coincidence: the frequency was read once, not tied to the core that ran the benchmark. Fits the S223 data: 7.2x CPU for the same SpliceAI work, same cores busy, no extra sys/spin time. Pangolin on that node (killed at 2.0-2.9 cores) fits too.
- The node is still slow 5 days after S223 -> not transient.
- Unexplained: loadavg ~3 with nothing visible (possible tasks in uninterruptible wait); why c0051-cn1 was slower in S223.

### Not done / bug in my script
- The real-SpliceAI part did not run on any node: `/usr/bin/time` does not exist on the compute nodes (script bug, fixed: bash `time`). Not needed for the conclusion.
- Conclusion for the fix: the c0053-cn1 exclusion (already applied) addresses the cause; keep it until the admins fix the node. Message for the admins is in `TODO.md`.

## Why is S223 Pangolin so variable? (2026-10-02) - it ran with the MOUSE database
Script `src/tools/analyze_pangolin_chunk_runtime.py` (read-only: work dirs of the 48 Pangolin jobs in the two qacct TSVs), output `RUNS/transcript_mapping/misc17_pangolin_runtime.txt` (gitignored).
### FACTS
- **All 48 tasks ran `pangolin_grch38.db`, the mouse annotation db of `MISC-12`** (chr1-chr19 only; `.command.sh`). `.command.err` has one `[Line N] WARNING, skipping variant` per unscored variant ("not contained in a gene body", or "format not supported"). So a human variant is "scored" only when it coincides with a mouse gene body; S223's Pangolin times describe a wrong workload.
- The instant chunks (96 s, 3 of 37) contain only chr20-chr22: those contigs are absent from the db, so every variant is skipped without a lookup. Two other chunks also skipped all 8,000 variants (0 scored) but took 8.4 h and 12.6 h with 0.35-0.42 CPU-h (0.03-0.04 cores): their variants are on chr3 / chr8-9, which ARE in the db, so each lookup found nothing but took ~4-6 s of waiting, not CPU. Chunks with CPU-bound behaviour (3.9 cores, 0.5-5 h) had many variants "scored".
- Among the 32 finished chunks with >= 5,500 variants on chr1-19 (wall 1.8-35 h) neither the number of chr1-19 variants nor of scored variants explains wall time (Spearman 0.37 and 0.47); only CPU time tracks it (0.78, trivially). Node does not explain it either (same node hosted both instant and 12 h chunks).
### INFERENCE
- The 14-24 h "healthy Pangolin baseline" and the 24 h limit question come from a run with the wrong db; the real per-chunk time with `gencode.v45.ensembl_canonical.grch38.db` (human, scores every in-gene variant) is UNKNOWN. The wait-dominated chunks suggest slow per-variant db lookups (sqlite file on network storage); not tested.
- So the Pangolin `h_rt` (e.g. 36 h) should not be decided now: wait for the S223 Pangolin re-run with the correct db (already pending, `MISC-12`) and measure with `analyze_misc17_qacct.py` / `analyze_pangolin_chunk_runtime.py`. The `c0053-cn1` exclusion is unaffected (its slowness shows on SpliceAI and in the node benchmark).

## Update 2026-10-02: Pangolin `h_rt` raised (user decision)
`nextflow.config`: `withName: 'PANGOLIN_ANNOTATE' { time = { 36.h * task.attempt } }` (was 24 h from `process_long`; SpliceAI stays 24 h). **36 h is a guess** (1.5x the longest healthy S223 chunk, 23.7 h), not derived. **The run-time variance of Pangolin chunks is UNEXPLAINED** (previous section: mouse db in S223; instant chr20-22 chunks explained, the 1.8-35 h spread of the rest not), so this raises the ceiling without understanding what sets it. Cost: a genuinely stuck Pangolin chunk now wastes up to 36 h before its retry (72 h on attempt 2). Checked: `nextflow config` shows the override; a local `-stub-run` trace (`trace.fields = name,attempt,time`) lists `PANGOLIN_ANNOTATE` 1d 12h, `SPLICEAI_ANNOTATE` 1d, others unchanged. To do: after the S223 Pangolin re-run with the human db, re-measure (`analyze_misc17_qacct.py`, `analyze_pangolin_chunk_runtime.py`), then tighten or raise the value.
