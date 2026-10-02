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
