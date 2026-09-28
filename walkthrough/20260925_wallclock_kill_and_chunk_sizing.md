# Wall-clock kill on real WGS run + gene-restricted-tier chunk sizing (S223, MISC-11)

**2026-09-25**

## Context

While `S223` (real WGS run, `spontaneous_saha`) was in progress, user asked for a progress check.
Found 256 "terminated for an unknown reason -- likely by the external system" retries on
`PANGOLIN_ANNOTATE`/`SPLICEAI_ANNOTATE` specifically — after ~16h, only 1/47 Pangolin chunks and
0/47 SpliceAI chunks had actually succeeded, while VEP finished 255/255 cleanly under the same
shared `process_medium` label.

## FACT — root cause confirmed via `qacct`, not guessed

Initial hypothesis was OOM. Wrong — `qacct -j 5010792` on a killed Pangolin task showed:
`maxvmem 955.984MB` (out of a 34G budget), `failed 37: qmaster enforced h_rt, h_cpu, or h_vmem
limit`, `exit_status 137`, `ru_wallclock 28802s`, and `category` showing `h_rt=28800` — exactly
`4h * task.attempt(2)` from `process_medium`'s own `time` directive. The job ran its full allotted
8 hours and was hard-killed right at that ceiling, using almost no memory. A **wall-clock** problem,
not a memory one.

Per-chunk runtime is highly variable even at fixed chunk size: one 20000-variant Pangolin chunk
that *did* succeed finished in 12.5 minutes (`.command.begin` to `.exitcode`, real timestamps);
another, same chunk size, exceeded 8 hours and was still incomplete when killed.

## Fix

Two complementary changes, not one:

1. **New `process_long` label** (`nextflow.config`) for `PANGOLIN_ANNOTATE`/`SPLICEAI_ANNOTATE`
   only — `time = 24h * task.attempt`, same cpu/memory scaling as `process_medium` (memory was
   never the constraint, no reason to touch it). VEP/Branchpointer stay on `process_medium`,
   untouched — both already ran clean at real scale under it.
   **`maxRetries` overridden down to 2** (from the global 5) for this label specifically, per
   Mario 2026-09-25: "time shouldn't be a problem by itself" — the point of a generous 24h/48h
   ceiling is that a real chunk should already fit inside it, so escalating through 5 attempts
   (up to 240h/10 days total under the original 16h-base design) buys nothing once the base
   itself is generous; it would just delay surfacing a chunk that's stuck for some other reason.
2. **Gene-restricted-tier-specific chunk size**, separate from the broad tier's `chunk_size`.
   Required a real code change, not just a config value: `CHUNK_VCF` (`modules/local/chunk.nf`)
   previously read `params.chunk_size` directly inside its own script — both tiers were physically
   forced to share one value. Added `val(chunk_size)` as a real process input; each tier's own
   workflow (`broad_pass_whole.nf`, `gene_restricted_whole.nf`) now passes its own value at the
   call site.

**Value chosen for `gene_restricted_chunk_size`, with user input**: reducing from the shared
`20000` down to `3000` would maximize the worst-case-runtime reduction but pushes the
gene-restricted tier to ~313 chunks / ~944 total jobs (6.7x the original 47/146) — a real
scheduling-overhead cost on a shared cluster. User chose a moderate middle ground instead: `8000`
→ ~118 chunks / ~360 jobs (~2.5x more jobs than before, for a ~2.5x reduction in per-chunk size).
Full job-count math (broad tier unchanged at 255 chunks/513 jobs):

| chunk_size | gene-restricted chunks | gene-restricted jobs | full fresh-run total |
|---|---|---|---|
| 20000 (original) | 47 | 146 | 662 |
| 3000 (considered, rejected) | 313 | 944 | 1460 |
| **8000 (chosen)** | **118** | **360** | **876** |

Neither the `24h`/`maxRetries=2` time policy nor the `8000` chunk size is backed by real profiling data beyond
the two data points above (12min vs. >8h on the same old chunk size) — flagged explicitly as the
same class of untrustworthy-until-measured resource number `TODO.md` already warns about for
every other `h_vmem`/thread value in this pipeline. `MISC-11` in `BUG_TRACKER.md`.

## Verification

Static only — no `nextflow` invocation was run in this repo directory while `S223` was still
live (confirmed still running at the time of this fix; the user was explicitly asked before any
local verification, given the earlier `rm -rf work` incident this same session). Checked:
- `nextflow.config`, `modules/local/{chunk,pangolin,spliceai}.nf`,
  `workflows/{broad_pass_whole,gene_restricted_whole}.nf` all read back correctly after editing.
- Both `CHUNK_VCF` call sites now pass an explicit chunk-size argument; no other call sites exist
  anywhere in the codebase (`grep -rn "CHUNK_VCF("`).
- A real `-stub-run` verification is still owed once `S223` is confirmed killed.

## Unrelated finding caught while doing this: a real cross-repo ID collision

While flagging the separate Pangolin `heart_lv`/`heart_aa` misrepresentation
(`walkthrough/20260925_pangolin_tissue_of_origin.md`) as a `carve-platform` `PRED-#` finding,
`git status` surfaced uncommitted changes to this repo's own `TODO.md`/`BUG_TRACKER.md` from what
looks like separate, concurrent work (a non-coding-predictor gap analysis, unrelated to anything in
this session) that had already claimed `PRED-13` for an unrelated UTRAnnotator finding. Not touched
— left exactly as found. The Pangolin tissue finding was renumbered to `PRED-14` before it was
committed anywhere.

## Follow-up (2026-09-28): the `maxRetries=2` assumption was wrong, fixed live

Three of 118 `SPLICEAI_ANNOTATE` chunks on the real `S223` run ran past 46 hours on attempt 2
(confirmed via `qstat -j`: `h_rt=172800`, i.e. `24h * task.attempt(2)`), with `maxRetries=2`
meaning no further retry if killed at the 48h ceiling — despite each chunk actively producing real
inference output the whole time (`.command.log` updating live), not stuck. Removed the
`maxRetries: 2` override in `nextflow.config` (falls back to the global default of 5) for future
relaunches. For the three already-running jobs, `qalter` rescued them live without losing progress
— SGE requires the *entire* resource list restated, not just the field being changed:

```bash
qalter -l h_rss=24576M,h_rt=259200,h_vmem=34G,mem_free=24576M 5012705 5012707 5012712
```

(`h_rt=259200` = 72h). Confirmed applied via a second `qstat -j`, job still accumulating real CPU
time afterward. No progress lost, no relaunch needed.
