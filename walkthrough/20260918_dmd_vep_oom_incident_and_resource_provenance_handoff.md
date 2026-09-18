# Walkthrough — DMD VEP OOM Incident, Run-Provenance Gap, and Resource-Sizing Handoff

**Date**: 2026-09-16 to 2026-09-18
**Session focus**: literature review (SAI-10k-calc integration feasibility, SpliceAI `-D` window), a live production incident (DMD's VEP job stuck 19 days), and the provenance/reproducibility gaps it exposed.

---

## 1. Where things stand right now (read this first)

- **DMD's VEP job is still running** as of session end. The original attempt (job `4975105`, `-l h_vmem=15G -pe smp 4`) was SIGKILLed after 2 days at `maxvmem=56.885G`; it was resubmitted unchanged (job `4979674`) and ran 19 days pinned at the same ceiling with zero output. The user manually killed/resubmitted it outside `main.sh` with a larger allocation (`h_vmem=40G -pe smp 4` = 160GB total, confirmed `h_vmem` is per-slot on this cluster). **Its current status was not re-checked at session end.**
- **The user is going to run `resources/pull_qacct_snapshot.sh`** at some point after this session — do not re-derive job IDs or re-write this script; check whether `resources/run_more_genes_20260817_1143_qacct.tsv` already exists before doing anything else with real resource data.
- **Nothing in `TODO.md`'s Priority 0 section has been implemented yet.** All three items there (provenance manifest, SpliceAI `-D` window change, per-predictor resource scaling) are documented only, by explicit user instruction ("track it, don't implement it").

## 2. Next steps, in order

1. **Check DMD's VEP job.** `qstat -j 4979674` (or whatever job number the user's manual resubmit produced — get it from them if not obvious) to see if it finished, and if so:
   ```bash
   zcat RUNS/run_more_genes_20260817_1143/annotation/DMD.annVEP.vcf.gz | grep -vc "^#"
   ```
   must equal **533174** before trusting the file — `is_valid_file()` in `main.sh` only checks existence + non-zero size (`ORCH-6`), so a truncated/partial bgzip from an interrupted job would look "done" and get silently consumed by `merge` otherwise.

2. **Once DMD's VEP output is verified**, resume the rest of that gene's pipeline (merge → vcf2parsed → filter/plot/report/clinical_report) with:
   ```bash
   bash main.sh --raw-dir RUNS/predictors_050826/input/by_gene_more_genes --run-name run_more_genes_20260817_1143 --gene DMD
   ```
   The raw-dir was **not recorded anywhere** for this run and had to be recovered forensically (matched `RUNS/predictors_050826/input/by_gene_more_genes/` to this run by gene-count (154) and then confirmed by exact row-count match on `DMD.pq`, 533,174 rows). This is now the confirmed value — don't re-derive it.

3. **Once `resources/run_more_genes_20260817_1143_qacct.tsv` exists** (user-run, not yet done as of session end): analyze it to derive real, evidence-based `h_vmem`/thread tiers per predictor (VEP, SpliceAI, Pangolin, SPiP, Branchpointer), joined by gene against the real variant counts already available locally (`RUNS/predictors_050826/input/by_gene_more_genes/<GENE>.pq` row counts, or `_tmp/<GENE>.vcf.gz`). **Caveat carried over from the user**: the SpliceAI portion of this snapshot reflects the *current* `-D 10000` setting — if the `-D 4999` change (item below) has landed by the time you do this analysis, SpliceAI's tiers need a fresh snapshot; VEP/Pangolin/SPiP/Branchpointer tiers are unaffected by that change and remain valid from this snapshot.

4. **Then implement, in this order** (per `TODO.md` Priority 0 — do not skip ahead without checking with the user first, this session's pattern was strictly "propose → confirm → implement"):
   - Run provenance/reproducibility manifest (`run_manifest.json`/`.md` per `RUNS/<RUN_NAME>/`, capturing resolved CLI args, resolved `env.sh` binary paths, git commit hash, per-predictor tool/model/database versions).
   - Per-predictor resource scaling (`main.sh`'s per-gene qsub calls for VEP/SpliceAI/Pangolin/SPiP/Branchpointer — 2-3 tiers from real `qacct` data, not a fitted formula).
   - SpliceAI `-D 10000` → `4999` (this is a pipeline-wide behavior change affecting all future runs' SpliceAI scores, not just the SAI-10k-calc idea that originally prompted the literature review — treat it with matching care).

## 3. What was found this session (context for why the above matters)

### 3a. SAI-10k-calc (Canson et al. 2023, Bioinformatics btad179) — feasibility review
Reviewed the paper (`/home/mruizp/Downloads/btad179.pdf`, not yet copied into the repo) for a tool that turns raw SpliceAI delta scores into splicing-aberration-type/size/frame/AA-sequence predictions. Verdict: belongs in this annotation pipeline (deterministic, per-variant, same shape as SPiP/Pangolin/Branchpointer), not the `clinical_variant_prioritization` app. Logged as `TODO.md` Phase 5. Blocked on the SpliceAI `-D` window question below.

### 3b. SpliceAI `-D` window literature review
Confirmed via the actual SpliceAI source (`Illumina/SpliceAI` GitHub) that the official CLI hard-caps `-D` at `range(0, 5000)` — our pipeline's `-D 10000` (`src/python/annotate_spliceai.py`) only works because it calls `get_delta_scores()` directly, bypassing that cap. More importantly, `-D` changes the reported delta *score*, not just the delta *position* (`cov = 2*dist_var+1`, argmax taken over that window) — so existing `-D 10000` output cannot be reused for anything calibrated at `-D 4999` (Moles-Fernández 2021, SAI-10k-calc, ClinGen SVI Splicing Subgroup/Walker 2023 all use `4999`). Logged as a Priority 0 TODO item; **not implemented**.

### 3c. DMD VEP incident (`ORCH-10` in `BUG_TRACKER.md`)
User reported a job running since Aug 28 via a pasted `qstat -j`. Investigation found: `main.sh:442` gives every gene's VEP job the identical `-l h_vmem=15G -pe smp 4` regardless of variant count. DMD (533,174 variants — 2.3x the panel's next-largest gene, CACNA1C) was SIGKILLed once already at the same `maxvmem=56.885G` ceiling, then resubmitted unchanged and ran 19 days memory-thrashing with zero output. Every other gene in the run finished in hours under the same allocation, confirming this is a sizing bug specific to variant-count outliers, not a general hang. User later confirmed `h_vmem` is per-slot on this cluster and expanded the fix's scope to all predictor tools, not just VEP.

### 3d. Run-provenance gap (Priority 0 in `TODO.md`)
Triggered by trying to answer "what parameters did a run from a month ago actually use." Confirmed by direct inspection: no file anywhere records CLI args/thresholds, tool/model versions, git commit hash, or which auto-detected `env.sh` binary path was actually used for a given run. This gap became concrete twice in one session: (1) DMD's original `--raw-dir` had to be recovered forensically rather than looked up, and (2) there is no way to know in advance whether a future resource-sizing change actually took effect for a given run.

### 3e. `ORCH-11` — dormant timing/memory instrumentation
While looking for real historical memory data to calibrate resource tiers from (instead of guessing), found that `src/hpc/logger.sh`'s `/usr/bin/time -v` wrapper has **never once executed** in this pipeline's history: `main.sh` never exports `$LOG_DIR` (only a differently-named `$ERROR_LOG_DIR`), so the wrapper's activation check always fails. Confirmed by finding zero `.time` files and zero `pipeline_report.md` files anywhere under `RUNS/`. This is why the resource-tier research had to fall back to `qacct` (SGE's own independent accounting) instead of the pipeline's own (broken) instrumentation. One-line fix, not yet applied.

### 3f. `qacct` data-gathering detour
`qacct -o <user>` (broad owner query, no job ID) returned nothing on this cluster — likely an implicit date-window default, not confirmed. `qacct -j <jobnumber>` works fine for a known ID. Worked around by extracting all 1,996 real job IDs for every (gene, step) in `run_more_genes_20260817_1143` from local `_log/*/*.out`/`*.err` files (they embed `active_jobs/<jobnumber>.1/...` paths) into `resources/run_more_genes_20260817_1143_job_ids.tsv`, and wrote `resources/pull_qacct_snapshot.sh` to loop `qacct -j` over that list. **User has not run this yet as of session end.**

### 3g. `C7` downgraded
User confirmed compound-het/inheritance phasing now lives in `clinical_variant_prioritization`, not this pipeline — `C7` in `BUG_TRACKER.md` downgraded from Critical to Low (left open, not resolved-by-design, since the dead genotype-request code itself hasn't been removed).

## 4. Files touched this session
- `TODO.md` — new Priority 0 section (provenance manifest, SpliceAI `-D` window, per-predictor resource scaling with `qacct` progress notes); new Phase 5 (SAI-10k-calc).
- `BUG_TRACKER.md` — `ORCH-10` (new, Critical), `ORCH-11` (new, High), `C7` (downgraded Critical → Low).
- `resources/run_more_genes_20260817_1143_job_ids.tsv` — new, 1,996 rows (gene, step, jobnumber).
- `resources/pull_qacct_snapshot.sh` — new, not yet run.
- This file.

## 5. Things NOT to redo
- Don't re-derive DMD's `--raw-dir` — it's `RUNS/predictors_050826/input/by_gene_more_genes` (confirmed by exact row-count match, see §2.2).
- Don't re-extract job IDs from logs — `resources/run_more_genes_20260817_1143_job_ids.tsv` already has all 1,996, cleaned of the `.1` PE-task-suffix parsing bug.
- Don't implement any Priority 0 `TODO.md` item without checking with the user first — this session's explicit, repeated instruction was to track/scope, not implement, until they say go.
