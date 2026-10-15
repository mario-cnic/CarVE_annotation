# Performance review of `annotate_vcf.nf`: where the wall-clock goes

**2026-10-14T16:30+02:00.** Read-only review. No pipeline code or `nextflow.config` was edited, and nothing was run on the cluster. Repo at `05b21a3` (clean apart from the new report folder). Goal: shorter end-to-end wall-clock on SGE with no change to scientific output.

Artifacts: [`reports/perf_review_20261014/`](../reports/perf_review_20261014/) contains:
- `parse_nf_runs.py`: per-task, per-process and critical-path tables from traces or Nextflow logs.
- `simulate_fanout.py`: what-if replay of the predictor fan-out.
- `tasks.tsv`, `per_process.tsv`, `critical_path.tsv`, `run_summary.tsv`, `scenarios_F2005.tsv`.
- `F2005_union_job_ids.tsv`: the job list for `qacct`.
- `manifest.json`.

`*.tsv` is gitignored (`.gitignore:26`), so only the scripts and the manifest would be tracked.

Labels: **FACT** = measured or read now. **INFERENCE** = reasoned from facts. **PROPOSAL** = not applied. **UNVERIFIED** = needs a cluster run.

## Runs analysed and their limits

| Run | Source | Size | Wall-clock |
|---|---|---|---|
| `F2005_union` | `.nextflow.log` + work dirs (no trace; trace was not enabled at commit `48233a3`) | 264 broad chunks, 124 restricted chunks; 23,686,438 table rows × 614 columns | **56.6 h** |
| `S223` | `.nextflow.log.8` + work dirs + qacct subsets | 255 / 118 chunks; mouse Pangolin db, SpliceAI `-D 10000` | 93.9 h |
| `panel8_trio_20261008_NMD2` | trace | 1 chunk per tool | 0.70 h |
| `clinvar_af_fixture_20261014` | trace (first of the three launches) | 1 chunk per tool | 0.11 h |

Limits:
- F2005 and S223 have no `%cpu` or `peak_rss`. Only the S223 qacct subset has CPU and memory figures.
- F2005 ran with `workDir` on `/data_lab_PGP` and before the 2026-10-07 VEP field cleanup. The tail I/O, and the 50.8 GB intermediate TSV, may be smaller at HEAD.
- `/mnt/scratch` (the current `workDir`, `nextflow.config:234`) is **not mounted locally**, so the work dirs of the traced runs could not be read.
- **FACT:** the `hostname` column of both traces is `-`. The trace cannot show node placement, so the TODO's expectation that it does is not met. Node placement here comes from `.command.log`.
- Realtimes from logs are unreliable for tasks shorter than about 1 min: Nextflow detects the start late. `%cpu` is not summarised for tasks under 60 s.

## What counts as "output-neutral"

Byte identity between two runs is already impossible:
- bcftools writes `##bcftools_*Command=...; Date=` header lines.
- The Parquet footer carries a stamp with `launched_at`.

Test used for every proposal below:
- **VCF:** record bodies identical, and header lines identical apart from bcftools command/date and other provenance lines.
- **Parquet:** column names, types and order identical; values identical (NaN = NaN); row order identical; footer metadata keys identical.

Row-group boundaries are reported separately. They change file bytes but not content.

## Critical path (FACT, `critical_path.tsv`)

**F2005, 56.6 h**:

| Step | Time | Note |
|---|---|---|
| Pangolin `chunk_0045` submitted | at 1.07 h | Already 63 min after chunking, because of the queue cap |
| Same chunk, attempt 1 on **c0051-cn1** | **36.1 h** | Killed at `h_rt` |
| Retry on c0051-cn3 | 4.5 h | Fan-out end at **41.7 h** |
| Concat + merge | 0.3 h | |
| `VCF_TO_TABLE` | **3 failed attempts (1.6 + 1.7 + 4.7 h)**, then success (5.9 h) | |
| `TAG_TRANSCRIPT_PRIORITY` | 5 quick failures, then 0.2 h | |
| `ADD_WIDE_GENOTYPES` | 1 failure, then 0.2 h | |

Decomposition of the 56.6 h:

| Component | Hours |
|---|---|
| Successful runtime | 11.1 |
| **Failed attempts** | **44.4** |
| Waiting between steps | 1.1 |
| SGE queue wait | 0.06 |

Failed attempts make up 78% of the critical path.

**S223, 93.9 h**:
- SpliceAI `chunk_0029` was killed at 24 h and then at 48 h on c0053-cn1, and finished in 7 h elsewhere.
- `CONCAT_SPLICEAI` and `MERGE` first attempts hung for 2 h each (MISC-16).
- `VCF_TO_TABLE` failed 2× (3.1 h) before succeeding (6.6 h), and TAG failed 4×.

**Small runs:** one Pangolin chunk dominates (panel8: 39 min of 42 min, about 800 variants). Every other step costs about 15 s of polling and queueing.

**Barrier slack (FACT, `run_summary.tsv`):** in F2005 the broad tier (VEP/Branchpointer concat) finished at 17.1 h, SPiP at 14.4 h, SpliceAI at 27.4 h and Pangolin at 41.7 h. `MERGE_ANNOTATIONS` (join, `workflows/merge_subworkflow.nf:20-24`) waits for the slowest tool, which is the correct dependency: the table needs every tool. The broad tier is never on the critical path.

## Ranked findings

The savings are F2005 wall-clock hours from the replay in `scenarios_F2005.tsv`, unless stated otherwise. The replay reproduces the observed run: 41.6 vs 41.7 h fan-out end, 56.6 vs 56.6 h wall. It assumes tasks start as soon as they are submitted. At 400 slots in use, the observed SGE wait was under 1 min. **Savings are not additive; use the combined rows.**

Output effect columns, using the definition above: **B** = file bytes, **O** = row/record order, **S** = schema (VCF header fields, Parquet columns and types, footer keys), **V** = scientific values. Values: "=" unchanged; "prov" = only provenance header lines differ, as they already do every run.

| # | Finding | Est. saving (F2005) | Effort | B | O | S | V | Validate |
|---|---|---|---|---|---|---|---|---|
| 1a | Exclude c0051-cn1 for SpliceAI/Pangolin | **~23 h** (A→C: 56.6→33.4) | config | prov | = | = | = | cluster run |
| 1b | Attempt-1 `time` limit of 20 h for SpliceAI/Pangolin | ~14 h (A→C3: 42.3); at 16 h only ~12 h (C2: 44.6) | config | prov | = | = | **may change on retried Pangolin chunks** unless threads are pinned first | cluster run |
| 2 | Table steps fail out-of-memory on attempt 1 and re-run the parser each time | **~8.3 h** (A→B); with 1a: 56.6→**25.1 h** (D) | config + process split | = | = | = | = | stub run + cluster |
| 3 | `TAG_TRANSCRIPT_PRIORITY` loads 23.7 M × 614 into pandas (needs 72 GB) | removes 5 retries (~0.4 h incl. queue) + ~10 min/run | small | differ | = | = if footer metadata copied, else `pandas` key lost | = | local fixture + cluster |
| 4 | `filter_variants.py` row-wise `apply` returning `pd.Series` | **~3 h** of the 4.3 h filter stage (INFERENCE) | medium | = | = | = (dtypes at risk) | = (text formatting at risk) | per-function equality on fixtures |
| 5 | Nextflow queue cap of 100 tasks × 4 slots = 400 slots | ≤ 4 h (D→F, **upper bound**) | config | = | = | = | = | needs cluster capacity |
| 6 | Restricted chunk size 8000 | ~2 h after 1a (D→H/I); small runs ~40→~15 min | config | prov | = | = | = (per-variant scoring) | chunk-invariance diff |
| 7 | Concat/merge re-compression and full-decompression parity checks | ~0.2-0.5 h | small | prov | = | = | = | local fixture |
| 8 | Annotation cache (`PLT-040`) | **95% of restricted-tier compute** on family re-runs | high | design | | | | separate design |

### 1. Slow chunks only caught at `h_rt` (resource fit / MISC-17 interaction)

**FACT:** with the human Pangolin db (F2005), healthy chunks are tight:

| Tool | Healthy nodes | c0051-cn1 |
|---|---|---|
| Pangolin | median 6.5-6.7 h, max 7.8 h per 8,000-variant chunk (6 nodes) | median **12.6 h**, max 22.4 h, plus one chunk killed at 36 h |
| SpliceAI | median 4.6-5.0 h, max 7.9 h | max **18.3 h** |

S223 showed c0051-cn1 at 1.1-2.2× too (MISC-17).

**FACT:** the Pangolin limit is `36.h * task.attempt` (`nextflow.config:214-216`) and `process_long` is `24.h * task.attempt` (`:211`). The 36 h was set from mouse-db runtimes (BUG_TRACKER MISC-17, "a guess").

PROPOSAL options. Decision needed.
- **(a)** Add `c0051-cn1` to `long_task_excluded_hosts` (`nextflow.config:171`; mechanism already in `:238-241`). Scheduling only.
  - Replay C/E: fan-out end 41.7→18.4 h, wall 56.6→33.4 h.
  - Excluding the node also removes its slots. E (7 × 64 slots) gives the same result, because the queue cap binds, not the slots.
- **(b)** A lower attempt-1 `time` for Pangolin and SpliceAI as a safety net for nodes nobody has flagged yet. The replay models it with `--h-rt-hours`: an attempt over the limit holds its slot until the limit, then re-runs at the process median.
  - 20 h (2.5× healthy max): wall 42.3 h (C3, −14 h). Only the 36 h Pangolin straggler is cut.
  - 16 h: 44.6 h (C2, −12 h). This is *worse* than 20 h, because it also kills an 18.3 h SpliceAI chunk on c0051-cn1 that would have finished.
  - **Precondition:** pin Pangolin threads first (output note below). (b) creates more retries, and retries currently change the thread count.
- **(c)** Both: (a) gives the saving, and (b) at 20 h limits the damage from the next slow node.

**Output note (FACT + INFERENCE):** Pangolin sets `OMP_NUM_THREADS=${task.cpus}` (`modules/local/pangolin.nf:26`), and `process_long` cpus scale with the attempt (`nextflow.config:209`). A retried chunk therefore runs with 8 or 12 threads; F2005 `chunk_0045` succeeded with 8. Different thread counts can change float summation order. Pangolin rounds to 2 decimals (`src/external/Pangolin-main/pangolin/pangolin.py:192-214`), so a flip is rare but possible.
- This exists **today** and grows with every retry.
- Pinning the thread count to 4 regardless of attempt would make retries deterministic. That is a code change and untested.
- SpliceAI pins TF to 2 + 2 threads (`src/python/annotate_spliceai.py:17-18`) and is not affected.

### 2. First-attempt memory for the table steps, and the parser re-run on each retry

**FACT, F2005 work dirs:**
- `VCF_TO_TABLE` attempt 1 ran `vcf_parser_pysam.py` for 1.6 h (08:26 → 10:01, writing a **50.8 GB** `scratch.tsv`). Then `filter_variants.py` died in `pd.read_csv` with "C error: out of memory" (12 GB).
- Attempts 2 and 3 repeated the 1.6 h parse. Attempt 4 (48 GB) parsed in 1.6 h and filtered in 4.3 h.
- TAG failed at 12/24/36/48/60 GB:
  - attempt 1: `ArrowMemoryError: realloc of size 738197504 failed` (22G h_vmem);
  - attempt 2: `malloc of size 70301348032 failed`;
  - attempt 3: "Unable to allocate 40.6 GiB for an array with shape (230, 23686438) and data type object".
  - It succeeded with 72 GB (6th attempt).
- `ADD_WIDE_GENOTYPES` needed 24 GB.
- S223: VCF_TO_TABLE 2 failures, TAG 4.

PROPOSAL, no output change:
- **(a)** First-attempt memory: `VCF_TO_TABLE` 48 GB, `ADD_WIDE_GENOTYPES` 24 GB, `TAG_TRANSCRIPT_PRIORITY` 72 GB (or 12 GB after #3). Use `withName` blocks so `cpus` stays at 4.
  - The h_vmem pairing (`nextflow.config:241`, `max(12, mem+10)`) then gives 58G/34G/82G automatically.
  - These numbers are F2005-sized (23.7 M rows). A bigger cohort needs more. #3 and #4 are the durable fix.
- **(b)** Split `VCF_TO_TABLE` (`modules/local/vcf_to_table.nf:62-78`, two programs in one task) into a parse process and a filter process, so a filter retry does not repeat the 1.6 h parse.
  - Table bytes unchanged. The TSV becomes a cached task output on scratch.
  - In F2005 this alone would have saved 3 × 1.6 h = 4.8 h of the 8 h.
- Already a TODO item ("Then, in this order: (1) first-attempt memory…"), waiting on Mario's OK.

### 3. `TAG_TRANSCRIPT_PRIORITY`: avoid the pandas round-trip

**FACT** (`src/python/tag_transcript_priority.py`):
- The script reads the whole Parquet into pandas (`:46`), copies it (`:76`), runs a row-wise `df.apply(..., axis=1)` (`:77`) and converts back with `pa.Table.from_pandas` (`:54`).
- The out-of-memory messages above come from the Arrow→pandas object conversion, so vectorizing alone will **not** fix memory.

**FACT, local** (50,750-row fixture: the panel8 output replicated 50×, scratchpad only):
- A vectorized tier gives **identical tiers** at 2.3 vs 50 µs/row.
- A streaming PyArrow version (1 Mi-row batches, `ParquetWriter`) gave the same schema (types, order), identical values in all 585 columns and identical tiers, in 0.8 s vs 5.0 s.
- **Difference:** the `pandas` footer key is dropped and only `ARROW:schema` remains.

PROPOSAL: Arrow streaming version. Peak memory becomes about one batch (INFERENCE: a few GB instead of 72 GB).
- The input footer metadata (`pandas` key plus the provenance stamp) must be copied into the writer, or the footer changes.
- Row groups would follow batch boundaries (now 23 row groups of 1 Mi rows; the same if `batch_size=1<<20`).
- Validate:
  - local: fixture equality including the footer key set;
  - cluster: full F2005-size run, comparing schema and values against the current table.
- Matches TODO item 3.

### 4. `filter_variants.py` row-wise passes (the 4.3 h filter stage)

**FACT, local cProfile** (50,750 rows × 510 columns, 201 s): about 85% of the time is in row-wise passes:

| Function | Location | Time |
|---|---|---|
| `build_priority_tier` | `scoring.py:14`, incl. `build_acmg_criteria` `acmg.py:151` `.iloc` loop | 52 s |
| `build_signed_intron_offset` | `splicing.py:519` | 27.5 s |
| `parse_splicevault` | `splicing.py:290` | 26.5 s |
| `parse_branchpointer` | `splicing.py:409` | 23.3 s |
| `build_qc_status` | `io_qc.py:221` | 22.9 s |
| `parse_maxentscan` | `splicing.py:464` | 21.1 s |
| `analyze_pedigree_inheritance` | `pedigree.py:249` `iterrows` | 6.9 s |

(`splicing.py`, `io_qc.py`, `pedigree.py` are in `src/python/parsing/modules/`.)

Don't extrapolate the 201 s linearly: that gives ~26 h for F2005 against the 4.3 h observed, because WGS rows are mostly empty. Use only the shares.

**FACT, micro-benchmark** on the same frame: a row function returning `pd.Series([...])` costs **153 µs/row**, the same function returning a tuple costs 35 µs/row, and a tuple over 5 columns costs 5 µs/row.
- `build_qc_status` with only its 5 columns gave identical output but was barely faster (8.3 → 7.9 s). The cost is the per-row `pd.Series` construction (`io_qc.py:171,217,219`; `splicing.py:260-287,406,450-461,488-517`), not the row width.

INFERENCE: at 23.7 M rows each such pass costs about 1 h. Returning tuples over the needed columns, or vectorizing, would cut the filter stage from 4.3 h to roughly 1 h.

Risk: result dtypes and string formatting can change when `pd.Series` returns become tuples or when numeric columns are upcast in a subset. Each function needs its own equality test (values and dtypes) on fixtures that include missing values, multiallelic-split rows and indels. Effort: medium.

Second part of TODO item 3: drop the 50.8 GB TSV between parser and filter. It saves I/O plus a CSV parse of 50 GB, but needs a schema-preserving reader. This was not measured.

### 5. Queue cap (scheduler overhead)

**FACT:**
- No `executor` block in `nextflow.config`. The log shows "capacity: 100; pollInterval: 5s".
- F2005 kept 100 tasks in SGE for 14.3 h (S223: 23.7 h). Peak slots in use were exactly 400 in both runs: 56-60 per host on 8 hosts (F2005), 64 per host on 7 hosts (S223).
- Restricted chunks were submitted up to 12 h after chunking (median 4.6-5.4 h); VEP chunks were submitted 13.6 h after chunking.
- SGE queue wait itself was 5-35 s per task.

INFERENCE: we may already use most of the queue (8 × 64 = 512 slots if the nodes have 64 slots, which is unknown). Replay F (queue 128, 512 slots) gives 25.1 → 21.2 h; G (800 slots) gives 18.5 h. Both are **upper bounds** until capacity is known.

Cheaper lever, the same throughput with fewer slots:
- SpliceAI requests 4 slots and uses 2.63 cores (median, S223 qacct). Requesting 3 frees 25% of SpliceAI slot-hours (~580 slot-h in F2005).
- Branchpointer requests 4 slots and 12 GB for seconds of work: 264 tasks in F2005, `modules/local/branchpointer.nf:16` (`process_medium`).

Risk: **h_vmem.** The panel8 trace shows SpliceAI `peak_vmem` 29.6 GB under `h_vmem=22G`. F2005 TAG failed a 70 GB malloc under `h_vmem=34G` with 8 slots. These two facts conflict on whether `h_vmem` is per slot or per job, so do **not** lower SpliceAI slots or memory before `qconf -sc` shows how `h_vmem` is consumed (commands below).

### 6. Restricted chunk size (chunking)

**FACT:**
- `gene_restricted_chunk_size = 8000` (`nextflow.config:36`), `chunk_size = 20000` (`:29`).
- `chunk_threshold = 15000` (`:30`) is **not referenced anywhere** in `annotate_vcf.nf`, `workflows/`, `modules/` or `src/python/`: a dead parameter.
- `spip_max_lines = 22000` (`:53`) is larger than the 8,000-record chunk, so SPiP's internal batching never splits a chunk (INFERENCE from the option name; not read in the R script).
- Per-chunk variance is low on healthy nodes (Pangolin p95/median 1.2 in F2005 excluding c0051-cn1). VEP chunks vary 7 min (median) to 58 min (p95), but VEP is never on the critical path.

INFERENCE: once #1 is applied, halving restricted chunks gains about 2 h on a WGS-size run. Replay H/I (2 or 4 min fixed cost per task): 25.1 → 23.0-23.3 h. The fixed cost per near-empty task is 38-58 s for Pangolin and 44 s-4 min for SpliceAI (fixture traces).

The larger effect is on **small validation runs**. panel8's single Pangolin chunk took 39 min; splitting the restricted set into at least 4 chunks would bring the run to about 15 min (INFERENCE).

Neutral: each variant is scored independently, and `CONCAT_CHUNKS` sorts (`modules/local/chunk.nf:106`). Validate with a record-body diff of the predictor VCFs at two chunk sizes on `panel8_trio`.

### 7. Redundant work in concat/merge (tool level and I/O)

**FACT:**
- `CONCAT_CHUNKS` re-indexes every chunk (`chunk.nf:101-103`), even though VEP/Pangolin/SPiP already wrote a `.tbi` that the wiring drops (`workflows/broad_pass_whole.nf:60`).
- It then runs `concat -a` and a full `bcftools sort` (`:104-106`), and counts records by decompressing both whole files (`:109-110`).
- `MERGE_ANNOTATIONS` runs 4 `bcftools annotate` passes, each fully compressing and indexing the 3 GB VCF (`modules/local/merge.nf:50-60`), with no `--threads`.
- Measured: `CONCAT_VEP` 12 min, `MERGE` 17 min on F2005. MERGE is on the critical path.

PROPOSAL:
- Pipe the merge with `-Ou` and compress once with `--threads`.
- Count records with `bcftools index -n` on the existing indexes.

The records would be identical; only the bcftools command header lines change, as they already do every run. Saving about 15-25 min per WGS run (INFERENCE). Validation: local fixture, record bodies and non-provenance header lines.

### 8. Annotation cache (`PLT-040`)

**FACT:** 936,412 of F2005's 984,004 restricted-tier variants (**95.2%**) are exactly the S223 restricted set, compared by CHROM/POS/REF/ALT. S223 is the WGS member of the union (1 shared sample column). The family run therefore re-scored its own member's variants with SpliceAI/Pangolin, about 5,600 slot-hours.

Not reusable today: S223 was scored with SpliceAI `-D 10000` and the mouse Pangolin db.

INFERENCE: a cache would need a key of (CHROM, POS, REF, ALT, tool version, model weights, `-D`, FASTA, Pangolin db, gene-restriction BED) and the invalidation work of TODO item 6. Cross-individual reuse between unrelated samples was **not measured**. Report and design separately; do not fold it into this pass.

## Other dimensions (one line each)

- **DAG shape:** no stage waits on something it does not need. The only real barrier is the `MERGE` join. Splitting the merge per chunk so parsing starts early would be high effort for little gain, because Pangolin decides the end anyway.
- **Retry strategy:** `maxRetries = 5` with linear scaling (`nextflow.config:190-191`). Memory retries go up in 12 GB steps, so TAG needed 6 attempts. Time retries double `h_rt`, which doubles the waste on a slow node (S223: 24 h + 48 h on c0053-cn1).
- **MISC-16:** F2005 `CONCAT_VEP` attempt 1 hung for 1 h on c0053-cn1 (not on the critical path there). In S223 two 2 h hangs were on the critical path. The 1 h cap holds; the cause is still unknown.
- **VEP:** `--fork ${task.cpus}`, default `--buffer_size`, offline cache, gnomAD v4.1 `--custom` with about 80 fields (`modules/local/vep.nf:73-101`). Plugin cost cannot be profiled locally (no cache mirror). VEP is not on the critical path, so this is low priority.
- **SpliceAI/Pangolin internals:**
  - SpliceAI makes 10 Keras `predict()` calls per allele and gene (`/data_lab_PGP/shared/utils/conda_envs/spliceai_env/lib/python3.10/site-packages/spliceai/utils.py:159-160`).
  - Pangolin runs 24 forward passes per variant and strand, ref and alt separately (`pangolin.py:54-59`).
  - Batching would change batch shapes and possibly float bits (scores are rounded to 2 dp). **Out of scope** unless an equality test on a fixture passes.
  - The VEP SpliceAI plugin's precomputed scores use a different distance, so using them to skip the model would **change values**. Out of scope.
- **SPiP/Branchpointer:** SPiP median 23 min/chunk, 8.7 GB RSS (panel8). Branchpointer takes seconds. Neither is on the critical path.
- **I/O:**
  - `workDir` is now on `/mnt/scratch` (`nextflow.config:234`). F2005's four `VCF_TO_TABLE` attempts each wrote a 50.8 GB TSV under `/data_lab_PGP`.
  - `TMPDIR=/data_tmp` is unreliable (BUG_TRACKER, CONCAT sort temp dir), and `bcftools sort -T ./` is fixed.
  - `publishDir mode: 'copy'` copies about 15 GB at the end. That is needed for `-resume`; keep it.
  - Parquet: SNAPPY, 1 Mi-row groups. bgzip is single-threaded except in Pangolin.
- **Resume:**
  - The stamp path `pipeline_info/TABLE_STAMP.json` is fixed per `run_id` (`run_annotate_vcf.sh:113`), but its content has `launched_at`. It is passed as a string (`modules/local/tag_transcript_priority.nf:26-32`, and similarly in `genotypes.nf`).
  - INFERENCE: on `-resume` the cached TAG/genotype outputs keep the **previous** launch's stamp. That is a provenance mismatch, not a speed problem. Add it to TODO item 6.
  - `meta.run_id` is part of every task hash, so resuming under another `--run_id` re-runs everything.
  - `collectForConcat` uses `groupTuple()` with no `sort` (`workflows/broad_pass_whole.nf:30-31`, `workflows/gene_restricted_whole.nf:33-34`), so each `CONCAT_*` gets its chunk list in completion order.
    - INFERENCE (Nextflow documents `groupTuple`'s `sort` option for exactly this): a resumed run whose chunks finish in a different order changes the input hash, so the `CONCAT_*` steps and everything downstream re-run. `walkthrough/20260925_resume_cache_invalidation_and_fix.md` does not cover this; not verified here.
    - PROPOSAL: `groupTuple(sort: 'hash')`. A name sort does not work, because the chunk outputs share one file name. Output: records unchanged (`CONCAT_CHUNKS` sorts); only the bcftools command header line lists the files in another order.
- **MISC-17 (open, not re-analysed):**
  - c0053-cn1 is excluded for long tasks (applied 2026-10-02). There is still no trace evidence that the exclusion took effect, since `hostname` = `-`.
  - The IT answer is pending, and the Pangolin `h_rt` is still the mouse-db guess.
  - c0051-cn1 is not excluded and is slow in both runs (see #1).

## Decisions needed (Mario)

1. Slow node (#1):
   - exclude c0051-cn1 (1a); and/or
   - a 20 h attempt-1 limit (1b), only after Pangolin threads are pinned to 4 on every attempt. Pinning is a code change; without it, retried chunks may get different scores.
2. First-attempt memory for the 3 table steps plus splitting `VCF_TO_TABLE` (#2). This is the existing TODO item; confirm.
3. TAG Arrow rewrite: keep the `pandas` footer key (copy the metadata) or accept dropping it (#3).

## Cluster commands (UNVERIFIED; for the login node, nothing below has been run)

```bash
# Capacity: slots per host, queue limits, how h_vmem/h_rss/mem_free are consumed (per slot = YES, per job = JOB)
qhost -q
qstat -g c
qconf -sq BIG_NEW.q | grep -E '^(hostlist|slots|h_vmem|h_rt|s_rt)'
qconf -sc | grep -E 'h_vmem|h_rss|mem_free|slots'
qquota -u "$USER"
# Per-slot or per-job h_vmem, directly (5 min, 4 slots, 4G): prints the address-space limit the job sees
qrsh -pe smp 4 -l h_vmem=4G,h_rt=00:05:00 -P BIGN -A PGP bash -c 'ulimit -v'
# CPU and memory of every F2005 task (925 jobs, several minutes), then rebuild the tables with it
bash resources/pull_qacct_snapshot.sh reports/perf_review_20261014/F2005_union_job_ids.tsv RUNS/perf_review_20261014/F2005_qacct.tsv
# (local, afterwards) the run folder's copy of the log is identical to the rotating repo-root .nextflow.log used here
python3 -I reports/perf_review_20261014/parse_nf_runs.py --log F2005_union nf_work/annotation_out/F2005_union/pipeline_info/nextflow_2026-10-04_23-05.log --qacct RUNS/perf_review_20261014/F2005_qacct.tsv --outdir RUNS/perf_review_20261014/with_qacct
```

## Proposed `TODO.md` entries (not applied)

- [ ] **Long-task straggler control** (perf review 2026-10-14 #1): c0051-cn1 exclusion; optionally an attempt-1 `time` of 20 h for SpliceAI/Pangolin, after Pangolin threads are pinned. Healthy human-db Pangolin max is 7.8 h per 8,000 variants. Verify on the next WGS-size run: no long task over 20 h on a healthy node, `.command.run` carries the exclusion.
- [ ] **Table steps: first-attempt memory and split `VCF_TO_TABLE`** into parse + filter (perf review #2). Merge with the existing "first-attempt memory" item.
- [ ] **TAG via PyArrow streaming** (perf review #3; part of item 3): footer metadata copied, local fixture equality (schema, values, row order, footer keys), then F2005-size cluster comparison.
- [ ] **`filter_variants.py`: replace per-row `pd.Series` returns** in `io_qc.build_qc_status` and `splicing` parsers, plus the `acmg`/`pedigree` loops (perf review #4; part of item 3). One function at a time, each with a values+dtypes equality test.
- [ ] **Cluster capacity and h_vmem semantics** (perf review #5): run the `qconf`/`qrsh` commands, then decide queueSize, SpliceAI slots (3) and Branchpointer label.
- [ ] **Restricted chunk size / minimum chunk count** (perf review #6): chunk-invariance diff on `panel8_trio`; remove the unused `chunk_threshold` param.
- [ ] **Concat/merge streaming** (perf review #7): `-Ou` pipe, `--threads`, `bcftools index -n` parity.
- [ ] **Resume keeps a stale table stamp; `groupTuple` order breaks the CONCAT cache** (perf review, resume): add to item 6; `groupTuple(sort: 'hash')`.
- [ ] **Pangolin threads independent of attempt** (perf review #1 output note).

## Implementation — 2026-10-14T17:12:05+02:00 (Mario's decisions after the review)

Decisions taken:
- exclude c0051-cn1 (1a, no time-limit change);
- first-attempt memory and the `VCF_TO_TABLE` split (#2);
- TAG via PyArrow with the footer metadata copied (#3).

Also asked for: the `filter_variants` loops (#4), `groupTuple(sort: 'hash')`, concat/merge streaming (#7), and removal of `chunk_threshold`.

**Not done** (no decision or out of scope): pinning Pangolin threads, the 20 h limit, queueSize/slot changes, chunk size, the stale-stamp-on-resume issue.

### What changed

| Area | Change |
|---|---|
| `nextflow.config` | `long_task_excluded_hosts = ['c0053-cn1', 'c0051-cn1']`; `TSV_TO_TABLE` 48 GB (+24 GB per retry), 12 h × attempt; `ADD_WIDE_GENOTYPES` and `TAG_TRANSCRIPT_PRIORITY` 24 GB (+12 GB per retry); `chunk_threshold` removed. The `TSV_TO_TABLE` override also applies to the legacy per-gene `main.nf` (48 GB per gene task; that path is idle) |
| `modules/local/vcf_to_table.nf` | `VCF_TO_TABLE` is now a workflow running `VCF_TO_TSV` (parser) → `TSV_TO_TABLE` (filter). Same commands and environment; a filter retry reuses the cached TSV. Callers (`annotate_vcf.nf`, `main.nf`) are unchanged. Trace names: `VCF_TO_TABLE:VCF_TO_TSV`, `VCF_TO_TABLE:TSV_TO_TABLE` |
| `workflows/{broad_pass_whole,gene_restricted_whole}.nf` | `groupTuple(sort: 'hash')` in `collectForConcat` |
| `modules/local/chunk.nf` | `concat -a … -Ou \| sort -Oz` (no intermediate compressed file); parity counts from `bcftools index -n`, falling back to a full count; `set -o pipefail` |
| `modules/local/merge.nf` | intermediates are `-Ob1 -W` BCF, indexed on write, instead of VCF.gz + `tabix`. A pipe is not possible because `annotate -a` needs an indexed target (`Could not retrieve index file for '-'`, bcftools 1.21 and 1.24) |
| `src/python/tag_transcript_priority.py` | Parquet → Parquet streams one row group at a time with PyArrow. The tier is computed with Arrow compute using the same rule. The footer is the `pandas` entry `from_pandas(preserve_index=False)` would write, plus the stamp. The TSV path uses a vectorized `tag_priority` |
| `src/python/parsing/modules/{io_qc,splicing,pedigree,acmg,scoring,disease_hpo}.py` | New helpers `io_qc.row_records` / `apply_rows`. Row functions get a dict of only the columns they read, boxed as `apply(axis=1)` boxes them, return tuples, and results get the same `infer_objects` step. The `.iloc[i]` loops read `.to_numpy()` arrays. Row logic is unchanged |

### Checks (local only)

**Concat** — FACT, 78-record panel8 VEP output split into 4 shuffled chunks:
- old vs new commands give identical record bodies (also identical to the input), and identical headers apart from `##bcftools_*` lines;
- holds with bcftools 1.21 (pysam 0.23 bundle, the cluster version) and 1.24;
- a lost chunk still fails parity.

**Merge** — FACT, panel8 predictor outputs:
- identical record bodies, equal to the published `panel8_trio_20261008_NMD2.annotated.vcf.gz`; headers differ only in `##bcftools_annotateCommand`; same on 1.21 and 1.24;
- on a 200× synthetic (15,600 records): 13.6-14.6 s → 8.0-9.0 s.

**TAG** — FACT, compared with the HEAD script on 4 inputs × with/without stamp (50,750-row panel8 replica, `tests/fixtures/sample_variants.pq`, a 2.5 M-row synthetic with 3 row groups and null/NaN/whitespace/`.`/version edge values, a float `MANE_SELECT` with no `Feature`):
- all 8 have identical footer bytes, schema with metadata, row groups and values;
- 2.5 M rows: 22.7 s / 1.58 GB → 3.0 s / 0.38 GB.
- An ordered categorical (`PRIORITY_TIER`-like, one row group with a single category, one category unused) keeps its categories and order. The output is cut into 1 Mi-row groups exactly as `pq.write_table` cuts it, whatever the input layout.
- **Full scale** (the published F2005 table, 23,686,438 rows, read-only, output written to `/dev/null`):
  - the recomputed tier differs from the stored column on 0 rows; counts {1: 38,535, 2: 2,144,530, 3: 21,503,373} equal the old run's log;
  - the streaming writer peaks at **11.2 GB RSS** in 9.1 min. One 1 Mi-row group is 5.5-7.7 GB in Arrow, and keeping the old row groups needs one in memory, so TAG gets 24 GB, not 12.

**`filter_variants.py`** — FACT, old (HEAD snapshot) vs new on 5 inputs (ClinVar/AF fixture, panel8 without GT, with one and three synthetic `GT_` columns, the 50,750-row replica):
- identical schema, footer, column order, row order and values;
- 50,750 rows: 72.9 s → 22.0-25.6 s.

**Note (FACT):** production parses without FORMAT fields, so the table has no `GT_` columns; the multi-sample inheritance branch is only covered by the synthetic input.

**Resume and `groupTuple`:** that `sort: 'hash'` keeps `CONCAT_*` cached is INFERENCE from Nextflow's documented `groupTuple` behaviour. The local stub resume below would also have passed without the change, because stub chunks finish in a stable order.

**Nextflow** — local `-stub-run` (Nextflow 24.04.2):
- `annotate_vcf.nf` with 20-record chunks: 4 chunks per tool, real concat + parity = 78 for all five tools; a `-resume` keeps all tasks cached; resources resolve as set (`TSV_TO_TABLE` 48 GB, `ADD_WIDE_GENOTYPES` 24 GB);
- legacy `main.nf` (MYBPC3) runs through the new workflow, with one harmless warning that no process matches the `ADD_WIDE_GENOTYPES` selector.

**Tests:** new `tests/python/test_tag_transcript_priority.py` and `test_row_apply.py`. All modules pass except the two that already failed before these changes from local environment issues: `test_combine_valencia_hic.py` (`genebe` missing) and `test_write_run_manifest_releases.py` (sqlite symbol).

### Not verified (needs a cluster run)

Before a WGS run, check the two-host exclusion string on the login node (it has never been rendered by the `standard` profile or submitted): `qsub -w v -P BIGN -A PGP -l h='!c0053-cn1&!c0051-cn1' -b y /bin/true`. **Done 2026-10-15T10:00:43+02:00 (Mario, login node): `verification: found suitable queue(s)`.** Still to see: the string as rendered in a real chunk's `.command.run`.

Output check: compare against the **same input run at `d5e6461`** (code before these changes), not against an older run such as F2005, which predates the VEP field cleanup and the ClinVar/AF changes.
- Cheapest: two `panel8_trio` or ClinVar-fixture runs (about 8 min each), one per commit, compared with `reports/perf_review_20261014/compare_tables.py A.pq B.pq --ignore-footer-key carve.provenance` (schema, footer, row groups, column and row order, values).
- At WGS scale: rerun only `TSV_TO_TABLE` and TAG with both code versions on one cached `scratch.tsv`.

`MERGE` now writes level-1 BCF intermediates, probably larger than the old VCF.gz ones; this matters only for scratch quota.

**Cluster panel run 2026-10-15T10:19:49+02:00 (FACT):** `nf_work/annotation_out/panel8_perf_20261015`, launched from `afd1abf` on a clean tree.
- Exit 0, 20 tasks, 0 retries. `VCF_TO_TABLE:VCF_TO_TSV` and `VCF_TO_TABLE:TSV_TO_TABLE` ran; `TSV_TO_TABLE` requested 48 GB, TAG and `ADD_WIDE_GENOTYPES` 24 GB.
- Compared with `panel8_trio_20261008_NMD2` (commit `1bcf165`; the only other code change since then is `e1f8fa2`):
  - main table, `withGT` table and genotype table identical (schema, footer without the stamp, row groups, column and row order, values; `compare_tables.py`);
  - annotated VCF and the 5 predictor VCFs have identical records (78 each); headers differ only in timestamps, work-dir paths and the order of the VEP version list;
  - `e1f8fa2` changed no ClinVar value on this input.
- Host exclusion confirmed 2026-10-15T10:54:46+02:00 (Mario, login node): the SpliceAI and Pangolin chunks' `.command.run` (`96/5f8a41`, `57/0f1a93`) carry `-l h_vmem=22G -l h=!c0053-cn1&!c0051-cn1` after the unchanged `-v TMPDIR=/data_tmp -S /bin/bash -P BIGN -A PGP`.

On the next WGS-size run, check:
- `.command.run` of SpliceAI/Pangolin chunks has `-l h=!c0053-cn1&!c0051-cn1`;
- `TSV_TO_TABLE` succeeds on attempt 1 within 48 GB and 12 h;
- TAG fits 24 GB;
- `MERGE` and `CONCAT_*` times drop and their parity messages look right;
- the `-W` CSI write works with the cluster's bcftools 1.21 binary (the local test used the pysam 1.21 build);
- the final table is unchanged against the previous run of the same input: schema, footer keys, values, row order.

Expected saving on an F2005-size run, per the replay (scenario D): 56.6 → about 25 h. INFERENCE until measured.
