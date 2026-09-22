# Nextflow migration — Phase 1 (SpliceAI spike) + Phase 2 (Pangolin, SPiP)

**2026-09-21** — filename kept from Phase 1 (same day, same migration); this file now covers both
Phase 1 and Phase 2, see the "Phase 2" section below.

## Context

Started from a request to annotate a WGS VCF sample. Investigation (this session) found `main.sh`
is architected around a single-gene-per-file model (`gene_name` parsed from the input filename,
`main.sh:363-364`) with no genome-wide/WES execution path. Follow-up research established:

- SpliceAI and Pangolin have no usable genome-wide precomputed score alternative at this
  pipeline's target window (~4999-10000bp) — verified via web search: Illumina's retired SpliceAI
  precomputed table and the public Pangolin Zenodo set (DOI 10.5281/zenodo.15649338) both used a
  50bp window. Both tools must stay gene/region-restricted. Branchpointer/LaBranchoR
  (`annotate_branchpointer_vars.sh:43`, a lookup against a precomputed genome-wide catalog) and
  VEP's whole plugin stack are broad-safe by design.
- No `h_vmem`/thread number anywhere in this pipeline is a trustworthy cost signal — confirmed by
  the user (2026-09-21 correction), these were all set manually with no profiling.
- `carve-platform/ROADMAP.md`'s D-3 ("migrate Module 2 to Nextflow, phased") has a reviewer's lean
  but a **blank Verdict column** (`carve-platform/REVIEW-2026-09-18.md:115`) and is listed as
  **still open** in `carve-platform/STATUS.md:17` as of this writing — this work is the first
  concrete step toward that lean, not an already-approved platform decision.

Full architecture design and the 8-phase migration sequence: see
`/home/mruizp/.claude/plans/scalable-wibbling-snowflake.md` (the approved plan document; not
copied into this repo, referenced by path since it lives in the user's Claude Code plans
directory).

## What was built (Phase 1 only)

New files, none of them touching or replacing `main.sh`/`src/hpc/*.sh`:

- `nextflow.config` — `standard` profile targets the real SGE cluster (production default,
  confirmed by user: the pipeline always runs there, never locally); `local_dev` profile is a
  sandbox-only convenience pointing at a separate local miniforge3 install, **not** a supported
  deployment mode.
- `main.nf` — entry point; legacy per-gene adapter (mirrors `main.sh:363-364`'s filename-derived
  gene identity), requires `--input_vcf <gene>.vcf.gz` with a sibling `.tbi`.
- `modules/local/spliceai.nf` — `SPLICEAI_CHUNK`, `SPLICEAI_ANNOTATE`, `SPLICEAI_CONCAT`, wrapping
  the existing, unmodified `src/python/split_vcf_chunks.py` and `src/python/annotate_spliceai.py`.
- `workflows/spliceai_spike.nf` — wires the three processes with native `.flatMap()`/`.groupTuple()`
  fan-out/collection, replacing the hand-rolled bash worker-pool in
  `annotate_spliceai_vars.sh:120-152`.

## Verification — what was and wasn't actually checked

**Correction (2026-09-21, same session)**: this section originally claimed no GRCh38 FASTA exists
anywhere in this sandbox. That was wrong — a first `ls`/`find` against
`/home/mruizp/data_references/` came back empty (likely a transient NFS stall on first access),
and that empty result was reported as fact without a second check. The user caught it. Re-checked:
`/home/mruizp/data_references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta`
is real, 3.2GB, readable, with a valid `.fai`.

**Could not verify, for a different and now-accurate reason**: real SpliceAI model-score parity
against the known-good `test_data/test_run/annotation/MYBPC3.annSpliceAI.vcf.gz`. With the real
FASTA wired in, `SPLICEAI_ANNOTATE` was run for real (no `-stub-run`) and failed —
`AttributeError: 'Series' object has no attribute 'get_values'` inside the vendored `spliceai`
package, because this sandbox's local (non-NFS) miniforge3 mirror of `spliceai_env` has
`pandas 2.2.3`, incompatible with whatever older `spliceai` build is installed there. This is a
stale **local-only** environment mismatch, not a pipeline or Nextflow-port bug — confirmed by
checking a real prior production run
(`RUNS/run_20260813_1028/annotation/MYBPC3.annSpliceAI.vcf.gz`), which has genuine non-trivial
`SpliceAI=...` scores, proving the production `spliceai_env` (on the NFS-mounted, non-executable-
from-here `shared/utils/conda_envs/`) works correctly. Real model-output parity still has to be
checked when this runs on the actual cluster, against that working environment — `-stub-run`
remains the practical way to validate the DAG/chunking mechanics from this sandbox.

**Side finding from the failed real run**: `errorStrategy = 'retry'` (both here and in the sarek
pattern this design copies) retries *any* process failure with `task.attempt`-escalated resources,
including a deterministic code/dependency bug that can never succeed regardless of retries — 5
wasted retries (up to 72GB requested, exceeding this machine's 62.4GB) happened before the error
surfaced. Worth tightening later to only retry on OOM-like exit codes (e.g.
`errorStrategy = { task.exitStatus in [137,140,143] ? 'retry' : 'terminate' }`), not applied in
this Phase-1 config since it wasn't part of the original scope — noting it here so it isn't lost.

**Was verified for real, locally** (using a separate local, non-NFS miniforge3 install —
`/home/mruizp/data_lab_PGP` is mounted `noexec`, confirmed this session, so the production conda
envs under `shared/utils/conda_envs/` cannot execute here at all): `SPLICEAI_CHUNK` and
`SPLICEAI_CONCAT` ran their real scripts; `SPLICEAI_ANNOTATE` ran via Nextflow's native `stub:`
block (`-stub-run`) in place of the unavailable model/FASTA. Against
`test_data/test_run/_tmp/MYBPC3.vcf.gz` (22 records):
- Single-chunk path (default `--chunk_size 20000`): 22 in, 22 out, record-count parity confirmed.
- Forced multi-chunk path (`--chunk_size 5` → 5 chunks): 22 in, 22 out, parity confirmed — proves
  the fan-out/collect/sort/concat mechanics, not just the trivial single-chunk case.

## Three real bugs the spike caught (the point of doing this as a spike first)

1. **Scalar-vs-list Path unwrapping.** A `path(glob)` process output emits a scalar `Path`, not a
   `List`, when the glob matches exactly one file — and `java.nio.file.Path` implements
   `Iterable<Path>`, iterating its path *segments*. The first run (single-chunk MYBPC3 case) fanned
   out into 13 bogus "chunks" named after directory components (`home`, `mruizp`, `work`, ...)
   instead of the one real chunk. Fixed in `workflows/spliceai_spike.nf` by normalizing to a list
   before `.collect()`. This is a general Nextflow gotcha, not specific to this pipeline — worth
   remembering for every later phase that fans out over a glob output.
2. **Retry cascade hit a local CPU ceiling.** The bogus fan-out above meant 13 failing tasks
   retrying with `cpus = 4 * task.attempt`, which at attempt 4 requested 16 cores against a
   12-core machine and hard-errored ("Process requirement exceeds available CPUs"). Capped
   `process_medium`'s cpus at `Math.min(4 * task.attempt, 12)` in `nextflow.config` as a defensive
   fix independent of bug 1 (a real gene could legitimately need enough retries to hit a similar
   ceiling on a real node).
3. **`bcftools concat -a` needs per-input indices.** Initially matched
   `annotate_spliceai_vars.sh:166`'s use of `-a`/`--allow-overlaps`, which requires each input
   tabix-indexed — none of the chunks were indexed, since nothing downstream needs random access
   into them. Since chunks are non-overlapping by construction (sequential record-count split, not
   a region split), dropped `-a` rather than adding indexing that serves no purpose for this data
   shape — plain `bcftools concat` doesn't require input indices.

## What's next (not done in Phase 1)

Per the plan's migration sequencing: Phase 2 (Pangolin + SPiP at `restricted` tier default) is the
next slice, still validated on the legacy per-gene adapter before any WGS-entry-point work starts
(Phase 7). Real model-score parity for this Phase 1 slice should be checked the next time this
runs on the actual cluster. `carve-platform`'s D-3 Verdict and a Module 3 decision-record note (per
this repo's own CLAUDE.md rule) are flagged in the plan as follow-up, not done in this session.

## Phase 2 — Pangolin + SPiP, gene-restricted tier complete

**2026-09-21, same session**

### What changed vs. Phase 1

Unlike SpliceAI, `main.sh` runs Pangolin and SPiP directly on the whole per-gene VCF, with no
chunking stage (`main.sh:420-429` for SPiP, `main.sh:458-463` for Pangolin — both take the same
`$vcf_gz` that SpliceAI's chunker also starts from). So Phase 2 is two new single-process modules,
not two more chunk/fan-out subworkflows:

- `modules/local/pangolin.nf` — `PANGOLIN_ANNOTATE`, ported from
  `src/hpc/annotate_pangolin_vars.sh`. Wraps `python3 -m pangolin.pangolin` against the vendored
  `src/external/Pangolin-main` (via `PYTHONPATH`, matching the bash script's own approach) and the
  ~1GB `pangolin_grch38.db` under `shared/utils/pangolin_db/`. `-d 10000` kept verbatim (non-goal,
  same treatment as SpliceAI's `-D`).
- `modules/local/spip.nf` — `SPIP_ANNOTATE`, ported from `src/hpc/annotate_spip_vars.sh`. Wraps
  `Rscript SPiPv2.1_main.r` (an external script under `/data_lab_PGP/resources/annotation/SPiP/`,
  not vendored inside this repo) and reproduces the bash script's post-processing step that strips
  SPiP's non-standard `##SPiP output v2.1` header line before bgzip/tabix. `--maxLines 22000` kept
  verbatim (non-goal — internal detail of the R script, not orchestration-level chunking).
- `workflows/gene_restricted.nf` — new `GENE_RESTRICTED_SUBWORKFLOW`, replacing `main.nf`'s direct
  call to `SPLICEAI_SPIKE`. Runs SpliceAI, Pangolin, and SPiP **in parallel off the same
  `meta_vcf_ch`**, matching `main.sh` submitting all three as independent qsub jobs with no
  dependency between them (not a pipeline chain). SPiP is gated behind
  `params.spip_tier == 'restricted'` (the default) — a future `'broad'` placement would move
  `SPIP_ANNOTATE` into `BROAD_PASS_SUBWORKFLOW` (Phase 3, not built yet) instead.
- `nextflow.config` — added `pangolin_python`/`pangolin_repo`/`pangolin_db`/`pangolin_distance`,
  `spip_rscript`/`spip_script`/`spip_max_lines`, and `spip_tier` (default `'restricted'`). All
  cluster-absolute paths confirmed real on the NFS mount this session (`pangolin_grch38.db` is
  961MB at `shared/utils/pangolin_db/`, `SPiPv2.1_main.r` is real at
  `/data_lab_PGP/resources/annotation/SPiP/`). `outdir` default renamed
  `spliceai_spike_out` → `gene_restricted_out` since it now holds all three predictors' output.

### Verification

Same constraint as Phase 1, one step further: **neither Pangolin nor SPiP has any local execution
path at all**, not even a version-mismatched one — unlike SpliceAI's `spliceai_env`, this sandbox's
local miniforge3 mirror has no `pangolin_env` or `spip_env` (confirmed by listing
`/home/mruizp/apps/miniforge3/envs/`). Both processes are `-stub-run`-only in this sandbox; real
model-output parity can only be checked on the cluster.

Ran `nextflow run main.nf -profile local_dev -stub-run --input_vcf test_data/test_run/_tmp/MYBPC3.vcf.gz`
twice (once before, once after the `outdir` rename). Both runs: exit 0, all three predictors
(`SPLICEAI_SPIKE`, `PANGOLIN_ANNOTATE`, `SPIP_ANNOTATE`) completed `1 of 1 ✔` in parallel, and all
three published outputs landed under `nf_work/gene_restricted_out/MYBPC3/`:
`MYBPC3.annSpliceAI.vcf.gz(.tbi)`, `MYBPC3.annPangolin.vcf.gz(.tbi)`, `MYBPC3.annSPiP.vcf.gz(.tbi)`.
This confirms the DAG/channel wiring (parallel fan-out from one shared input channel, `spip_tier`
conditional gating, `publishDir` per-gene layout) — not model-score correctness, which remains
cluster-only exactly as in Phase 1.

The stub run verified DAG wiring, parallel fan-out from one shared input channel, `spip_tier`
conditional gating, and per-gene `publishDir` layout — it did not execute either new `script:`
block (both stubs are a plain `cp`), so it says nothing about whether those command lines actually
run or write the filenames the `output:` blocks declare. Checked that separately by reading the
two external tools directly (no execution needed):
- `src/external/Pangolin-main/pangolin/pangolin.py:251` treats its 4th positional arg as an output
  *prefix*, appending `.vcf` only if the given name doesn't already end in `.vcf`
  (`annotate_pangolin_vars.sh`'s own `if [ -f "$RAW_OUT" ]` guard around bgzip is a tell that this
  ambiguity bit someone before). `PANGOLIN_ANNOTATE`'s `raw_out` already ends in `.vcf`, so no
  double-extension — confirmed by reading the source, not by running it.
- `SPiPv2.1_main.r:318` opens `outputFile` (the literal `--output` argument) with `file(...)` and
  writes straight to it — no derivation or extension-appending. `SPIP_ANNOTATE`'s `--output` arg
  matches its declared `output:` path exactly.

Neither real `script:` block has been executed anywhere yet, in this sandbox or elsewhere — that
still has to happen on the cluster.

## What's next (not done in Phase 2)

Per the plan's migration sequencing: Phase 3 (broad tier — VEP containerized, Branchpointer) is
next, still on the trivial single-gene-sized "region" case. Real model-output parity for both
Phase 1 (SpliceAI) and Phase 2 (Pangolin, SPiP) still needs to be checked together the next time
this runs on the actual cluster — there is no local way to validate any of the three predictors'
actual scores, only their orchestration.

## Phase 3 — VEP (containerized) + Branchpointer, broad tier

**2026-09-21, same session**

### Course correction before this phase: WGS split-by-gene walked back

Between Phase 2 and Phase 3, revisited the plan's Phase 7 draft (`GENE_SUBSET_AND_SPLIT`), which
had WGS/WES input split into one file **per gene** via `bcftools view -R` — at
`Complete_gene_list_V5` scale (6,517 genes) that's 6,517 per-gene jobs, the same job-explosion
problem this whole migration exists to move away from. Corrected in the plan: the default is now
**one `bcftools view -R <gene-list BED>` subset**, producing a single filtered VCF, chunked by
record count (reusing `SPLICEAI_CHUNK`'s existing mechanism) — not by gene identity. Per-gene
splitting stays available only as an explicit opt-in for when a per-gene file set is genuinely the
wanted deliverable, not the pipeline's default WGS/WES shape. Phase 7 itself isn't built yet;
today's Phase 1-3 modules are unaffected either way since they just take whatever single VCF
they're given.

### What changed vs. Phase 1/2

Like Pangolin and SPiP, `main.sh` runs VEP and Branchpointer directly on the whole per-gene VCF
(`main.sh:437-446` for VEP, `main.sh:488-497` for Branchpointer) — no chunking. Two more
single-process modules:

- `modules/local/vep.nf` — `VEP_ANNOTATE`, ported from `src/hpc/annotate_vep_vars.sh`. **First
  containerized process in this migration.** The bash version already ran VEP inside `vep.sif` via
  a manual `singularity exec --bind ...` call; this process instead uses a native `container =
  params.vep_sif` directive, so the `script:` block is a plain `vep --fork ...` invocation with no
  singularity wrapper — Nextflow's own singularity integration injects it. Every plugin/custom-file
  flag (CADD, dbNSFP, AlphaMissense, gnomAD, MaxEntScan, UTRAnnotator, SpliceVault, SpliceVarDB,
  REVEL, Mastermind, LoFtool, pLI, VEP's own precomputed-lookup SpliceAI plugin, SpliceRegion) kept
  verbatim — non-goal, this is broad-tier plumbing, not a re-tuning of VEP's annotation set. Note:
  VEP's `--plugin SpliceAI` is a precomputed hg38 score lookup, unrelated to and not a duplicate of
  this pipeline's own live `SPLICEAI_ANNOTATE` (gene-restricted tier) — `main.sh` already ran both
  together. Two deliberate, documented behavior changes from the bash original, not verbatim ports:
  `--fork ${task.cpus}` (main.sh submits `-pe smp 4` but never exports `$THREADS`, so the bash
  version's `--fork ${THREADS:-1}` actually ran with `--fork 1` — a latent bug, not an intentional
  throttle; this port uses real parallelism matching the allocated slots), and indexing with the
  **container's own** `tabix` (htslib 1.9, confirmed present in `vep.sif` by running `singularity
  exec vep.sif which tabix` in this sandbox) rather than the host conda `tabix` the bash script
  calls after `singularity exec` returns — this script runs entirely *inside* the container, and a
  host binary reached only via the bind mount isn't guaranteed to load correctly there.
- `modules/local/branchpointer.nf` — `BRANCHPOINT_ANNOTATE`, ported from
  `src/hpc/annotate_branchpointer_vars.sh`. Broad-safe by design (lookup against a precomputed
  genome-wide BED, `labranchor_grch38_top.bed.gz`, not live model inference — established during
  this migration's earlier research). Not containerized, matching the bash source's own choice
  (plain `spliceai_env` python, no container).
- `workflows/broad_pass.nf` — new `BROAD_PASS_SUBWORKFLOW`, runs `VEP_ANNOTATE` and
  `BRANCHPOINT_ANNOTATE` in parallel off the same `meta_vcf_ch`, matching `main.sh`'s independent
  qsub jobs. Mirrors `GENE_RESTRICTED_SUBWORKFLOW`'s `spip_tier` gate in reverse: `SPIP_ANNOTATE`
  only fires here if `params.spip_tier == 'broad'` (not the default, so a no-op today).
- `main.nf` — now calls both `GENE_RESTRICTED_SUBWORKFLOW` and `BROAD_PASS_SUBWORKFLOW`, all five
  predictors running off one shared input channel.
- `nextflow.config` — added `vep_sif`/`vep_dir`/`vep_plugins_dir`/`vep_cache_version`,
  `branchpointer_python`/`branchpointer_script`/`labranchor_bed`, and — only inside the `standard`
  (cluster) profile — `singularity.enabled`/`autoMounts`/`runOptions`. The bind-mount set mirrors
  `annotate_vep_vars.sh:33`'s manual `singularity exec --bind ... -e` flags, including `--cleanenv`
  (bash's `-e`) — without it, host `PERL5LIB`/`PYTHONPATH` leak into the container and can break
  VEP plugin loading in confusing ways. `local_dev` leaves `singularity.enabled` at its default
  `false`, so `VEP_ANNOTATE`'s `container` directive is simply ignored there and `-stub-run`
  executes it directly on the host. Also renamed `outdir`'s default a second time
  (`gene_restricted_out` → `annotation_out`, since it's now not even gene-tier-specific) — chose a
  name that shouldn't need renaming again as more predictors land.
- Retroactive fix to Phase 2's `pangolin.nf`/`spip.nf`: neither ported the bash originals'
  `OMP_NUM_THREADS`/`MKL_NUM_THREADS`/`OPENBLAS_NUM_THREADS` exports (Pangolin: `${task.cpus}`,
  matching the bash `${THREADS:-4}` pattern; SPiP: unconditional `1`, matching the bash script's
  own pin — SPiP does its own internal parallelism via `-t`, so this prevents double-parallelizing,
  not a mistake to "fix away"). Missing on a cluster job means numpy/torch can oversubscribe
  threads past the SGE-allocated slot count — a real gap in the Phase 2 commit, not reproducible or
  catchable in this local sandbox, caught on review before Phase 3 shipped.

### Verification

**VEP**: no local execution path at all — even beyond the missing container pull, the full VEP
cache/plugin data tree under `/references` (many plugin files, not just the one FASTA SpliceAI
needed) isn't mirrored locally. `-stub-run` only; real output must be checked on the cluster.

**Branchpointer**: the one predictor ported so far with a genuine local real-run path — its only
dependency is `pysam` against the LaBranchoR BED (checked by reading
`src/python/annotate_branchpointer.py`), both present in the local sandbox. Ran the exact command
`BRANCHPOINT_ANNOTATE` invokes directly (outside Nextflow, same script/args) against
`test_data/test_run/_tmp/MYBPC3.vcf.gz` and diffed the result (`grep -v "^##"`, so VCF header lines
were not compared, only the column header + 22 data records) against the known-good
`test_data/test_run/annotation/MYBPC3.annBranchpoint.vcf.gz`: **all 22 records identical**. Of
those, only 1/22 actually carried a Branchpointer/LaBranchoR annotation (the rest have no variant
within range of a branchpoint) — real, but thin, evidence; a gene with more annotated hits would be
a stronger check. This is the first real (non-orchestration-only) output-parity check in this
migration, not just a stub/DAG check. Not re-run through Nextflow itself in this pass, since
`-stub-run` is an all-or-nothing switch for the whole pipeline and the other four processes have no
local real-run path — the standalone check validates the same script/argument construction the
Nextflow process uses, which is what actually differs from the bash original.

Full `nextflow run main.nf -profile local_dev -stub-run` (both subworkflows together): exit 0, all
five predictors (`SPLICEAI_SPIKE`, `PANGOLIN_ANNOTATE`, `SPIP_ANNOTATE` (restricted), `VEP_ANNOTATE`,
`BRANCHPOINT_ANNOTATE`) completed `1 of 1 ✔` in parallel; `SPiP (broad)` correctly emitted nothing
(`spip_tier` still `'restricted'`). No new DAG/channel bugs surfaced — both new processes and the
new subworkflow reused Phase 1/2's established shape directly.

## What's next (not done in Phase 3)

Per the plan's migration sequencing: Phase 4 (merge/rejoin subworkflow, still the trivial 1:1
single-gene case) is next. Real model-output parity for VEP, and Nextflow-orchestrated (not just
standalone-script) parity for Branchpointer, still need to be checked on the actual cluster.

## Phase 4 — merge/rejoin subworkflow

**2026-09-21, same session**

### What it does

Combines all five predictors' output into one final annotated VCF per gene, ported from
`src/hpc/merge_vep_spip.sh`. VEP is the authoritative record set; SPiP, Pangolin, SpliceAI, and
Branchpointer are layered onto it in that order via `bcftools annotate -a <file> -c <tag>`, each
transferring exactly the INFO tag(s) that predictor writes. Confirmed the exact tag names by
reading the known-good reference headers under `test_data/test_run/annotation/`: SPiP writes one
INFO tag literally named `SPiP`, Pangolin one named `Pangolin`, SpliceAI one named `SpliceAI`, and
Branchpointer five explicit fields (`Branchpointer_prob`, `Branchpointer_U2_energy`,
`Branchpoint_disrupted`, `LaBranchoR_score`, `LaBranchoR_acc_dist`) — all match the bash original's
`-c` arguments exactly, and the final `MYBPC3.annotated.vcf.gz` reference header confirms this
order (CSQ from VEP, then SPiP, then Pangolin, then SpliceAI, then the five Branchpoint fields).

Still the trivial single-gene case (Phase 4 scope per the plan) — the plan's
`AMBIGUOUS_GENE_SOURCE` overlapping-gene dedup logic isn't needed yet, since no cross-gene join can
happen until Phase 7 introduces multi-gene "regions."

### New files

- `modules/local/merge.nf` — `MERGE_ANNOTATIONS`, the single process running the four sequential
  `bcftools annotate` calls.
- `workflows/merge_subworkflow.nf` — new `MERGE_SUBWORKFLOW`, joins the five predictor output
  channels (`vep`, `branchpoint`, `pangolin`, `spliceai`, `spip`) by `meta` via `.join()` — same
  operator Phase 1's `spliceai_spike.nf` already used for meta-keyed joins — into the single tuple
  `MERGE_ANNOTATIONS` needs.
- `main.nf` — wires it in. Since `spip_tier` gates SPiP into exactly one of
  `GENE_RESTRICTED_SUBWORKFLOW.out.spip` / `BROAD_PASS_SUBWORKFLOW.out.spip` (the other emits
  nothing, see Phase 2/3), `.mix()`es the two into one real stream before handing it to
  `MERGE_SUBWORKFLOW` — the merge doesn't need to know which tier was active.

### Two deliberate departures from the bash original (documented, not silent)

1. **`large_sv_vcf` handling dropped entirely.** Confirmed dead code (`BUG_TRACKER.md` `SV-1`): no
   script anywhere in this repo ever produces a `*.large_svs.vcf.gz` file (`grep -rn
   "large_sv" src/` outside `merge_vep_spip.sh` itself returns nothing), so the bash branch that
   concatenates it back in has always been a no-op in this pipeline as it stands today. Reviving
   SV/CNV handling is `PLT-102` (future, separate) — not something to silently recreate as an
   always-skipped branch here.
2. **No per-predictor soft-continue**, unlike the bash original's `if bcftools annotate ...; then
   ... else echo Warning ...; fi` around each step, which lets the merge "succeed" with a predictor
   silently missing from the output. That, combined with `ORCH-3` (the bash script's
   `CMD_EXIT_CODE=$?` is captured *after* its cleanup `rm` loop, not after the actual merge), means
   a real merge failure today can go completely unnoticed downstream. This port drops the
   soft-continue — any `bcftools annotate` failure fails the process. This is a channel/logic
   decision, not a shell-flag one: confirmed by inspecting a real `.command.sh` in this sandbox
   that Nextflow already runs every script under `#!/bin/bash -ue` by default, so a
   missing/malformed predictor file failing the task is Nextflow's native behavior, nothing this
   script opts into. `ORCH-3` is closed because exit status is now determined natively, not from a
   manually captured `$?` after unrelated cleanup — the failure-propagation fix the plan names as a
   byproduct of the migration (see plan Context, `ORCH-1`/`ORCH-2`).

### Verification

**Strongest verification in this migration so far.** Ran the exact `bcftools annotate` command
sequence `MERGE_ANNOTATIONS` executes directly (outside Nextflow) against the five known-good
`MYBPC3` reference predictor outputs already present under `test_data/test_run/annotation/`
(`MYBPC3.annVEP.vcf.gz`, `.annSPiP`, `.annPangolin`, `.annSpliceAI`, `.annBranchpoint.vcf.gz`), and
diffed the result against that same directory's `MYBPC3.annotated.vcf.gz`: **all 22 data records
identical** (header lines not compared). This confirms the merge order, `-c` tag lists, and
`bcftools` invocation are correct — not just that the channel wiring is right. Content-parity
against known-good production output has now been confirmed for `SPLICEAI_CHUNK`/`SPLICEAI_CONCAT`
(Phase 1, record-count parity), `BRANCHPOINT_ANNOTATE` (Phase 3, full content parity), and now
`MERGE_ANNOTATIONS` — three of seven ported process types. `SPLICEAI_ANNOTATE`, `PANGOLIN_ANNOTATE`,
`SPIP_ANNOTATE`, and `VEP_ANNOTATE` still need their real model-output parity checked on the
cluster; those four are the ones with an actual model/cache dependency this sandbox can't run.

`bcftools`/`tabix` for this check came from the local miniforge3 `genomics` env mirror (same
binaries `local_dev`'s `params.bcftools`/`params.tabix` already point at) — no new local
dependency needed.

Full `nextflow run main.nf -profile local_dev -stub-run` (all three subworkflows together): exit 0,
all five predictors plus `MERGE_ANNOTATIONS` completed `1 of 1 ✔`; the 5-way `.join()` correctly
produced exactly one merge task per gene (no duplication, no hang) with the `spip_tier`-gated
`.mix()` resolving cleanly. No new DAG/channel bugs surfaced.

## What's next (not done in Phase 4)

Per the plan's migration sequencing: Phase 5 (`VCF_TO_TABLE` + reports — mostly plumbing, already
gene-agnostic Python/R) is next. Phase 6 (full 217-gene panel parity run vs. real `bash main.sh`
output) is the gate before any WGS work (Phase 7) can start.

## Phase 5 — VCF_TO_TABLE (reports deliberately out of scope)

**2026-09-22, continuation of the same session**

### Scope correction mid-phase

The plan's Phase 5 description bundles `VCF_TO_TABLE` with "reports." Mid-phase, the user
explicitly steered: "Do not focus too much on reports as this will be done by module 3." This
matches existing memory (`project_acmg_tiering_deprioritized`: real classification now happens in
`clinical_variant_prioritization`, not here). Scope narrowed to `VCF_TO_TABLE` only — the four
downstream scripts (`filter_and_summarize.py`, `plot_annotation_results.R`,
`generate_interactive_report.py`, `generate_clinical_prioritization_report.py`) are **not ported**
this phase, deliberately, not forgotten.

### What it does

`modules/local/vcf_to_table.nf` — `VCF_TO_TABLE`, ported from `src/hpc/vcf2parsed.sh`. Two stages
in one process (matches the bash original, which never persists its intermediate TSV):
`vcf_parser_pysam.py` flattens the merged VCF's INFO/CSQ fields into a scratch TSV, filtered to one
transcript (`--filter_by Feature --gene_set <ENST>`); `filter_variants.py` then applies status
contracts/splicing metrics and writes the final table. This is the actual deliverable this module
produces for Module 3 to consume.

`meta.transcript` is now resolved for real in `main.nf` (previously always `null`, unused before
this phase) — ported from `main.sh:378`'s exact lookup: whole-word match against
`resources/gene_transcript_mapping.txt`'s `Gen,NM,ENST` columns, 3rd field. Deliberately the
read-only lookup only, not `query_new_transcripts.py --auto-append`'s live write path (`ORCH-9` —
unlocked file write, non-MANE API fallback — stays an explicit non-goal). One deliberate departure
from a verbatim port, caught in review: `main.sh` defaults an unresolved gene to the literal string
`'UNKNOWN'` and carries on; here, `'UNKNOWN'` reaching `VCF_TO_TABLE`'s `--gene_set` matches no VEP
CSQ block, so every variant gets skipped and `filter_variants.py` writes an empty-but-valid table —
exit 0, published, no warning. `main.nf` now fails the run instead when the lookup can't find a
transcript, alongside the existing `spip_tier` guard.

**Not ported**: the bash script's `/data_tmp` fast-scratch staging for its intermediate TSV
(PID/timestamp-unique filename + trap-based cleanup). That existed to avoid writing a large
intermediate file to *shared, multi-job-concurrent* scratch storage — a concern that doesn't apply
inside a Nextflow task's own isolated work directory. Noted as a candidate for Nextflow's native
`scratch` directive if cluster profiling ever shows it's needed, not silently dropped.

### Verification

Same pattern as Phase 3/4: can't run end-to-end through Nextflow in this sandbox (needs a genuinely
VEP-annotated VCF; every upstream predictor here is `-stub-run`-only), but the exact command
sequence was verified for real outside Nextflow. Ran `vcf_parser_pysam.py` + `filter_variants.py`
(local `vcf_parser`/`datasci` env mirrors, both confirmed present) against the known-good
`test_data/test_run/annotation/MYBPC3.annotated.vcf.gz`, using the real resolved transcript
(`ENST00000545968`, from `resources/gene_transcript_mapping.txt`), and diffed the result against
`test_data/test_run/results/MYBPC3.parsed.clean.pq`:

- Row count: **17/17, exact match**.
- Of 559 columns common to both outputs, **557 identical**. The only two that differed —
  `PRIORITY_TIER` and `VARIANT_PRIORITY_SCORE` — are ACMG/scoring-tier columns, exactly the domain
  already deprioritized to Module 3 (memory: `project_acmg_tiering_deprioritized`). Cause of the
  difference not investigated further per the scope steer — could be scoring-logic drift since the
  reference was generated, or a different original invocation; either way it's Module 3's surface,
  not a finding about this port's correctness.
- The new output also carries 17 columns the reference doesn't (`ACMG_CRITERIA`,
  `DISEASE_PHENOTYPE_MATCH`, `HPO_MATCH`, `TIER_RATIONALE`, etc.) — same domain, not investigated
  further per the user's scope steer.

(First diff attempt used `x is not True` to check a numpy `.all()` result, which is always `True`
in Python regardless of the boolean's actual value since numpy bools are never identity-equal to
the `True` singleton — flagged all 559 columns as "differing" falsely. Caught and fixed before
trusting the result; noting it since it's an easy mistake to repeat.)

Full `nextflow run main.nf -profile local_dev -stub-run`: exit 0, all five predictors plus
`MERGE_ANNOTATIONS` plus `VCF_TO_TABLE` completed `1 of 1 ✔`. No new DAG/channel bugs surfaced.

## What's next (not done in Phase 5)

Per the plan's migration sequencing: Phase 6 (full 217-gene panel parity run vs. real `bash
main.sh` output) is the gate before any WGS work (Phase 7) can start. The four report/filter
scripts skipped this phase remain unported — revisit only if Module 3 doesn't end up owning that
surface after all.
