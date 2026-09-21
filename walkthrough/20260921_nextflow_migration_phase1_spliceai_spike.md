# Nextflow migration — Phase 1 SpliceAI spike

**2026-09-21**

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

## What's next (not done here)

Per the plan's migration sequencing: Phase 2 (Pangolin + SPiP at `restricted` tier default) is the
next slice, still validated on the legacy per-gene adapter before any WGS-entry-point work starts
(Phase 7). Real model-score parity for this Phase 1 slice should be checked the next time this
runs on the actual cluster. `carve-platform`'s D-3 Verdict and a Module 3 decision-record note (per
this repo's own CLAUDE.md rule) are flagged in the plan as follow-up, not done in this session.
