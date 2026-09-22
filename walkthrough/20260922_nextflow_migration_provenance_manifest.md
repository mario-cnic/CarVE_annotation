# Nextflow migration — run-id fix + provenance manifest (item 12)

**2026-09-22, continuation of the same session**

## Context

Follows `walkthrough/20260922_nextflow_migration_whole_vcf_pipeline.md` (items 7-10, the
whole-VCF pipeline itself). This entry covers the two limitations flagged when that work landed
(`BUG_TRACKER.md` `MISC-7`, and `TODO.md`'s long-standing Priority-0 provenance gap) — plan item
12: provenance manifest. Dual-mode config (the other half of item 12) was already built in
Phase 1 (`nextflow.config`'s `local_dev`/`standard` profiles), so this item narrowed to the
manifest alone.

## What was built

### `run_id` fix (closes `MISC-7`)

`annotate_vcf.nf` used `workflow.sessionId` (Nextflow's own internal session UUID, fresh and
unrelated to the input on every invocation) as `meta.partition_id`. Changed to derive `run_id`
from the input filename by default — the exact same `.baseName.replaceAll(/\.vcf$/, '')` pattern
`main.nf` already uses for `gene_name` — with a new `--run_id` override for inputs whose filename
isn't a meaningful identifier on its own. Verified: `--input_vcf MYBPC3.vcf.gz` with no override
lands output under `nf_work/annotation_out/MYBPC3/`; `--run_id patient042_wgs` overrides it to
`nf_work/annotation_out/patient042_wgs/`.

### Provenance manifest (closes `TODO.md`'s Priority-0 item, for the Nextflow path)

- `src/python/write_run_manifest.py` — same schema `sarek_pipeline/src/write_run_manifest.py`
  already committed to for its own `PLT-012` (a `launch`/`completed` two-stage JSON manifest); a
  same-schema sibling per `DEC-0002` (repos stay separate, no shared import), not a copy-paste of
  sarek's file — adapted for this repo's own resolved paths (no `config/env.sh` here, everything
  comes from `nextflow.config`), and captures a few things sarek's copy doesn't need: SpliceAI/
  Pangolin's `-d`/`--distance` value actually used (the specific gap `TODO.md` named — "the
  SpliceAI `-D` question... there was nothing to check"), the VEP container's own file hash (VEP
  is containerized here, not a local binary, so `vep --version` isn't a meaningful query — the
  `.sif`'s hash is the real provenance artifact), and hashes of the two curated resource files
  (`gene_transcript_mapping.txt`, `v5_genes_loc.bed`) so a manifest can prove which version of
  each was in effect for a given run.
- `run_annotate_vcf.sh` — provenance-wrapped launch, mirroring `sarek_pipeline/run_sarek.sh`'s
  launch/run/completed pattern: writes the `launch` record, runs `nextflow run annotate_vcf.nf`
  with every argument passed through unmodified, writes the `completed` record, then rescues
  `RUN_MANIFEST.json` past `nf_work/`'s default-deny `.gitignore` with `git add -f` (not a
  `.gitignore` negation rule — git cannot re-include a file whose *parent directory* is itself
  excluded, same documented limitation sarek's own copy already worked around the same way).

## Verification

Ran `bash run_annotate_vcf.sh -profile local_dev -stub-run --input_vcf
test_data/test_run/_tmp/MYBPC3.vcf.gz` end to end (the script itself needed `bash`, not `./`, per
this sandbox's now-familiar `noexec` NFS mount — same constraint as every other script here).
Inspected the resulting `RUN_MANIFEST.json`:
- `run_id` resolved to `MYBPC3` (confirming the filename-derived default works through the wrapper
  too, not just `annotate_vcf.nf` directly).
- Real, correct values captured: `bcftools 1.24` (the actual local-mirror binary used),
  `vep_cache_version: "111"`, `spliceai_distance`/`pangolin_distance: "10000"`, `nextflow version
  24.04.2.5914`, real SHA-256 hashes for the input VCF and both curated resource files.
- `vep_sif_sha256: null` — correct, not a bug: no local mirror of `vep.sif` exists in this
  sandbox (`local_dev` profile), so there's nothing to hash; a cluster run would get a real value.
- `annotation_pipeline_repo.dirty: true` with the exact list of files this session had modified
  at the time — proves the dirty-tree detection works, and is itself a small demonstration of why
  this matters (this manifest is honest that its own run didn't match a clean committed state).
- `completed.workflow_stats` correctly parsed from `.nextflow.log`
  (`succeededCount=16; failedCount=0; ...`).
- `git status` confirmed `git add -f` successfully staged the manifest despite `nf_work/` being
  gitignored.

The verification run's own manifest was not committed (it's local test output, not a real run
worth preserving in history) — `nf_work/` was cleaned up after confirming the mechanism works.

## What's next

Plan item 12 is now complete. Item 11 (full-run verification against real multi-gene input and
`bash main.sh`'s output) remains the real parity gate, and remains cluster-only — this sandbox
cannot do real model inference at any scale, which is what that item actually needs.
