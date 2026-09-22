# Nextflow migration — whole-VCF pipeline (items 7-10)

**2026-09-22, continuation of the same session**

## Context

Follows the course-correction and `resources/v5_genes_loc.bed` build recorded in
`walkthrough/20260922_nextflow_migration_whole_vcf_redesign.md`. This entry covers building the
actual new primary entry point (`annotate_vcf.nf`) — plan items 7-10: the whole-VCF entry point,
generalized chunking reused by both tiers, the tier subworkflows rewired to fan out/reassemble per
chunk, and `VCF_TO_TABLE` reworked with the new transcript-priority tagging step.

## What was built

- `modules/local/chunk.nf` — `CHUNK_VCF` (generalizes Phase 1's `SPLICEAI_CHUNK`, no tool-specific
  naming) and `CONCAT_CHUNKS` (generalizes `SPLICEAI_CONCAT` with an explicit `suffix` param so one
  process serves all five predictors).
- `modules/local/gene_subset.nf` — `GENE_SUBSET`, one `bcftools view -R
  resources/v5_genes_loc.bed` per run (not per gene — the correction from the previous entry).
- `workflows/broad_pass_whole.nf` / `workflows/gene_restricted_whole.nf` — the two tiers, rewired
  to chunk → fan out the *existing, unmodified* predictor processes (`VEP_ANNOTATE`,
  `BRANCHPOINT_ANNOTATE`, `PANGOLIN_ANNOTATE`, `SPIP_ANNOTATE`, `SPLICEAI_ANNOTATE`) per chunk →
  reassemble via `CONCAT_CHUNKS`.
- `annotate_vcf.nf` — the new top-level entry point. `main.nf` (legacy per-gene adapter) is
  untouched, exactly as the plan specified.
- `modules/local/{vep,branchpointer,pangolin,spip}.nf` — one line changed each: `publishDir`'s
  plain form replaced with a `saveAs` closure that skips publishing when
  `meta.partition_type == 'chunk'` (whole-VCF mode). Reused as-is otherwise.
- `modules/local/vcf_to_table.nf` — `--gene_set`/`--filter_by` now conditional on `meta.transcript`
  being present; `publishDir` also skips whole-VCF-mode runs (superseded by
  `TAG_TRANSCRIPT_PRIORITY` below).
- `src/python/tag_transcript_priority.py` + `modules/local/tag_transcript_priority.nf` — the new
  transcript-priority tagging step (tier 1: 217-panel curated, tier 2: VEP's own `MANE_SELECT`,
  tier 3: everything else, kept not dropped). Only ever called from the whole-VCF path.

## Two real bugs found and fixed while building this

### 1. Contig-naming mismatch — `v5_genes_loc.bed` used `chr11`, pipeline VCFs use `11`

`GENE_SUBSET`'s `bcftools view -R` against a fresh `v5_genes_loc.bed` silently returned **zero
records** for a real MYBPC3 test VCF, even though MYBPC3 is unambiguously in `Complete_gene_list_V5`
and had a correct-looking BED entry. Root cause: `resources/v5_genes_loc.bed` was built from the
MANE GTF, which uses `chr`-prefixed contigs (`chr11`) — but every VCF actually flowing through this
pipeline uses bare chromosome names (`11`), confirmed by checking both a pre-annotation test file
and a real production VEP output (`test_data/test_run/annotation/MYBPC3.annVEP.vcf.gz`). The
reference FASTA itself *does* use `chr`-prefixed contigs; VEP bridges that gap internally via
`--synonyms chr_synonyms.txt`, but `bcftools view -R` requires an exact string match and has no such
bridging. `bcftools view -R` doesn't error on a contig-naming mismatch, it just returns nothing —
exactly the kind of silent-empty-tier failure this session has caught twice before (`spip_tier`,
`VCF_TO_TABLE`'s `UNKNOWN` transcript guard). Fixed in
`src/python/build_gene_restriction_bed.py`: strip the `chr` prefix before writing the BED.
`resources/v5_genes_loc.bed` regenerated with the fix — same 5,742/5,877 gene coverage as before,
now genuinely usable.

### 2. `CONCAT_CHUNKS` input/output filename collision

Took most of this session's remaining time to isolate — logged in detail because the investigation
itself is worth remembering. Every whole-VCF run failed at every `CONCAT_*` step with `bcftools`
reporting `Input is not detected as bcf or vcf format` against a 0-byte file, `exit 255`, retried 5
times, always failing identically. Manually re-running the exact failing task's own `.command.run`
moments later always succeeded — which read as a filesystem timing race (this sandbox's home
directory is genuinely NFS-mounted, confirmed via `mount`), and several plausible-looking races were
ruled out empirically before the real cause surfaced:
- Moving `workDir` to local `/tmp` — no change.
- `-process.maxForks 1` — no change (this only limits concurrency *per process type*, not
  globally — a real gap in that test, caught only after re-reading the flag's actual semantics).
- `-qs 1` (true global serialization) — no change. This is what should have ruled out timing
  entirely, and in hindsight did — but the investigation kept going a bit past that point.
- A minimal isolated repro replicating the same chunk→annotate(stub)→concat topology, including
  two independent chunk/concat chains sharing a `meta.id`, the real `nextflow.config` resource
  labels, and the conditional `publishDir` — repeatedly succeeded, never reproducing the failure.

Called `advisor` after this — the actual cause was immediately visible from the failing command
line, which the investigation had been reading past: `bcftools concat X.annVEP.vcf.gz -Oz -o
X.annVEP.vcf.gz` — **input and output are the identical filename.** Every per-chunk `ANNOTATE`
process (reused unmodified from the legacy path) names its output
`${meta.partition_id}.annVEP.vcf.gz`; in whole-VCF mode `meta.partition_id` is the run ID, constant
across every chunk, and `CONCAT_CHUNKS` declared its own output under that exact same convention.
With a single chunk, Nextflow staged that one file into `CONCAT_CHUNKS`'s task directory under its
original name — the same name `bcftools concat -o` was about to write to — so `bcftools` truncated
its own only input to 0 bytes before ever reading it. This is the identical class of bug already
caught once this session, in `TAG_TRANSCRIPT_PRIORITY` (`modules/local/tag_transcript_priority.nf`,
via `stageAs`) — just not applied to `CONCAT_CHUNKS` too. Every "later, unmodified, it worked" replay
was simply re-running against the *unstaged, un-truncated* source file in the upstream task's own
directory, not the actual failure path.

Fixed with the same mechanism: `path(annotated_chunks, stageAs: 'chunk_??/*')`, staging each chunk
into its own numbered subdirectory so it can never collide with the declared output filename. That
introduced one more instance of the already-documented scalar-vs-List `Path` gotcha (Phase 1,
`java.nio.file.Path` implements `Iterable<Path>`): with exactly one chunk, `annotated_chunks` bound
to a single two-segment `Path` (`chunk_01/<file>`), and `.join(' ')` on a bare `Path` iterates its
path *segments*, producing two bogus arguments instead of one real path. Fixed with the same
`instanceof List` normalization Phase 1 already established.

A `sync chunks/` was added to `CHUNK_VCF` partway through this investigation, based on a race
hypothesis that turned out to be wrong — removed once the real cause was found, so a wrong
explanation doesn't sit in the codebase as if it were confirmed.

## Verification

Full `nextflow run annotate_vcf.nf -profile local_dev -stub-run --input_vcf
test_data/test_run/_tmp/MYBPC3.vcf.gz`: **exit 0**, all five predictors, `MERGE_ANNOTATIONS`,
`VCF_TO_TABLE`, and `TAG_TRANSCRIPT_PRIORITY` completed `1 of 1 ✔`. Published output under one run
directory (`nf_work/annotation_out/<run-id>/`): five predictor VCFs, one merged
`<run-id>.annotated.vcf.gz`, one `<run-id>.parsed.clean.pq` — no chunk-level clutter (confirming the
conditional `publishDir`/`saveAs` fix on the four reused predictor processes works correctly).

`MERGE_ANNOTATIONS` has no stub (runs its real `bcftools annotate` chain even under `-stub-run`,
same as Phase 4) — its own record-count parity held: 22/22, matching the MYBPC3 test input exactly,
confirming the merge step genuinely works correctly when fed whole-run (not one-gene) files, not
just that the DAG wires up.

`VCF_TO_TABLE`/`TAG_TRANSCRIPT_PRIORITY` themselves are stub-only in this full-pipeline run (their
real execution needs a genuinely VEP-annotated input, unavailable when every upstream predictor is
stubbed) — but their real logic was already verified standalone, outside Nextflow, against the
known-good MYBPC3 reference earlier this session (see the previous walkthrough entry and this
session's `tag_transcript_priority.py` verification: 238 unfiltered rows, tier distribution
`{1: 17, 2: 9, 3: 212}`, tier 1's count exactly matching the legacy path's single-transcript-filtered
row count).

Also explicitly checked for regression on the legacy path: the `publishDir`/`saveAs` change applied
to `vep.nf`/`branchpointer.nf`/`pangolin.nf`/`spip.nf`/`vcf_to_table.nf` controls *whether files
get published*, which a stub run's exit code alone says nothing about (both paths already exit 0
whether or not the change is correct). Ran `main.nf -profile local_dev -stub-run` on the same
MYBPC3 input and confirmed the same seven artifacts Phase 5 produced are still published under
`nf_work/annotation_out/MYBPC3/`: five predictor VCFs, `.annotated.vcf.gz`, `.parsed.clean.pq` —
`partition_type == 'gene'` still publishes correctly, no regression from the whole-VCF changes.

## Known limitations, not addressed this round

- **`CONCAT_CHUNKS`'s record-count parity check is per-tool, not end-to-end.** Each tool's
  reassembled output is checked against its own pre-chunk input (the subsetted VCF for the
  gene-restricted tier, the full VCF for the broad tier) — nothing compares the broad tier's count
  against the restricted tier's, and nothing separately re-checks `MERGE_ANNOTATIONS`'s output
  count (it doesn't assert one itself, matching Phase 4). "Parity confirmed" in the logs means
  chunk-reassembly parity per tool, not a claim about end-to-end record parity across the whole
  run.
- **`meta.partition_id` is a Nextflow session UUID** in whole-VCF mode, so output lands under
  `nf_work/annotation_out/<uuid>/` — functional, but every rerun gets a fresh, unrelated directory
  name with no link back to the input sample. Module 3 will be consuming these tables; a UUID
  directory is not a handoff contract anyone can plan around, and this needs a real naming scheme
  (sample ID, run label, or similar) before this path is used for anything beyond local
  verification — not solved by this round.

## What's next

Per the plan: item 11 (full-run verification against real multi-gene input and the equivalent
`bash main.sh` output) is the real parity gate — only meaningfully runnable on the actual cluster,
since it needs real model inference at real scale, which this sandbox cannot do. Item 12
(provenance manifest + dual-mode config) remains unstarted, parallelizable with anything else — the
run-directory-naming gap above is a natural fit for that work, since a real provenance manifest
would need a meaningful run identifier anyway.
