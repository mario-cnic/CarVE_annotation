# Nextflow migration — multi-gene test data + chunk-reassembly-order bug (item 11 prep)

**2026-09-22, continuation of the same session**

## Context

Preparing plan item 11 (full-run verification against real multi-gene input and `bash main.sh`'s
output — cluster-only, since it needs real model inference at real scale). Two corrections landed
first:

- **`RUNS/run_more_genes_20260817_1143` is real production data, not test data.** Initially
  proposed reusing its real per-gene inputs (which happen to include `EMD`/`LEMD2`/`TNNI3`/
  `TNNI3K`) as a "test" case — the user corrected this immediately. Not used for anything.
- **`ORCH-5` is not a genomic-coordinate-overlap issue.** Corrected in the plan file directly (see
  its item 11 for the full explanation): `ORCH-5` is `main.sh:390`'s unanchored
  `grep -q "$gene_name" "$GENES_BED_FILE"` falsely substring-matching `EMD` inside an existing
  `LEMD2` BED line (and `TNNI3` inside `TNNI3K`) — a bash string-matching bug, confirmed unrelated
  to physical position (`EMD`/chrX vs. `LEMD2`/chr6; `TNNI3`/chr19 vs. `TNNI3K`/chr1 — different
  chromosomes). Doesn't exist in the Nextflow pipeline at all (no per-gene "already in BED, or
  fetch" check anywhere in it) and stays open for `main.sh`, left for later per the user.

## What was built: `test_data/raw_vcfs/panel7_test.vcf.gz`

Purely synthetic, same style as the pre-existing `test_data/raw_vcfs/MYBPC3_test.vcf.gz` (confirmed
by diff: that file is literally the pre-`variant_converter` source that `test_data/test_run/_tmp/
MYBPC3.vcf.gz` — used for every single-gene check so far this session — was converted from).
7 genes, 72 synthetic variants, real GRCh38 gene-body coordinates (from `resources/
v5_genes_loc.bed`, minus its 10kb margin): `MYBPC3` (reused as-is), `ACTC1`, `TNNT2`, `PKP2`,
`DSP`, `SCN5A` (all in the 217-panel — tier-1 priority-tagging coverage), and `KCNQ1` (deliberately
*not* in the 217-panel — tier-2/3 fallback coverage). Spans 6 different chromosomes. Fake rsIDs/
sample-free, matching the existing file's synthetic convention — no real patient data anywhere.
Not tracked in git (`test_data/` is gitignored entirely, same as the pre-existing file already
was — confirmed via `git ls-files`, corrected an earlier wrong assumption that it was tracked).

Built with a scratch Python script (record generation) + `bcftools concat`/`sort`/`tabix` (file
assembly) — not committed, this was a one-off generation step; the resulting `panel7_test.vcf.gz`
is the artifact that matters.

## Real bug found: `CONCAT_CHUNKS` chunk-reassembly order (`BUG_TRACKER.md` `MISC-8`)

Running this multi-gene, multi-chromosome input through `annotate_vcf.nf` with a small
`--chunk_size 15` (forcing 5 chunks instead of 1, to actually exercise multi-chunk reassembly)
failed immediately: `bcftools concat` refused with `"The chromosome block 11 is not contiguous,
consider running with -a"`.

**Root cause**: the chunk-reassembly sort (`chunk_files.sort { it.name }` in each tier workflow's
`collectForConcat`, added in the previous phase specifically to counter `groupTuple()`'s
undocumented-order behavior) only works when each chunk's annotated output has a unique filename.
`SPLICEAI_ANNOTATE` names its output after the input chunk (`${chunk.baseName}.annSpliceAI.vcf.gz`
— genuinely unique per chunk). `VEP_ANNOTATE`/`BRANCHPOINT_ANNOTATE`/`PANGOLIN_ANNOTATE`/
`SPIP_ANNOTATE` all name theirs `${meta.partition_id}.<suffix>.vcf.gz` — identical across every
chunk of the same run (this is what the earlier `CONCAT_CHUNKS` `stageAs` fix, `MISC` from the
previous session's advisor review, was already working around at the *input-staging* level). So
for these four predictors, `chunk_files.sort { it.name }` compares identical strings — a no-op —
and chunks reassemble in whatever order `groupTuple()` happened to collect them (task completion
order, not submission order). With single-chunk test data (every check before this one), there was
only ever one element to "sort," so this never mattered. With 5 real chunks and 6 chromosomes,
`chr11`'s 32 records ended up split across non-adjacent chunks, and `bcftools concat` correctly
refused to silently produce a malformed file.

**Fix** (`modules/local/chunk.nf`, `CONCAT_CHUNKS`): rather than solve chunk-order preservation
for four processes' worth of channel wiring, made the *reassembly* robust to arbitrary input order
instead. Each chunk is now tabix-indexed individually right before concatenation, `bcftools concat
-a` (order-tolerant, but requires those per-input indices — confirmed directly: without them it
fails with "Could not retrieve index file") reassembles them, and `bcftools sort` runs on the
result regardless, as a final correctness guarantee independent of whatever order concat received
its inputs in. The existing record-count parity check is unaffected — it still catches any real
record loss, it just no longer depends on getting reassembly order right first.

## Verification

Ran `nextflow run annotate_vcf.nf -profile local_dev -stub-run --input_vcf
test_data/raw_vcfs/panel7_test.vcf.gz --chunk_size 15` (forcing 5 chunks) **three independent
times** to confirm the fix is stable, not lucky ordering: all three exit 0, all five predictors'
`5 of 5 ✔` chunk fan-out plus their `CONCAT_*` reassembly succeeded every time. Merged output
(`panel7_test.annotated.vcf.gz`) has exactly 72 records, and the per-gene breakdown matches the
constructed input exactly (`ACTC1: 6, DSP: 9, KCNQ1: 10, MYBPC3: 22, PKP2: 10, SCN5A: 6, TNNT2: 9`).
This is also the first time this migration has tested: multiple genes in one whole-VCF run, a gene
outside the 217-panel flowing through `GENE_SUBSET`/tiering correctly, and genuine multi-chunk
reassembly for all five predictors (not just SpliceAI, which Phase 1 already covered).

## What's next

`panel7_test.vcf.gz` is ready as the input for a real cluster run of `annotate_vcf.nf`
(`bash run_annotate_vcf.sh -profile standard --input_vcf test_data/raw_vcfs/panel7_test.vcf.gz`)
— the actual execution needs real SGE job submission, which this sandbox cannot do
(`CLAUDE.md` bars `qsub`/`qstat`/`qdel` here); this is the user's own cluster session to run.
