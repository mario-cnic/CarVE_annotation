# Predictor score cache for SpliceAI and Pangolin

Date: 2026-10-22T17:23:46+02:00. Code: `src/python/annotation_cache.py`, `modules/local/annotation_cache.nf`,
`workflows/gene_restricted_whole.nf`, `run_annotate_vcf.sh`. Tests: `tests/python/test_annotation_cache.py`.

## Why

SpliceAI and Pangolin are 89% of a WGS family run's task-hours and recompute variants already scored
in earlier runs. The cache scores each variant once per predictor version and looks it up afterwards.
Scope is these two predictors only; it is not PLT-040's store of all site annotations.

## Evidence (measured 2026-10-22)

FACT unless labelled. Commands were run on the published outputs; temporary files in the session scratchpad.

| | Value |
|---|---|
| Task-hours, F1830_union (`reports/perf_review_20261019_F1830/per_process.tsv`) | Pangolin 809, SpliceAI 570, VEP 86, SPiP 72, Branchpoint 0.4; total 1,546 |
| Model cost per gene-region variant (4 slots) | Pangolin ~3.0 s, SpliceAI ~2.1 s (task-hours ÷ 975,496 sites) |
| Lookup of 975,496 sites against a 984k-site cache (`bcftools annotate -a`) | 6 s, 32 MB |
| F1830_union_avi sites already in F2005_union (CHROM:POS:REF:ALT) | 608,244 / 975,496 = 62% |
| Same site in both runs, predictor INFO compared as text | Pangolin 608,244/608,244 identical; SpliceAI 608,244/608,244 identical |
| SPiP, same comparison | 49,933 identical as text; all 608,244 equal after removing padding and comparing numbers (output formatting depends on the chunk) |
| Shared sites computed without a score | Pangolin 61,101; SpliceAI 170,090 (the cache must keep these) |
| Cache size | ~22 MB per million sites (Pangolin, sites-only) |

INFERENCE: at a 62% hit rate, a family like F1830 needs ~850 fewer task-hours, and the hit rate grows
as families accumulate. Wall time should also fall (fewer Pangolin chunks, less of the queue throttling
seen on F1830); to be measured on the next WGS run.

Not cached: VEP (cheap, and its resources change often), Branchpoint (~0 h), SPiP (5%, non-reproducible
text, RefSeq database rebuild pending). Pre-scoring every gnomAD variant in the gene BED (~1.2 Gb) is not
feasible at ~3 s per variant. The Illumina precomputed SpliceAI files already loaded by VEP use GENCODE
v24 and a short window, so they cannot replace `-D 4999`.

## Format

bgzip + tabix sites-only VCF, one append-only folder per (predictor, version). This is how precomputed
scores are distributed (Illumina SpliceAI, CADD, dbNSFP, AlphaMissense, AlphaGenome AVI), it is read
natively by bcftools and VEP, and append-only files are safe on the shared filesystem (SQLite with
concurrent writers on NFS is not). Other designs seen: on-demand caches keyed per variant (Broad
SpliceAI-lookup, NIH `spliceai-wrapper` with SQLite, `cached_spliceai` with a key-value store) and
echtvar's compressed archives (support for composite string fields not verified; not chosen).

## Design

```
<cache_dir>/<predictor>/<version_id>/VERSION.json
<cache_dir>/<predictor>/<version_id>/segments/<run_id>.<UTC stamp>.vcf.gz(.tbi) + .json
```

- A record per evaluated site, INFO limited to `SpliceAI` or `Pangolin`/`PangolinTissue`. A record without
  the field = computed, no score; no record = not in the cache.
- Minimal header: fileformat, contigs, the predictor's INFO lines, `##annotation_cache_*` metadata.
- The `.json` sidecar (size, sha256, record count, run, version) is written last. A segment whose file does
  not match its sidecar is skipped, so its sites are recomputed rather than misread.
- `version_id` = first 16 hex of the sha256 of the identity: FASTA size and sha256 or `.fai` sha256,
  `pangolin_db` sha256 (Pangolin), the conda env fingerprint, the distance, and the scoring code at the
  run's commit (git blob ids of `annotate_spliceai.py` or `Pangolin-main/pangolin/` without hidden files,
  plus the sha256 of the `SPLICEAI_ANNOTATE`/`PANGOLIN_ANNOTATE` `script:` block). Taken from
  `RUN_MANIFEST.json`; a missing field or uncommitted scoring code stops the launch.
- Pipeline (`--annotation_cache_dir` set): `CACHE_SPLIT` splits the gene-subset VCF into hits (INFO copied
  with `bcftools annotate`; records are biallelic, so its default allele matching is exact; the cluster's
  bcftools 1.21 has no `--pair-logic`) and misses; only misses are chunked and scored; the hits
  file (header-only when empty) joins the chunks at `CONCAT_CHUNKS`, whose record-count parity check covers
  both; `CACHE_WRITE` publishes the computed misses as one segment. Report per predictor:
  `pipeline_info/<run_id>.<predictor>.cache_split.json`.
- The launcher writes the segments a run reads to `pipeline_info/annotation_cache_segments.<predictor>.txt`
  and passes that list to Nextflow. A `-resume` reuses the list of the original launch, so the run's own
  new segment does not re-run `CACHE_SPLIT` and everything after it; a new launch with the same run_id
  lists the folder again and reads that segment (all hits).
- Default on (Mario, 2026-10-22): `annotation_cache_dir` = `/data_lab_PGP/resources/annotation_cache` under
  `standard`, `false` under `local_dev`; `--annotation_cache_dir false` turns it off. `CACHE_WRITE` is skipped
  under `-stub-run`, whose predictor outputs carry no scores.
- Not done: compaction of many segments into one (`build_lookup` merges them at each split; cost grows
  with the number of segments). The version ids are passed to Nextflow after the manifest is written,
  so the manifest's recorded command line does not show them.

Version ids of the existing runs (computed 2026-10-22): F1830_union_avi SpliceAI `ca2807064ffcf846`,
Pangolin `b21e77aef22ff4e8`. F2005_union: Pangolin `b21e77aef22ff4e8` (identical); SpliceAI
`57350be21714b929`, whose only difference is `annotate_spliceai.py` (an unused `import sys` removed between
`48233a3` and `fd79df6`).

## Checks (2026-10-22)

- `tests/python/test_annotation_cache.py`: 19 passed (version identity, split incl. unscored sites and
  two alleles at one position, damaged segment, segment header, segment list, verify conflicts, seed).
- Full suite: only `test_provenance_stamp` and `test_tag_transcript_priority` fail, as at `fd79df6`.
- Local `-stub-run` (`local_dev`, `clinvar_af_fixture.vcf.gz`, cache in the scratchpad): run 1 with an empty
  cache 5/5 misses and one segment per predictor; run 2 5/5 hits, no misses chunking, no SpliceAI/Pangolin
  task, no cache write, CONCAT/MERGE ran, predictor VCF records identical to run 1; run without
  `--annotation_cache_dir` unchanged; `--annotation_cache_dir` without version ids stops with an error.
  With the frozen segment list: run, then `-resume` (every task cached, MERGE and table included), then a
  new launch of the same run_id (5/5 hits, no new segment).
- Unit tests also pass with the cluster's bcftools 1.21 (`/data_lab_PGP/shared/utils/conda_envs/genomics`).
- Seeding dry run on the real outputs into the scratchpad (not the shared cache): F1830_union_avi 975,496
  sites per predictor (~8 s each); F2005_union adds 375,760 sites with 608,244 shared sites and 0 conflicts
  (Pangolin directly; SpliceAI refused without `--override-reason`, which names the one differing file).
  Two-family cache: 17 MB (Pangolin), 19 MB (SpliceAI); `verify` 0 conflicts. `split` of F1830's 975k
  sites (3 samples) against it: 65 s, 32 MB, all hits.
- Not run: the cluster. Stub runs prove the wiring, not real scores.

## Cluster check (2026-10-22)

Test cache `/mnt/scratch/mruizp/annotation_cache_test`, input `test_data/raw_vcfs/clinvar_af_fixture.vcf.gz`.

- `cache_fixture_1`, launched 17:37 from `fd79df6` with the cache code uncommitted (`RUN_DIRTY` patch =
  the content committed as `038b679`/`040557a`): exit 0; SpliceAI and Pangolin 5/5 misses; one segment each.
- A first `cache_fixture_2` launched while run 1 was still running read an empty segment list and was
  stopped. Relaunched 17:43 from `040557a` (clean) after run 1 finished: exit 0 in under 3 min; 5/5 hits for
  both predictors; no `CHUNK_*_MISSES`, `SPLICEAI_ANNOTATE`, `PANGOLIN_ANNOTATE` or `CACHE_WRITE` task;
  `CONCAT_*` and `MERGE_ANNOTATIONS` ran; SpliceAI, Pangolin and PangolinTissue values identical to run 1.
- Final tables: 199 rows x 586 columns in both; all values equal except `5UTR_annotation`, which holds the
  same key=value pairs in a different order (VEP UTRannotator output order varies between runs; VEP is not
  cached).
- Runs launched together do not see each other's segments: each reads the list frozen at its launch.

## Seeding (2026-10-22)

`/data_lab_PGP/resources/annotation_cache/` created by `annotation_cache.py seed`: SpliceAI
`ca2807064ffcf846` and Pangolin `b21e77aef22ff4e8`, each with F1830_union_avi (975,496 sites) and F2005_union
(375,760 new sites; 608,244 shared with F1830, 0 conflicts). F2005 SpliceAI entered with
`--override-reason` (recorded in its segment sidecar). `verify`: 0 conflicts. 36 MB in total.

## Pending

1. Next WGS family run with the cache: hit rate, task-hours, wall time vs F1830 (23.7 h). A re-run of
   F1830 or F2005 will be ~100% hits and does not measure a new family's hit rate.
2. Segment compaction once segments accumulate.
