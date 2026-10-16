# Dead code and comment cleanup

**Date:** 2026-10-15T15:48:51+02:00 · **Base commit:** `6aa0607` · **Requested by:** Mario ("remove unnecessary, redundant or bloat from the code")

No behaviour change intended for `annotate_vcf.nf`. The legacy per-gene path was removed with Mario's approval.

**Commits:** `625c7e3` (untracked junk, six tests), `17ce662` (Python dead code), `f870a8a` (test imports), `92f3584` (wrapper comments), `89c71c8` (legacy per-gene path, `.nf`/config comments), then the docs commit.

## Done

- **Unused imports** removed from 16 files in `src/python/` and `src/tools/`, and from 12 retained tests. Found with an AST scan of tracked `.py` files; each hit was checked by hand. `from __future__ import annotations` is kept: it is a directive, not an unused name.
- **Unreferenced functions** removed: `audit_run_results.parse_time_from_log`, `io_qc.sort_by` (unused since the 2026-10-06 row-order change), `vcf_parser_pysam.add_sample_header` (body was `pass`). `git grep` and `git log -S` found no callers. Also removed: a commented-out import and the now-unused `typing.Literal` import.
- **Comments trimmed** to say briefly what the code does, per the TODO "Comments and code hygiene" item: `annotate_vcf.nf`, `nextflow.config`, `modules/local/{chunk,vep,vcf_to_table,merge,spip,branchpointer,tag_transcript_priority,spliceai}.nf`, `workflows/{broad_pass_whole,gene_restricted_whole,chunk_fan_concat,merge_subworkflow}.nf`, `run_annotate_vcf.sh`, and one block in `parsing/modules/splicing.py`. Removed: phase/plan history, `/home/mruizp/.claude/plans/...` paths, `main.sh:NNN` pointers, tracker IDs, and records of past verification. That history stays in git and in the earlier walkthroughs.

## Checks

| Check | Result |
| --- | --- |
| Every `"""` script block in the edited `.nf`/config files is byte-identical to `HEAD`; non-comment code lines are identical; no `//` appears inside a script block (the MISC-15 failure mode) | pass. The only code-line difference is trailing comments on `take:`/`emit:` lines in `chunk_fan_concat.nf` |
| `local_dev -stub-run` of `annotate_vcf.nf` on `test_data/raw_vcfs/panel8_trio.vcf.gz --genotype_wide true`, before and after the edits, isolated in a scratch work dir and outdir | exit 0 both times; same 20 tasks; same published file list |
| `bash -n run_annotate_vcf.sh`, `py_compile` on every edited `.py` | pass |
| `tests/run_all_tests.sh` (local, before and after) | 229 passed, 1 failed, 1 skipped after the edits. The failure, `test_tag_transcript_priority.py::test_streaming_parquet_equals_pandas_path`, also failed before any edit: the pandas metadata `attributes` key differs between library versions |
| `-resume` of the HEAD stub run with the edited code (same work dir) | all 20 tasks `Cached`, including `CHUNK_VCF`/`CONCAT_*`, which have no stub and hash their real `script:` |
| Toy one-process pipeline run for real (not stub), Nextflow 24.04.2: comments edited before `input:`, before `script:`, between `def`s inside `script:`, and after the closing `"""` | `Cached` on `-resume`. Positive control: a one-character edit inside the script block resubmits |

Not run: a cluster run. The cache checks used local Nextflow 24.04.2; the cluster runs 25.10.2, so a resume there is not tested.

## Removals (Mario's decisions)

The session's tool permission blocked file deletion on `/data_lab_PGP` at first, in line with the global rule against `rm` there; Mario ran items 1-2, and approved item 3 on 2026-10-16.

1. **Tracked junk.** `0` (144 bytes of `1` lines, a stray redirect): removed by Mario with `git rm 0` on 2026-10-15 (staged). `src/external/Pangolin-main/pangolin/.fuse_hidden0000252700000002` and `nf_work/annotation_out/S223/RUN_MANIFEST.json` (S223's real manifest): untracked by Mario with `git rm --cached` on 2026-10-15. Both files are still on disk, and both are now ignored (`.fuse_hidden*` and `nf_work/` in `.gitignore`).
2. **Tests that never call repo code.** Each imports a module and then asserts only pandas or dict behaviour: `test_filter_gene_transcript.py`, `test_refseq_to_ensembl.py`, `test_combine_valencia_hic.py`, `test_eda_parsed_pq.py`, `test_noonan_string_audit.py`, `test_parse_raw_hic.py`. The last one takes about 42 s. Removed by Mario with `git rm` on 2026-10-15. After the removal, `tests/run_all_tests.sh` gives 222 passed, 1 failed (the same pre-existing `test_streaming_parquet_equals_pandas_path`), 1 skipped.
3. **Legacy per-gene Nextflow path.** README called it "kept working but no longer developed". Removed (approved by Mario 2026-10-16T16:03+02:00): `main.nf`, `workflows/{broad_pass,gene_restricted,spliceai_spike}.nf`, `src/python/query_new_transcripts.py`, `tests/python/test_query_new_transcripts.py`; BUG_TRACKER ORCH-9 marked obsolete.
   - Also: `SPLICEAI_CHUNK`/`SPLICEAI_CONCAT` removed from `modules/local/spliceai.nf`; the `meta.transcript` (`--gene_set`/`--filter_by`) branch removed from `VCF_TO_TSV`; the `publishDir ... saveAs: partition_type == 'chunk' ? null` rules removed from `VEP_ANNOTATE`, `BRANCHPOINT_ANNOTATE`, `PANGOLIN_ANNOTATE`, `SPIP_ANNOTATE` and `TSV_TO_TABLE` (`annotate_vcf.nf` always sets `partition_type: 'chunk'`, so they never published); README line removed. `meta.partition_type` is kept, because `meta` is part of every task hash and removing the key would invalidate every cached task.
   - Check: `-resume` of the scratch stub run with these edits: exit 0, all 20 tasks `CACHED`, same published files. Under `-stub-run` the task hash uses the stub block, so this does not test the `VCF_TO_TSV` script change; on a real `-resume`, `VCF_TO_TSV` and its downstream tasks will run again once (INFERENCE).
   - After the deletions: no `.nf`/config/`.sh`/`.py` file references the removed files (the `main.nf` strings left in `test_provenance_stamp.py` and `test_write_run_manifest*.py` are fixture file names in temp dirs); `tests/run_all_tests.sh` gives 220 passed, 1 failed (same pre-existing failure), 1 skipped; a fresh isolated stub run exits 0 with 20 tasks `COMPLETED` and the same published files as the pre-edit run.
   - `--gene_set`/`--filter_by` stay in `vcf_parser_pysam.py` as CLI options. `resources/gene_transcript_mapping.txt` stays: other scripts read it.

Out of scope, left in place: `archive/deprecated_scripts/` (deliberately archived); one-off scripts in `src/tools/` and `reports/*/` (the provenance of documented analyses); the vendored Pangolin tree; `collectForConcat` duplicated in the two tier workflows (a refactor, not dead code).
