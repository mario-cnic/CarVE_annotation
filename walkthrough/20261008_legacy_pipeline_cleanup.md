# Legacy pipeline cleanup

Updated: 2026-10-08T11:57:54+02:00

## Scope

Mario authorized moving the reviewed old code/tests/container to trash, explicitly forbidding removal with `rm`. Twenty original paths (47 files, 2,003,382,695 bytes) were renamed within the same project filesystem into `trash/legacy_pipeline_20261008_115330/`. Every moved file was verified by size and inode after rename. [Manifest](../reports/legacy_cleanup_20261008/manifest.json) records original/current paths, metadata and pre-cleanup commit `107498d`.

Moved: `main.sh`, `config/`, `src/hpc/`, `src/R/`, `src/downstream/`, five legacy Python filter/report/dashboard scripts, their six test files, `tests/R/`, `tests/bash/`, the old readiness audit tied to shell wrappers, and `resources/containers/` (old image plus build definition). The old image alone was about 1.87 GiB; moving it does not reclaim disk space.

Kept: active Nextflow modules/workflows/launchers, vendored parser, predictor runners, active Pangolin code/weights, reference builders, existing raw/test data, references, outputs, logs, Nextflow work/cache and prior archives. `_tmp/pangolin/genes.gtf` remains pending provenance review. Module 3 and sarek tracked-source searches found no executable caller of the retired paths; one historical Module 3 plan mentioned `main.sh`. No sibling code was edited.

## Dependencies and documentation

Current Nextflow calls tools directly; references to `src/hpc/` in its source are historical comments. This corrects the old TODO assertion that those wrappers were shared with Nextflow.

`tests/python/test_population_frequency.py` retains active parser/scoring assertions; the retired standalone-filter/report dependencies were removed. A pre-edit snapshot is retained under the trash folder's `_snapshots/`. `test_vep_custom_columns.py` now checks the active Nextflow command instead of requiring the retired shell wrapper. The old manual trio dashboard script was archived; retained genotype/compound-het tests still cover their own logic.

`tests/run_all_tests.sh` runs retained Python modules independently to isolate their existing import paths. It no longer runs obsolete R/Bash tests. README points to active Nextflow usage and retires old launcher/container build instructions. TODO completion is archived; historical bug/provenance records have dated scope corrections rather than claims that old wrapper bugs were fixed.

## Validation and restoration

- `bash tests/run_all_tests.sh`: 28 modules successful, **208 tests passed, 1 skipped**.
- `bash -n tests/run_all_tests.sh`: passed.
- `git diff --check`: passed.
- All 47 files exist at recorded trash paths with matching size/inode, and original paths are absent.

No annotation workload or HPC command was run. Cluster fixture validation remains pending. Git reports deletions at original tracked paths because `trash/` is ignored; the files remain physically available in trash and in committed history.

Restore only selected manifest `moved_paths` to their original paths after checking for conflicts; never overwrite newer files. The trash directory includes `RESTORE.md`. Moving to trash retains disk usage and does not make deletion implicit.
