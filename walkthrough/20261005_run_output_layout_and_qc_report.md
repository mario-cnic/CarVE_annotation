# Run output sub-folders and run QC report (2026-10-05)

Item: `TODO.md` "NEXT SESSION - START HERE" (Mario, 2026-10-05). Layout and report format confirmed by Mario 2026-10-05 ("go with the proposals"): HTML + JSON + TSV, no MultiQC.

## What changed
- **Layout.** New runs publish into `<run_id>/{annotation,predictors,tables,genotypes,reports,pipeline_info}/`. The names are `params.out_subdir` in `nextflow.config`; each module's `publishDir` uses one of them. Old runs keep the flat layout (nothing moved). Modules shared with the legacy per-gene `main.nf` follow the same layout.
- **Nextflow reports.** `trace` (explicit field list incl. `hostname`, `peak_rss`), `timeline` and `report` are enabled, written to `pipeline_info/` with a launch timestamp. `pipeline_info_dir` falls back to the input filename without `.gz` / `.vcf` when `--run_id` is not given (the launcher's rule; `annotate_vcf.nf` uses `baseName`, which differs for `.vcf.bgz`, so a direct `nextflow run` of such a file would put `pipeline_info/` under another name than its outputs).
- **Launcher.** `run_annotate_vcf.sh` passes `-log <run>/pipeline_info/nextflow_<ts>.log` (the repo-root `.nextflow.log` is rotated by every launch), then runs `build_run_qc_report.py`; a QC failure does not change the exit code.
- **Manifest.** `write_run_manifest.py` keeps `--outdir` = run folder, writes `RUN_MANIFEST.json` / `RUN_DIRTY.*.patch` into `pipeline_info/`, and reads the check report and genotype table from the sub-folder first, then the flat layout (`run_files`).
- **QC report.** `src/python/build_run_qc_report.py`: run summary, record counts per stage (tabix index counts, Parquet metadata, column-limited reads for `*_status`, tier and distinct `Locus`), time / memory / nodes per process, failed attempts with the last stderr line, provenance. Task data come from the latest `execution_trace_*.txt`, else from a `.nextflow.log` (runtime and node from the work dirs; no memory). Missing sources print as `unavailable`, never 0. Output `reports/<run_id>.run_qc.<ts>.{html,json,processes.tsv}`, never overwritten.

## Checks (local, FACT)
- `-stub-run -profile local_dev` on the 3-sample test fixture before and after: the same 16 files, each in its sub-folder, plus the trace / timeline / report; with and without `--run_id`.
- `tests/python`: 112 passed, 2 modules not collected (`test_combine_valencia_hic.py` needs `genebe`; `test_predictor_attribution.py` imports `modules` from `shared/utils/src`, which shadows the vendored copy): both environment problems of the local `datasci` env, not caused by this change. New / updated: `test_build_run_qc_report.py`, `test_write_run_manifest*.py` (flat-layout case added).
- QC script on the finished family run (flat layout, log fallback, output to scratch): matched the figures already recorded (914 tasks, 11 failed attempts, 56.6 h, 5,263,518 input / merged records, 23,686,438 table rows, 15,790,554 genotype rows, the Pangolin chunk killed at 36 h, `CONCAT_VEP` killed at 1 h). 2.5 min on a cold read, 37 s warm, peak 2.5 GB RSS (with the distinct-`Locus` count; `--skip-distinct-locus` avoids reading that column). The wrapper runs it on the login node. Not written into the run folder.
- On the 3-sample trio run (flat, no manifest, no log) every missing section is `unavailable`.

## Not verified
- Not run on the cluster: the launcher end to end, the trace `hostname` / `peak_rss` fields under SGE (they are `-` locally), and `nextflow.config` under the cluster's Nextflow 25.10.2 (checked here with 24.04.2 only).
- Full `tests/python` suite does not collect in the local `datasci` env (missing `genebe`, vendored-parser import); unrelated to this change.
- Runs resumed with `-resume` from the flat layout would re-publish cached outputs into the sub-folders (duplicates of large files); do not do that.

## Evidence for the MISC-16 / MISC-17 checks (from the family-run process table, no new analysis)
- No SpliceAI or Pangolin task ran on `c0053-cn1`, so the node exclusion held; the one Pangolin attempt killed at the 36 h limit ran on `c0051-cn1`.
- `CONCAT_VEP` still had one hung first attempt, killed at 1.0 h and retried: the MISC-16 mitigation shortened the hang (2 h before) but did not prevent it.
