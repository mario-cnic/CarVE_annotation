# Population-frequency missingness

Updated: 2026-10-07T15:41:00+02:00

## Scope and finding

FACT: Module 2 scoring previously converted missing AF to zero and could select an arbitrary column containing `af`/`freq` when its requested frequency column was absent. This granted rarity points and PM2-like evidence to unknown frequency. The standalone filter preferred generic `AF` over custom gnomADv4; the clinical report could substitute old cache/MAX_AF and displayed observed zero as missing. Existing results were not regenerated.

User authorized the ordered improvements in [TODO](../TODO.md). This change implements item 1 in this repository. No thresholds were recalibrated and no sibling repository was edited.

## Implemented contract

`src/python/parsing/modules/population.py` resolves a valid custom `gnomADv4_AF_grpmax_joint` value first, then `gnomADv4_AF_joint` per row. An explicitly selected alternative column has no fallback. Finite AF must be in [0, 1]. Reported zero remains zero; missing/invalid values remain null. No generic AF or older release is substituted automatically.

Derived columns: `POPULATION_AF`, `POPULATION_AF_SOURCE`, `POPULATION_AF_STATUS`. Status is `available_frequency`, `missing_frequency`, or `invalid_frequency`. These describe frequency availability, not population coverage or confirmed absence. Original annotation columns are preserved.

`modules/acmg.py` and `modules/scoring.py` use this resolution consistently. Unknown AF receives no PM2/BA1/BS1 frequency evidence or frequency score contribution; it cannot meet the existing frequency requirement for Tier 1. Other evidence can still prioritize the variant. Existing frequency thresholds and tier rules otherwise remain unchanged; the previously different ACMG-tag and disease-score thresholds remain a separate calibration question.

Rarity filters retain unknowns with labels and reject a requested AF filter if the expected columns are absent. `filter_and_summarize.py --max-af` stays inclusive; parser `filter_freq` stays strict. The standalone filter accepts `--freq-col`; the parser retains `--freq_col`.

The clinical report displays zero as `0`, unknown as `-`, and identifies a joint-AF fallback. Python interactive and R static AF plots use custom v4 grpmax/joint values; log plots include only positive AF and say so in the title. Legacy-only tables no longer generate a misleading gnomADv4 AF plot.

## Validation

FACT: local synthetic checks, no patient data or annotation workload:

- `python -m pytest -q tests/python/test_population_frequency.py tests/python/test_filter_and_summarize.py tests/python/test_generate_clinical_prioritization_report.py tests/python/test_generate_interactive_report.py tests/python/test_multisample_trio.py tests/python/test_row_order.py tests/python/test_vep_custom_columns.py`: **31 passed**.
- `python -m pytest -q tests/python/test_predictor_attribution.py`: **26 passed** (separate invocation to avoid the existing top-level `modules` import collision).
- `Rscript -e 'invisible(parse(file="src/R/plot_annotation_results.R"))'`: syntax passed.
- `Rscript tests/R/test_plot_annotation_results.R`: could not generate plots because local R lacks `ggplot2`; dependencies were not installed or changed.
- `git diff --check`: passed.

Cases cover explicit zero vs unknown, placeholders, malformed/out-of-range/infinite values, per-row joint fallback, conflicting legacy AF, missing columns, explicit override, all-unknown retention, strict/inclusive boundaries, no unknown rarity bonus, no unknown Tier 1 promotion, and clinical display.

## Remaining work

HPC behavior is unverified: run a fresh small cluster fixture into a new output directory and inspect the three derived columns and unknown/zero scoring behavior. R plot generation needs its configured dependencies. Existing tables retain their old derived interpretations until deliberately regenerated.

Module 3 and shared-utils copies remain unchanged. Module 3 needs the same availability/source contract before it recomputes frequency-derived evidence; this is tracked separately in TODO item 1a and the handoff. Next ordered item: dated custom ClinVar interpretation (item 2).

## Scope correction — 2026-10-07T16:44:06+02:00

Mario clarified that legacy reporting/plotting and downstream analyses are outdated. The R plot changes and corresponding fixture update made during this AF fix have been withdrawn; existing scripts are retained. The earlier R checks above describe the inspected implementation at that time, not a current requirement. R plotting validation is removed from the active TODO. Work continues on the Nextflow annotation/parser contract; SPiP remains an active R dependency.

## Correction 2026-10-08T15:36:21+02:00: absent from gnomAD is PM2 evidence

FACT: the change above also treated alleles with no gnomAD v4 record as unknown. Under `type=exact` custom annotation both fields are empty for such alleles, so they lost PM2_Supporting (+15) and could not reach Tier 1 (`scoring.py` requires PM2 for Tier 1). In the existing F2005_union table, 12 of 50 Tier 1 rows (4 variants) had no gnomAD v4 value. Mario (2026-10-08): absence from controls is ACMG evidence, not unknown.

Implemented in `modules/population.py`: status `absent_from_gnomad` when `gnomADv4_AF_joint` is in the table, it and `gnomADv4_AN_joint` are empty for the row, and grpmax gives no valid or invalid value. `POPULATION_AF` stays null (absence is not an observed zero). `rarity_af` maps absent rows to 0 for thresholds only: PM2_Supporting and +15 in scoring, PM2 in ACMG tags, retained by `filter_freq`, `Candidate`-eligible in `classification_dominant`. A gnomAD record with `AN_joint` 0 stays `missing_frequency`. A table without the joint column cannot show absence and stays `missing_frequency`.

Read-only status counts on F2005_union (23,686,438 rows; `gnomADv4_AF_grpmax_joint`, `gnomADv4_AF_joint`, `gnomADv4_AN_joint`): `available_frequency` 23,179,961 rows; `absent_from_gnomad` 504,934 rows / 186,304 variants; `missing_frequency` 1,543 rows / 537 variants (record with `AN_joint` 0).

Limitation: absence at a site gnomAD did not cover is indistinguishable from true absence. ClinGen SVI asks for adequate coverage before applying PM2; that needs a gnomAD coverage resource (TODO).

Validation: `bash tests/run_all_tests.sh`: 211 passed, 1 skipped. New/updated synthetic cases: absent vs AN 0 vs no joint column vs invalid; absent scores equal to observed zero and reach Tier 1; filter retention; dominant `Candidate`; availability audit. Not run on the cluster.
