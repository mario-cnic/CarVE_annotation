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
