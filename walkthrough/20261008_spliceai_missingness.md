# SpliceAI two-source missingness

Updated: 2026-10-08T10:19:13+02:00

## Scope and change

User authorized fixing SpliceAI after checkpoint commit `4f3b7ab`. Scope is the active Module 2 parser and scoring; no HPC workload, reference, existing output, or sibling repository was changed.

FACT: VEP lookup fields (`SpliceAI_pred_*`) and custom INFO (`SpliceAI`, parsed `spliceai_custom_*`) are both exported. `SPLICING_COLUMNS` previously mixed their delta scores into `spliceAI_MAX`. The VEP maximum was then filled with zero, and source labels depended on whether a column existed rather than whether that row had a prediction. Two absent sources could therefore produce zero and a VEP/custom label.

Implemented in `src/python/parsing/modules/splicing.py`:

- VEP-only component maximum, with numeric coercion and finite 0–1 validation. Custom components cannot enter this maximum, even if passed in the caller's list.
- Missing maximum remains null. Genuine zero is scored. Known source errors are preserved.
- Select valid custom maximum first, including zero; otherwise valid VEP maximum. Two unavailable sources produce null and `Not Covered`.
- Raw source annotations, custom match reasons and both source-specific maxima remain exported. Custom gene/allele matching remains unchanged.

`scoring.py` exports selected-score `SPLICEAI_STATUS` (`scored`, `not_covered`, or explicit `error`). The existing `SPLICE_MAX_UNIFIED` and `SPLICEAI_METHOD` now preserve missingness/source consistently. The splicing composite ignores missing SpliceAI while retaining valid SPiP/Pangolin contributions, so the continuous score remains numeric. Unknown SpliceAI cannot satisfy the existing low-score requirements for BP4/BP7. Existing thresholds and score weights were not recalibrated.

## Validation

Local synthetic checks only:

- `python -m pytest -q tests/python/test_spliceai_missingness.py tests/python/test_population_frequency.py tests/python/test_clinvar_source_priority.py tests/python/test_multisample_trio.py tests/python/test_row_order.py tests/python/test_vep_custom_columns.py`: **60 passed**.
- `python -m pytest -q tests/python/test_predictor_attribution.py`: **26 passed**, separately to avoid existing top-level package collisions.
- `git diff --check`: passed.

Cases include both exports preserved, source maxima separated, genuine custom/lookup zeros, row-wise fallback, neither source/absent columns, numeric strings, malformed/nonfinite/out-of-range values, nondefault index, previously zero-imputed uncovered predictions, no unsupported benign evidence, other-predictor score preservation, explicit errors and valid custom rescue.

## Remaining validation

Run a fresh small cluster fixture into a new directory and inspect raw/source-specific/selected-score columns. Existing tables are not repaired automatically; a historical zero without status/raw components cannot retrospectively be distinguished from a real zero. Shared-utils and Module 3 copies remain unchanged; hand off this contract before independent interpretation recomputes scores. VEP lookup window labels are retained from the existing code; this fix did not independently validate the precomputed resource's window metadata.
