# Module 2 no longer writes variant scores, tiers or ACMG tags (SF-06, SF-07)

Updated: 2026-11-24T16:00:00+01:00

## What changed

New tables lose eight columns; every other column keeps its values.

| Removed | Why |
| --- | --- |
| `PRIORITY_TIER`, `VARIANT_PRIORITY_SCORE`, `TIER_RATIONALE`, `SCORE_EXPLANATION`, `ACMG_CRITERIA`, `HPO_MATCH` | Old scoring logic (ACMG-1/2/3: PP3_Strong from one predictor, two BA1/BS1 policies, an uncalibrated score); Module 3 recomputes and overwrites all of them. Part 1 of DEC-0008. |
| `Pangolin_heart_lv_score`, `Pangolin_heart_aa_score` | Copies of `Pangolin_max_score`, not tissue scores (PRED-1). |

Code: `src/python/parsing/modules/scoring.py` `build_priority_tier` is replaced by `build_interpretation_inputs`, which keeps the annotation the old function also wrote: gene curation (`apply_disease_phenotype_curation`), `POPULATION_AF*`, `SPLICE_MAX_UNIFIED` / `SPLICEAI_METHOD` / `SPLICEAI_STATUS`, and `QC_STATUS` when missing. `acmg.py` `build_acmg_criteria` and `PVS1_VETO_GENES` are removed (no other user). `src/tools/backfill_variant_prioritization.py` (wrote these columns into old tables through `shared/utils`) moved to `archive/deprecated_scripts/`.

Kept, although DEC-0008 lists them as Module 3's: `ACMG_SF_V3_2`, `INHERITANCE_MODEL`, `DISEASE_*`, `HIGH_RISK_DISEASE_GENE`, `NEW_IMPACT`, `QC_STATUS`, `QC_DETAIL`. Not decided here.

## Checks

- **Same table minus eight columns.** Stage 1 + stage 2 rebuilt locally from `nf_work/annotation_out/sf01_nmd_fixture/annotation/sf01_nmd_fixture.annotated.vcf.gz` with the code of `0e5716f` and with this change (`shared/utils` `vcf_parser` and `datasci` envs, flags of the pipeline): 1,015 rows each, same `Locus`/`Feature` order, exactly the eight columns above removed, no other column differs. The `0e5716f` build equals the cluster table of that run on every shared column.
- **Module 3 gives the same result without them.** Module 3 at `3266093`, `prepare_scored_frame` (the path of uploads and of candidate sets from a table), on both local builds: `PRIORITY_TIER`, `VARIANT_PRIORITY_SCORE`, `C_PATHSCORE`, `ACMG_CRITERIA`, `TIER_RATIONALE`, `SCORE_EXPLANATION`, `HPO_MATCH` identical on all 1,015 rows; no other output column differs except the two Pangolin heart columns, absent from the new one.
- **Module 3 readers** (read 2026-11-24): uploads, candidate sets, rescore and the demo builder all call `build_priority_tier` before reading these columns (`src/api/routers/variants.py:273,588`, `scripts/build_unified_demo_dataset.py:923`, `src/filter_variants.py:350`); a scoring failure returns HTTP 500 instead of defaults (`variants.py:276-278`); `engine/column_sets.py:63` fetches only columns present in the table; `io_qc.py:526-529` only warns on missing columns.
- Local suite: 462 passed, none skipped. Tests retargeted to the kept contracts (`POPULATION_AF_STATUS`, `CLINVAR_SOURCE`, `SPLICEAI_STATUS`); removed: `test_pvs1_veto.py` (4 tests, scoring only) and `test_unknown_cannot_promote_to_tier_one`.

## Pending

- Cluster fixture (covers SF-06 and SF-07): new table has none of the eight columns and equals `sf01_nmd_fixture`'s table on every other column.
- Tables built before this change (`F2005_union`, `F1830_union_avi`, their LAB copies) still carry the columns; Module 3 overwrites them on scoring, so nothing to redo.
