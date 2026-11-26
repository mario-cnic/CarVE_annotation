# Module 2 no longer writes variant scores, tiers or ACMG tags (SF-06, SF-07)

Updated: 2026-11-26T12:30:00+01:00

## What changed

New tables lose eight columns; every other column keeps its values.

| Removed | Why |
| --- | --- |
| `PRIORITY_TIER`, `VARIANT_PRIORITY_SCORE`, `TIER_RATIONALE`, `SCORE_EXPLANATION`, `ACMG_CRITERIA`, `HPO_MATCH` | Old scoring logic (ACMG-1/2/3: PP3_Strong from one predictor, two BA1/BS1 policies, an uncalibrated score); Module 3 recomputes and overwrites all of them. Part 1 of DEC-0008. |
| `Pangolin_heart_lv_score`, `Pangolin_heart_aa_score` | Copies of `Pangolin_max_score`, not tissue scores (PRED-1). |

Code: `src/python/parsing/modules/scoring.py` `build_priority_tier` is replaced by `build_interpretation_inputs`, which keeps the annotation the old function also wrote: gene curation (`apply_disease_phenotype_curation`), `POPULATION_AF*`, `SPLICE_MAX_UNIFIED` / `SPLICEAI_METHOD` / `SPLICEAI_STATUS`, and `QC_STATUS` when missing. `acmg.py` `build_acmg_criteria` and `PVS1_VETO_GENES` are removed (no other user). `src/tools/backfill_variant_prioritization.py` (wrote these columns into old tables through `shared/utils`) moved to `archive/deprecated_scripts/`.

Kept, although DEC-0008 lists them as Module 3's: `NEW_IMPACT`, `QC_STATUS`, `QC_DETAIL` (DEC-0008 item 2). The curation and inheritance columns went in SF-08, below.

## SF-08: gene curation and inheritance columns (2026-11-24)

Removed from new tables: `ACMG_SF_V3_2`, `DISEASE_GENE_EVIDENCE`, `DISEASE_CURATION_SOURCE`, `DISEASE_PHENOTYPE_MATCH`, `DISEASE_EVIDENCE_SCORE`, `HIGH_RISK_DISEASE_GENE`, `INHERITANCE_MODEL`, `SAMPLE_GENOTYPES_SUMMARY`. In Module 2 tables the last two were always `Unclassified` and empty (no genotype columns in the table).

- Code: `filter_variants.py` no longer calls `analyze_pedigree_inheritance` or the curation; `build_interpretation_inputs(data, freq_col)` keeps only frequency status, unified SpliceAI and QC. `disease_hpo.py` keeps only the `--genes_file` gene-priority filter (`get_dynamic_priority`, `prioritize_genes`, `filter_by_gene_priority`), which never used the curation; the curation loader, core panels, ACMG SF list and HPO panels are gone.
- The disease-gene list is no longer an input of Module 2: `--disease-gene-list`, `params.disease_gene_list` (`Complete_gene_list_V4.csv`), its entry in the `tsv_to_table` cache key and in the run manifest resources are removed. Module 3 curates from its own `Complete_gene_list_V5` (`clinical_variant_prioritization/src/config.py:49-51`).
- Checks: same local rebuild of the `sf01_nmd_fixture` VCF, against the SF-06 build: 1,015 rows, same order and column order, exactly these eight columns removed, nothing else differs. Module 3's `prepare_scored_frame` on the `0e5716f` build and on this one: every shared output column identical; the output lacks `SAMPLE_GENOTYPES_SUMMARY` (Module 3 writes it only from genotype columns; its frontend shows it only when non-empty, `frontend/src/components/FormattedCells.tsx:234`) and the Pangolin heart columns. Suite: 458 passed, none skipped (removed `test_disease_gene_curation.py`, `test_disease_core_panels.py` and the `mini_disease_gene_list.csv` fixture; `--disease-gene-list` is now a rejected option in `test_filter_variants_cli.py`). `nextflow config` has no `disease_gene_list`; `annotate_vcf.nf -preview` parses (stops at the launcher-only cache check).
- Cluster fixture **passed 2026-11-24**: `sf08_fixture` (commit `6d1efec`, clean tree, exit 0, 21/21 tasks); 1,015 rows in the same order, exactly the eight columns removed, every other column identical to `sf06_sf07_fixture`'s table; the run manifest no longer lists `disease_gene_list`.

## SF-09: impact and QC columns (2026-11-26)

Removed from new tables: `NEW_IMPACT`, `QC_STATUS`, `QC_DETAIL` (DEC-0008 item 2). Module 3 computes all three on ingest since `clinical_variant_prioritization@f5663ac` and reported its output identical with and without them (`sf08_fixture`, `utr5_multieffect_run1`).

- Code: `build_newImpact` (`acmg.py`), `build_qc_status` (`io_qc.py`), `utr5_high_impact_mask` (`predictors.py`) removed. `--classification dominant` (`classification_dominant`, `Variant_rank` from gene priority x `NEW_IMPACT` x frequency) removed with them: the pipeline never passed it and it needed `NEW_IMPACT`.
- `build_newImpact` also selected the ClinVar record as a side effect (`CLINVAR_SELECTED_CLNSIG`, `CLINVAR_SOURCE`, `CLINVAR_SELECTED_REVIEW`, `ClinVar_Stars`). `build_interpretation_inputs` now calls `annotate_clinvar` itself, so those columns stay. Found while writing the contract test: the first version of this change dropped them.
- Contract: the columns Module 3 reads to compute impact and QC are listed in `scoring.py` (`MODULE3_QC_INPUTS`, `MODULE3_IMPACT_INPUTS`) and in `resources/column_dictionary.md`; `test_module3_input_columns_present` fails if one leaves the end-to-end table (SPiP names checked as constants: the fixture has no SPiP INFO). All 18 are in the `sf08_fixture` cluster table.
- Checks: local rebuild of the same panel8 stage-1 TSV with this code against the SF-08 build: 1,015 rows, same row and column order, exactly the three columns removed, no other column differs. Module 3 at `f5663ac`, `prepare_scored_frame` on both builds: same 680 output columns, all identical. Suite: 457 passed (removed the `build_newImpact` tests in `test_utrannotator_status.py`, `test_branchpoint_predictors.py` and `test_clinvar_source_priority.py`, and the `classification_dominant` check; `--classification` is now a rejected option).
- Cluster fixture **passed 2026-11-26**: `sf09_fixture` (commit `1015fa2`, clean tree, exit 0, 21/21 tasks); 1,015 rows in the same order, exactly the three columns removed, every other column identical to `sf08_fixture`'s table and to the local build; all 18 Module 3 input columns present.

## Checks

- **Same table minus eight columns.** Stage 1 + stage 2 rebuilt locally from `nf_work/annotation_out/sf01_nmd_fixture/annotation/sf01_nmd_fixture.annotated.vcf.gz` with the code of `0e5716f` and with this change (`shared/utils` `vcf_parser` and `datasci` envs, flags of the pipeline): 1,015 rows each, same `Locus`/`Feature` order, exactly the eight columns above removed, no other column differs. The `0e5716f` build equals the cluster table of that run on every shared column.
- **Module 3 gives the same result without them.** Module 3 at `3266093`, `prepare_scored_frame` (the path of uploads and of candidate sets from a table), on both local builds: `PRIORITY_TIER`, `VARIANT_PRIORITY_SCORE`, `C_PATHSCORE`, `ACMG_CRITERIA`, `TIER_RATIONALE`, `SCORE_EXPLANATION`, `HPO_MATCH` identical on all 1,015 rows; no other output column differs except the two Pangolin heart columns, absent from the new one.
- **Module 3 readers** (read 2026-11-24): uploads, candidate sets, rescore and the demo builder all call `build_priority_tier` before reading these columns (`src/api/routers/variants.py:273,588`, `scripts/build_unified_demo_dataset.py:923`, `src/filter_variants.py:350`); a scoring failure returns HTTP 500 instead of defaults (`variants.py:276-278`); `engine/column_sets.py:63` fetches only columns present in the table; `io_qc.py:526-529` only warns on missing columns.
- Local suite: 462 passed, none skipped. Tests retargeted to the kept contracts (`POPULATION_AF_STATUS`, `CLINVAR_SOURCE`, `SPLICEAI_STATUS`); removed: `test_pvs1_veto.py` (4 tests, scoring only) and `test_unknown_cannot_promote_to_tier_one`.

## Pending

- ~~Cluster fixture (covers SF-06 and SF-07)~~ **Passed 2026-11-24** (`sf06_sf07_fixture`, commit `7fb2e13`, clean tree, exit 0, 21/21 tasks; SpliceAI/Pangolin from the score cache, 78/78 hits, so 6 fewer tasks than `sf01_nmd_fixture`): 1,015 rows in the same order, exactly the eight columns removed, every other column identical to `sf01_nmd_fixture`'s table.
- Tables built before this change (`F2005_union`, `F1830_union_avi`, their LAB copies) still carry the columns; Module 3 overwrites them on scoring, so nothing to redo.
