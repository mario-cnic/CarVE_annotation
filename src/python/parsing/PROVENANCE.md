# `src/python/parsing/` — VCF-to-table parsing code (vendored from `shared/utils`)

Verbatim copy, made 2026-10-01 (user-authorized), of the files the Nextflow `VCF_TO_TABLE` step
(and the legacy `src/hpc/vcf2parsed.sh`) run. From here on **this copy is the one this pipeline
versions**; it is expected to diverge from `/data_lab_PGP/shared/utils/src/`.

| | |
|---|---|
| Source | `/data_lab_PGP/shared/utils/src/` (its own git repo) |
| Source commit | `19925ed2edad8298e0d84c4384f1fc7454e0074e` (2026-09-10 22:19:38 +0200); no uncommitted changes to `src/` at copy time |
| Copy verified | sha256 of all 11 files identical to source (below); `py_compile` OK |

Files: `vcf_parser_pysam.py`, `logger.py`, `filter_variants.py`, `modules/{__init__,acmg,disease_hpo,io_qc,pedigree,predictors,scoring,splicing}.py`.
External dependency not vendored: the `vcf_parser` pip package inside the `vcf_parser` conda env.

```
0ceb26767750ce22f0cf779260e010c1269b6427484fc1bdb838edc1d295ccb7  vcf_parser_pysam.py
8d0a8163ae1e13830d8459004ac35dcb2a08ca78c76d1b7fd2ec07aaef79db32  logger.py
ac5717c73b69d9fb460d48aac4f814f1a84b37d005d2a7fe76f4b3e4068fabbb  filter_variants.py
5e26b95dd42883462ed633489271aaec34cec04e5d17d596672424acce85ab6e  modules/acmg.py
daf6044d4fc9510898dc76e7a0de501ac58edeb5bda425eeb232190e06ddf811  modules/disease_hpo.py
5a94fa72ffe6063d5cc20f1cb90ade6e5beeb34e0015cf5c29ecc9af3cb90d12  modules/__init__.py
975d16acd0943a7851b9d4bf5fc56a5b97f5dff9417695cc2406315c37ba8a44  modules/io_qc.py
1bf4de98537edbdca2f3f8a7e55f3f685e358a3a3007e9645b0261eb6663627d  modules/pedigree.py
0067555ff504ad06407389b3edb18b4569b0ebc45aaa0a3b48560834c8aaef5c  modules/predictors.py
2a0b23fa69e2198fdd93e2f471dec193bfb0dac99173564789865b65ae4a5e29  modules/scoring.py
67562aef8ab3760504c49dc5590b3fb2a98230b449b49e7e61f9bee7f4355efb  modules/splicing.py
```

## Known caveats (documented, not fixed here)

- **Module 3 has its own, already-diverged copy** (`clinical_variant_prioritization/src/filter_variants.py`,
  `src/modules/splicing.py`; e.g. the shared copy has `parse_maxentscan`, Module 3's does not).
  Three copies now exist (shared, Module 3, this one). Changes made here do not propagate.
- **Open risk — silent undo of gene/transcript-matched scores.** `parse_pangolin`, `parse_spip` and
  `parse_spliceai_custom` (`modules/splicing.py`) recompute scores from the raw `Pangolin`/`SPiP`/`SpliceAI`
  string columns whenever those are present (Pangolin's pre-computed `Pangolin_max_score` is only used when the
  raw string is empty). If Module 3 (or anything) re-runs its `parse_*` on a table from this pipeline, any
  attribution fix applied upstream would be recomputed from the raw strings and lost. Not yet checked whether
  Module 3's app does this. See `TODO.md` (predictor gene/transcript attribution item).
- **Legacy bash path not repointed:** `config/env.sh` (`VCF_PARSER_SCRIPT`, `FILTER_VARIANTS_SCRIPT`) and
  `src/tools/backfill_variant_prioritization.py` still point at `shared/utils/src`. Only the Nextflow
  config uses this copy.

## Local modifications since the verbatim copy (2026-10-01, this repo only)

The files above were first committed byte-identical to `shared/utils@19925ed` (`d5a915c`); the sha256 table
describes that state. Since then, gene/transcript attribution was added (plan:
`walkthrough/20261001_predictor_gene_transcript_attribution_audit.md`):

- new `modules/gene_identity.py` (symbol/accession -> Ensembl gene resolution, per-row match statuses);
- `modules/splicing.py`: `parse_spliceai_custom`, `parse_spip`, `parse_pangolin` now attach only the entry of the row's
  own gene/transcript and add `*_match` and `*_anygene_*` columns; `_extract_spliceai_details` uses VEP's
  per-transcript SpliceAI score only when the custom score is absent (a genuine 0.00 custom score is kept);
- `filter_variants.py`: new `--hgnc-table`, `--spliceai-symbol-map`, `--gene-transcript-mapping`; refuses to run
  with SpliceAI/SPiP columns present but no identity resources (no silent fallback).

Later changes (same day, `f805585`, `a0c1f4f`, `05a470a`):
- `modules/predictors.py`: `parse_dbnsfp_by_row_transcript` + `load_aligned_columns` - each row keeps the element of its own
  transcript for the 67 transcript-aligned dbNSFP columns (`resources/dbnsfp_transcript_aligned_columns.txt`, classified by
  `src/tools/derive_dbnsfp_dependency.py` from the dbNSFP readme + raw lines); adds `dbNSFP_match`, `dbNSFP_transcripts_all`,
  `<col>_anytranscript_max`; runs before `parse_missense`.
- `modules/splicing.py` `parse_spip`: accession match, else exon-structure match (`matched_by_exon_structure`), plus labelled
  `SPiP_samegene_*` columns.
- `modules/gene_identity.py`: `GeneIdentity(spliceai_symbol_map, hgnc_table, enst_spip_map)` - all three are required; the
  curated `Gen,NM,ENST` pairing is not used for SPiP.
- `filter_variants.py`: required arguments `--hgnc-table`, `--spliceai-symbol-map`, `--enst-spip-map`, `--dbnsfp-aligned-columns`
  (each only when the corresponding predictor columns are present).

`git diff d5a915c -- src/python/parsing` shows the exact delta. The shared and Module 3 copies are unchanged.

Row order (2026-10-06): `filter_variants.py` no longer sorts the table (it sorted by `NEW_IMPACT`, `SPiP_prediction`, frequency); rows keep the input's
coordinate order, checked with `modules/io_qc.py` `check_coordinate_order` (a warning is logged when not ordered). `sort_by` stays in `io_qc.py`, unused.
The sha256 list above describes the 2026-10-01 copy and is out of date for these two files.

Population AF update (2026-10-07): new `modules/population.py` centralizes custom v4 grpmax → joint per-row resolution; `modules/acmg.py` and `modules/scoring.py` preserve unknown AF rather than assigning zero. Unknown AF earns no frequency evidence/bonus and is retained with availability labels by rarity filters. Source columns remain unchanged. The shared-utils and Module 3 copies remain unchanged; see `walkthrough/20261007_population_frequency_missingness.md`.

ClinVar source priority (2026-10-08T10:03:28+02:00): new `modules/clinvar.py` and updates to `acmg.py`/`scoring.py` prefer custom VCF per row, then dbNSFP/cache and legacy aliases; review stars stay with the selected source. See `walkthrough/20261008_clinvar_source_priority.md`. Shared-utils and Module 3 copies are unchanged.

SpliceAI missingness (2026-10-08T10:19:13+02:00): `modules/splicing.py` keeps the VEP-only maximum distinct from custom and preserves unknown selected scores; `modules/scoring.py` exposes selected-score availability and retains other predictor contributions. See `walkthrough/20261008_spliceai_missingness.md`. Shared-utils and Module 3 copies remain unchanged.

Legacy cleanup (2026-10-08T11:53:30+02:00): `src/hpc/vcf2parsed.sh` and `config/env.sh` are now retained only under `trash/legacy_pipeline_20261008_115330/`. Active Nextflow continues to use this vendored parser. Historical source/copy statements above remain evidence of their dates.

## Local changes

- 2026-10-14: row-wise passes in `modules/{io_qc,splicing,pedigree,acmg,scoring,disease_hpo}.py` made faster (tuple results through `io_qc.apply_rows`, arrays instead of `.iloc` in loops); output identical to the previous version on local fixtures (`walkthrough/20261014_performance_review.md`).

Dead code removal (2026-10-15T15:48:51+02:00): `io_qc.sort_by` and `vcf_parser_pysam.add_sample_header` (unused) were removed, along with unused imports; see `walkthrough/20261015_code_bloat_cleanup.md`. The sentence above saying `sort_by` stays is superseded.
