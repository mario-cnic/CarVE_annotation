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

`git diff d5a915c -- src/python/parsing` shows the exact delta. The shared and Module 3 copies are unchanged.
