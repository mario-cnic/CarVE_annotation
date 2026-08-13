# Walkthrough - Objective 5: Local SpliceAI at `-D 10000` Integration

**Date**: 2026-08-12  
**Author**: Antigravity Assistant  
**Target Genome Build**: GRCh38 / Ensembl GENCODE v24  

---

## 1. Overview & Motivation

Pre-computed VEP SpliceAI tabix lookup tables are constrained to a 50 bp window around canonical splice junctions. To capture deep intronic splice-altering mutations and pseudo-exon creations across non-coding intervals, we deployed local deep learning SpliceAI annotation dynamically scoring variants within a **±10,000 bp window (`-D 10000`)**.

---

## 2. Implementation Summary

### A. Environment Architecture
- **Dedicated Environment**: `spliceai_env` (Python 3.10, TensorFlow 2.15.1, `spliceai==1.3.1`, `setuptools<70`, `pysam`, `htslib`).
- **Reproducible Specification**: Saved at [`resources/conda/spliceai_env.yml`](../resources/conda/spliceai_env.yml).

### B. Execution Architecture
- **Python Engine**: [`src/python/annotate_spliceai.py`](../src/python/annotate_spliceai.py)
  - Loads 5 ensemble model weights (`spliceai1.h5` through `spliceai5.h5`).
  - Evaluates variant delta scores across all 4 channels (Acceptor Gain, Acceptor Loss, Donor Gain, Donor Loss).
  - Emits standard VCF INFO tag: `SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`.
- **HPC Execution Wrapper**: [`src/hpc/annotate_spliceai_vars.sh`](../src/hpc/annotate_spliceai_vars.sh)
  - Manages thread allocation, paths, automatic `bgzip`, and `tabix` indexing.

### C. Downstream Filtering Contract
- **Parsing & Max Score**: [`shared/utils/src/filter_variants.py`](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py)
  - Evaluates `spliceAI_MAX` as maximum delta score across all 4 channels (`DS_AG`, `DS_AL`, `DS_DG`, `DS_DL`).
  - Sets standardized contract column `SpliceAI_status`:
    - `scored`: Successfully evaluated by SpliceAI model.
    - `not_covered`: Variant outside covered gene boundary or unannotated contig.

---

## 3. Verification & Validation Audit

### Full Test Run on `BAG3.vcf`:
- **Input Variants**: 20,650
- **Variants Successfully Annotated with SpliceAI**: 18,545
- **Output Artifact**: `_tmp/BAG3.spliceai.vcf.gz` (343 KB) + `_tmp/BAG3.spliceai.vcf.gz.tbi` (266 B).
- **Unit Test Suite**: 22/22 Python unit tests, R statistical test suite, and Bash CLI test suite passed 100% (`tests/run_all_tests.sh`).

---

## 4. Output Column Dictionary Reference

| Column Name | Type | Description |
| :--- | :--- | :--- |
| `spliceAI_MAX` | Float (0.0–1.0) | Maximum delta score across all 4 SpliceAI channels |
| `SpliceAI_status` | Enum (`scored`, `not_covered`) | Status contract tracking whether variant was evaluated by SpliceAI |
| `SpliceAI_pred_DS_AG` | Float (0.0–1.0) | Delta score for Acceptor Gain |
| `SpliceAI_pred_DS_AL` | Float (0.0–1.0) | Delta score for Acceptor Loss |
| `SpliceAI_pred_DS_DG` | Float (0.0–1.0) | Delta score for Donor Gain |
| `SpliceAI_pred_DS_DL` | Float (0.0–1.0) | Delta score for Donor Loss |
| `SpliceAI_pred_DP_AG` | Integer | Distance to predicted Acceptor Gain site |
| `SpliceAI_pred_DP_AL` | Integer | Distance to predicted Acceptor Loss site |
| `SpliceAI_pred_DP_DG` | Integer | Distance to predicted Donor Gain site |
| `SpliceAI_pred_DP_DL` | Integer | Distance to predicted Donor Loss site |

All columns are registered in [`resources/column_dictionary.md`](../resources/column_dictionary.md) and [`resources/all_but_old_gnomad_vep_cols.txt`](../resources/all_but_old_gnomad_vep_cols.txt).
