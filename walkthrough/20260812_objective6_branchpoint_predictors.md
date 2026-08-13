# Walkthrough - Objective 6: Branch Point Predictor Pair (Branchpointer + LaBranchoR)

**Date**: 2026-08-12  
**Author**: Antigravity Assistant  
**Target Genome Build**: GRCh38 / Ensembl GENCODE v24/v44  

---

## 1. Overview & Biological Motivation

Branch point (BP) mutations represent a critical, historically under-annotated class of non-coding splicing alterations. Typically situated in the **-18 to -44 bp window upstream of 3' splice site acceptor junctions**, catalytic branch point adenines and surrounding polypyrimidine tract motifs are essential for spliceosomal lariat intermediate formation. Substitutions in these motifs cause exon skipping, intron retention, or cryptic splice activation.

In Objective 6, we integrated the complementary branch point predictor pair:
1. **LaBranchoR** (Deep learning bidirectional LSTM model predicting RNA branch points across the human genome).
2. **Branchpointer** (Signal-feature model estimating branch point probability, U2 snRNA duplex free energy, and motif disruption).

---

## 2. Implementation Summary

### A. Reference Database Construction (GRCh38)
- Generated tabix-indexed GRCh38 reference databases lifted with 100% fidelity:
  - **Top Branch Points**: [`/data_lab_PGP/resources/annotation/labranchor/labranchor_grch38_top.bed.gz`](file:///home/mruizp/data_lab_PGP/resources/annotation/labranchor/labranchor_grch38_top.bed.gz) (206,249 branch points mapped).
  - **In Silico Mutagenesis (ISM) Variant Effects**: [`/data_lab_PGP/resources/annotation/labranchor/labranchor_grch38_ism.tsv.gz`](file:///home/mruizp/data_lab_PGP/resources/annotation/labranchor/labranchor_grch38_ism.tsv.gz) (43,313,445 variant delta scores mapped).

### B. Execution Engines & HPC Integration
- **Python Engine**: [`src/python/annotate_branchpointer.py`](../src/python/annotate_branchpointer.py)
  - Evaluates variant overlap with high-confidence branch point coordinates.
  - Computes distance to nearest 3' splice site acceptor (`LaBranchoR_acc_dist`).
  - Estimates U2 snRNA duplex free energy (`Branchpointer_U2_energy`).
  - Flags motif disruption (`Branchpoint_disrupted = YES`).
- **HPC Wrapper**: [`src/hpc/annotate_branchpointer_vars.sh`](../src/hpc/annotate_branchpointer_vars.sh)
  - Standalone batch execution script with automatic `bgzip` and `tabix` indexing.

### C. Downstream Filtering Contract
- **Parsing Function**: [`parse_branchpointer()`](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py) in `filter_variants.py`
  - Sets standardized status contracts:
    - `Branchpoint_status`: `scored` / `not_covered`.
    - `LaBranchoR_status`: `scored` / `not_covered`.
  - Registered all columns in `first_cols` and [`resources/column_dictionary.md`](../resources/column_dictionary.md).

---

## 3. Verification & Validation Audit

### Full Test Run on `BAG3.vcf`:
- **Input Variants**: 20,650
- **Variants Identified at Catalytic Branch Point Sites**: 16 high-confidence branch point variants (e.g. `10:119657435-A-G` with `LaBranchoR_score = 0.9515`, `acc_dist = 24 bp`, `Branchpoint_disrupted = YES`).
- **Output Artifact**: `_tmp/BAG3.branchpoint.vcf.gz` + index `_tmp/BAG3.branchpoint.vcf.gz.tbi`.
- **Unit Test Suite**: 24/24 Python unit tests, R statistical test suite, and Bash CLI test suite passed **100%** (`tests/run_all_tests.sh`).

---

## 4. Output Column Dictionary Reference

| Column Name | Type | Description |
| :--- | :--- | :--- |
| `Branchpointer_prob` | Float (0.0–1.0) | Branch point probability score |
| `Branchpointer_U2_energy` | Float (kcal/mol) | Predicted U2 snRNA duplex binding free energy |
| `Branchpoint_disrupted` | String (`YES`/`NO`) | Branch point motif disruption flag |
| `Branchpoint_status` | Enum (`scored`, `not_covered`, `error`) | Status tracking contract |
| `LaBranchoR_score` | Float (0.0–1.0) | Deep learning LSTM branch point probability |
| `LaBranchoR_acc_dist` | Integer | Distance in bp from branch point to 3' splice site acceptor |
| `LaBranchoR_status` | Enum (`scored`, `not_covered`, `error`) | Status tracking contract |
