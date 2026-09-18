# Pipeline Run Audit & Verification Walkthrough: `predictors_050826`

**Date**: 2026-08-07  
**Author**: Antigravity Assistant & Bioinformatics Pipeline Auditor  
**Target Run**: `RUNS/predictors_050826`  
**Genome Build**: GRCh38  

---

## 1. Overview & Execution Scope

A comprehensive audit was performed on the pipeline execution run `predictors_050826`, evaluating input variant files, step-by-step intermediate execution outputs, final parsed `.parsed.clean.pq` Parquet files, column schemas, annotation completeness, and error logs.

---

## 2. Audit Script & Methods

Created dedicated quality control & auditing engine:
- [src/python/audit_run_results.py](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/python/audit_run_results.py)

Executed via:
```bash
python3 src/python/audit_run_results.py \
  --run-dir RUNS/predictors_050826 \
  --raw-dir RUNS/predictors_050826/input/by_gene
```

---

## 3. Detailed Audit Results

### 📊 Executive Metrics
- **Input Genes**: 7 (`BAG3`, `DSP`, `FLNC`, `LMNA`, `MYBPC3`, `PKP2`, `TTN`)
- **Genes Fully Finished**: **7 / 7 (100.0%)**
- **Pending Genes**: **0**
- **Error Genes**: **0**
- **Total Input Variants**: **395,858**
- **Total Output Clean Variants**: **380,292**
- **Overall Variant Retention Rate**: **96.07%**
- **Column Schema Consistency**: **100% Consistent** across all 7 genes (**504 columns** in every file).

---

### 📋 Per-Gene Detailed Status

| Gene | Input Vars | Output Vars | Retention | Columns | Step Completion | Parquet Status | Logs Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **BAG3** | 20,650 | 20,650 | 100.0% | 504 | 6/6 | ✅ VALID | ✅ Clean |
| **DSP** | 29,784 | 29,784 | 100.0% | 504 | 6/6 | ✅ VALID | ✅ Clean |
| **FLNC** | 25,957 | 25,957 | 100.0% | 504 | 6/6 | ✅ VALID | ✅ Clean |
| **LMNA** | 41,098 | 25,533 | 62.1% | 504 | 6/6 | ✅ VALID | ✅ Clean |
| **MYBPC3** | 28,207 | 28,207 | 100.0% | 504 | 6/6 | ✅ VALID | ✅ Clean |
| **PKP2** | 63,777 | 63,777 | 100.0% | 504 | 6/6 | ✅ VALID | ✅ Clean |
| **TTN** | 186,385 | 186,384 | 99.99% | 504 | 6/6 | ✅ VALID | ✅ Clean |

---

## 4. Key Annotation Completeness (% Non-Null)

| Gene | REVEL | AlphaMissense | SPiP | SpliceVarDB | SpliceVault |
| --- | --- | --- | --- | --- | --- |
| **BAG3** | 4.9% | 4.9% | 100.0% | 0.0% | 0.0% |
| **DSP** | 13.3% | 13.3% | 100.0% | 0.1% | 0.3% |
| **FLNC** | 15.7% | 15.7% | 100.0% | 100.0% | 0.7% |
| **LMNA** | 4.3% | 4.2% | 100.0% | 0.0% | 0.1% |
| **MYBPC3** | 10.1% | 10.1% | 100.0% | 0.1% | 1.4% |
| **PKP2** | 2.3% | 2.2% | 100.0% | 0.0% | 0.1% |
| **TTN** | 27.3% | 26.9% | 100.0% | 0.0% | 1.4% |

---

## 5. Artifact & Report Locations

- **Full Markdown Audit Report**: [RUNS/predictors_050826/reports/audit_report_predictors_050826.md](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/predictors_050826/reports/audit_report_predictors_050826.md)
- **Final Parquet Results**: [RUNS/predictors_050826/results/](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/predictors_050826/results/)
- **Filtered Variant Datasets**: [RUNS/predictors_050826/filtered/](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/predictors_050826/filtered/)
- **Interactive Dashboards**: [RUNS/predictors_050826/reports/](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/predictors_050826/reports/)
- **Publication Figures**: [RUNS/predictors_050826/plots/](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/predictors_050826/plots/)
