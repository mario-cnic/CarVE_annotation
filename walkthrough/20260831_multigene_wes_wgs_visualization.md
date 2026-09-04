# 🌐 Multi-Gene / WES / WGS Genomic Visualization Walkthrough
**Date**: `2026-08-31 13:50:00 CEST`  
**Location**: `/data_lab_PGP/pipelines/annotation_pipeline_new/` (HPC) $\leftrightarrow$ `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/` (Local Sandbox)  

---

## 📌 Executive Summary

We resolved the multi-gene / Whole Exome (WES) / Whole Genome (WGS) plotting limitation by building an **Adaptive Genomic Visualization Engine**. 

When a multi-gene or WES/WGS dataset is loaded, the pipeline automatically detects that variants span multiple genes (`len(unique(SYMBOL)) > 1`) and seamlessly switches from the single-gene transcript model to a **Chromosomal Karyotype Manhattan & Gene Burden View**, while preserving an interactive **Per-Gene Drilldown Selector**.

---

## 🛠️ Key Implementation Details

1. **Adaptive Figure Engine (`create_multigene_manhattan_figure`)**:
   - Built [`create_multigene_manhattan_figure()`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/python/generate_clinical_prioritization_report.py#L290) in `generate_clinical_prioritization_report.py`.
   - **Subplot 1: Chromosomal Manhattan Scatter Plot**:
     - X-axis: Chromosomes (`chr1` to `chrX/Y/MT` ordered logically).
     - Y-axis: `VARIANT_PRIORITY_SCORE` (0-100).
     - Color-coded by Priority Tier (🔴 Tier 1, 🟠 Tier 2, 🟡 Tier 3, 🟢 Tier 4).
     - Hover metadata: Gene (`SYMBOL`), `Locus`, `HGVSc`, `Consequence`, `Priority Score`, `SpliceAI Δ`.
   - **Subplot 2: Top Prioritized Gene Burden Bar Chart**:
     - Ranks top 20 genes by Tier 1 and Tier 2 variant burden.

2. **Standalone Web Application Integration (`app.py`)**:
   - In [`src/web_app/app.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/web_app/app.py#L380):
     - Added an interactive **🔍 Select Gene for Detailed Transcript Inspection** dropdown selector in Tab 3.
     - Allows clinicians to view either the overall **Multi-Gene Cohort View (Manhattan)** or drill down into any specific gene (e.g. `MYBPC3`, `LMNA`, `TTN`) to render its isolated 5-tab report!

---

## 🧪 Verification & Unit Test Suite

- **Created Unit Test Suite**: [`tests/python/test_multigene_plotting.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/tests/python/test_multigene_plotting.py).
- **Execution Results**:
  - `test_multigene_figure_switch`: **PASSED** (Multi-gene datasets render Manhattan plot).
  - `test_singlegene_figure_rendering`: **PASSED** (Single-gene datasets render transcript map).
  - `test_multigene_report_generation`: **PASSED** (HTML dashboard generated cleanly).
- **Full Test Suite Status**: 100% Pass Rate across `test_compound_het.py`, `test_multisample_trio.py`, and `test_multigene_plotting.py`.
