# 🌐 Standalone Clinical Reporting Web Application Walkthrough & Handover
**Date**: `2026-08-31 11:10:00 CEST`  
**Genome Assembly**: `GRCh38` (Ensembl v111 / GATK Bundle v0)  
**Location**: `/data_lab_PGP/pipelines/annotation_pipeline_new/` (HPC) $\leftrightarrow$ `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/` (Local Sandbox)  

---

## 📌 Executive Summary

We developed and integrated a **Standalone Clinical Reporting Web Application** ([run_web_app.sh](../run_web_app.sh) $\rightarrow$ [`src/web_app/app.py`](../src/web_app/app.py)) into the annotation pipeline. 

Clinicians, bioinformaticians, and researchers can now launch a local or server-side web interface, drag-and-drop any variant file (`.xlsx`, `.parquet` / `.pq`, `.csv`, `.tsv`, `.vcf`, `.vcf.gz`), and instantly run ACMG/ClinGen tiering, search/filter variants, preview an embedded interactive 4-tab clinical report, and download standalone reports (`.html`), Excel workbooks (`.xlsx`), or binary Parquet tables (`.pq`).

---

## 🛠️ Key Components & Architecture

### 1. Web Application Core (`src/web_app/app.py`)
- **Universal File Parser**: Ingests `.xlsx`, `.pq`, `.csv`, `.tsv`, and raw `.vcf` / `.vcf.gz` (via `vcf_parser_pysam`).
- **Automated Tiering Engine**: Automatically runs `build_newImpact()` and `build_priority_tier()` from [`filter_variants.py`](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py).
- **Executive KPI Dashboard**: Live summary cards for Total Variants, Tier 1 Critical Candidates, Tier 2 Strong Candidates, Tier 3 VUS, and Tier 4 Benign.
- **Interactive Plotly Visualizations**: Priority Tier Donut Chart & Functional Impact Distribution.
- **Searchable & Filterable Data Grid**: Instant filtering by Gene symbol, Priority Tier, Splicing Predictor scores, and free text search across HGVS / Locus coordinates.
- **Embedded Live HTML Clinical Dashboard**: Renders the complete 4-tab interactive clinical prioritization dashboard rendered via [`generate_clinical_prioritization_report.py`](../src/python/generate_clinical_prioritization_report.py).
- **Export & Download Center**: One-click downloads for `.html` (standalone dashboard), `.xlsx` (formatted Excel table), and `.pq` (compressed binary Parquet).

---

## 🚀 How to Launch and Use the Web App

### Terminal Command:
```bash
./run_web_app.sh [PORT]
```
*(Default Port: `8501`)*

### Accessing the Web App:
- **Local Access**: Open your browser at `http://localhost:8501`
- **Server Access**: Access `http://<server-ip>:8501` or set up an SSH port-forwarding tunnel:
  ```bash
  ssh -L 8501:localhost:8501 mruizp@samwise
  ```

---

## ✅ Verification & Test Results

1. **Dependency Installation**: `streamlit (v1.42.0)`, `openpyxl (v3.1.5)`, and `pyarrow (v19.0.1)` installed via `mamba` in the `datasci` Conda environment.
2. **Syntax & Compilation Verification**: Verified `python3 -m py_compile src/web_app/app.py` executed cleanly.
3. **Module Interoperability**: Successfully tested `filter_variants` and `generate_gene_report` imports inside the Streamlit application runtime.
