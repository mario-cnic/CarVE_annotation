# Pipeline Setup & Initial Documentation Log

**Date**: 2026-08-05 17:41:42 CEST  
**Author**: Mario Ruiz - Data Lab PGP  
**Pipeline Location**: `/data_lab_PGP/pipelines/annotation_pipeline_new`

---

## 🎯 Objective

Initialize and deploy a new universal genomic variant annotation, filtering, and visualization pipeline in `/data_lab_PGP/pipelines/annotation_pipeline_new`. The pipeline inherits the core HPC annotation modules from `/data_lab_PGP/projects/enrichment_vars_220126` and extends them to support:
1. Multi-format variant inputs (cDNA transcript ENST/NM notation, hg19 or GRCh38 genomic coordinates, direct VCFs, Excel, TSV, CSV, Parquet).
2. Automated resolution and standardization to **GRCh38** genomic coordinates.
3. gnomAD v4.1, Ensembl VEP 111 (CADD, REVEL, AlphaMissense, SpliceAI), and SPiPv2.1 annotations.
4. Fully customizable command-line variant filtering thresholds.
5. Publication-ready static plots (PDF/PNG) and interactive HTML dashboards.

---

## 🛠️ Implementation Summary

### 1. File Structure Setup
- Created main pipeline directories: `src/hpc/`, `src/python/`, `src/R/`, `resources/`, `walkthrough/`.
- Copied core gene resources:
  - [resources/cardio_genes_loc.bed](../resources/cardio_genes_loc.bed)
  - [resources/gene_transcript_mapping.txt](../resources/gene_transcript_mapping.txt)

### 2. Universal Input Converter
- Implemented [src/python/universal_variant_converter.py](../src/python/universal_variant_converter.py):
  - Ingests cDNA notation (`c.123A>G`, `NM_000257.3:c.526C>T`, `ENST00000343260:c.100A>G`) using GeneBe.
  - Ingests hg19 / GRCh37 coordinates and converts them to GRCh38 via `pyliftover`.
  - Ingests GRCh38 locus columns or direct VCF files.
  - Outputs standardized GRCh38 VCF format.
- Updated [src/hpc/variant_converter.sh](../src/hpc/variant_converter.sh) to execute `universal_variant_converter.py`.

### 3. Customizable Filtering Engine
- Implemented [src/python/filter_and_summarize.py](../src/python/filter_and_summarize.py):
  - Command-line flags for `--max-af`, `--min-revel`, `--min-alphamissense`, `--min-spip`, `--min-cadd`, and `--consequences`.
  - Outputs filtered Parquet/TSV/Excel tables and `filtering_summary_metrics.tsv`.

### 4. Visualization & Reporting Engine
- Implemented static figure generator in [src/R/plot_annotation_results.R](../src/R/plot_annotation_results.R) generating PDF and PNG charts (consequence breakdown, REVEL vs AlphaMissense scatter, AF spectrum).
- Implemented interactive HTML report generator in [src/python/generate_interactive_report.py](../src/python/generate_interactive_report.py) creating standalone Plotly HTML dashboards.

### 5. Master Pipeline Orchestrator
- Implemented [main.sh](../main.sh) CLI orchestrator supporting all customizable parameters and managing chained SGE cluster job submission (`qsub`, `-hold_jid`).

---

## 🔬 Verification & Syntax Checks

- Validated Python scripts: `python3 -m py_compile src/python/*.py` (Passed).
- Validated R scripts: `Rscript -e "parse(file='src/R/plot_annotation_results.R')"` (Passed).
- Validated Bash scripts: `bash -n main.sh`, `bash -n src/hpc/*.sh` (Passed).

---

## 📌 Next Steps

- Execute test run on sample input datasets across cDNA, hg19 coordinate, and VCF input formats.
- Verify generated `.parsed.clean.pq`, filtered tables, and HTML reports.
