# 5' UTR Annotator Resolution & Clinical Dashboard Enhancements

**Date**: 2026-08-18  
**Author**: Antigravity  
**Target Genome Build**: GRCh38 / Ensembl 111  

---

## 1. Executive Summary & Root Cause

During quality review of the multi-tab Clinical Prioritization Dashboard, two primary questions were investigated:
1. **Header Alignment & ClinVar Visibility**: Column headers in the Splicing tab were horizontally shifted by two columns due to a hardcoded HTML `<thead>` mismatch, and ClinVar classifications appeared as `"-"` because the script looked exclusively for `CLNSIG` instead of the pipeline's `clinvar_clnsig` / `CLIN_SIG` columns.
2. **Missing 5' UTR Annotator Fields**: Despite the VEP job executing with `--plugin UTRAnnotator,file=.../uORF_5UTR_GRCh38_PUBLIC.txt`, the columns `5UTR_annotation` and `5UTR_consequence` were absent from downstream Parquet tables.

### Root Cause of Missing 5' UTR Annotator:
* In the VCF parser configuration file [`resources/all_but_old_gnomad_vep_cols.txt`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/all_but_old_gnomad_vep_cols.txt), the column names were formatted with a leading underscore:
  ```text
  _5UTR_annotation
  _5UTR_consequence
  ```
  *(likely inherited from R's default behavior of prefixing numeric column names with `_` or `X`)*.
* When [`vcf_parser_pysam.py`](file:///home/mruizp/data_lab_PGP/shared/utils/src/vcf_parser_pysam.py) parsed the VCF's CSQ metadata (`5UTR_annotation` and `5UTR_consequence` without leading underscore), exact regex matching failed and logged:
  ```text
  WARNING: _5UTR_annotation not found in vep_columns
  ```
  resulting in the silent omission of all 5' UTR Annotator fields during the VCF $\rightarrow$ TSV stage.

---

## 2. Solutions Implemented

### A. Dynamic Table Header Rendering & Unified ClinVar Extractor
* **File Modified**: [`src/python/generate_clinical_prioritization_report.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/python/generate_clinical_prioritization_report.py)
* **HTML `<thead>` Generation**: Now rendered dynamically directly from the column definitions list (`render_table_html`), ensuring 100% mathematical alignment between table headers and cell values across all 4 tabs.
* **Unified ClinVar Extractor**: Automatically cascades through `clinvar_clnsig` $\rightarrow$ `CLIN_SIG` $\rightarrow$ `CLNSIG` and applies color-coded clinical significance badges.

### B. Column Dictionary & Regex Matching Fix
* **Files Modified**:
  1. [`resources/all_but_old_gnomad_vep_cols.txt`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/all_but_old_gnomad_vep_cols.txt): Corrected entries to `5UTR_annotation` and `5UTR_consequence`.
  2. [`/home/mruizp/data_lab_PGP/shared/utils/data/vcf_all_columns.txt`](file:///home/mruizp/data_lab_PGP/shared/utils/data/vcf_all_columns.txt): Removed deprecated `_5UTR_*` aliases.
  3. [`/home/mruizp/data_lab_PGP/shared/utils/src/vcf_parser_pysam.py`](file:///home/mruizp/data_lab_PGP/shared/utils/src/vcf_parser_pysam.py): Made column matching resilient to leading underscores:
     ```python
     clean_col = user_col.lstrip("_")
     pattern = re.compile(rf"^(_)?({re.escape(clean_col)}|{re.escape(user_col)})(\.\d+)?$")
     ```

---

## 3. Verification & Results

1. **VCF Parser Validation**:
   * Tested on `BAG3.annotated.vcf.gz`.
   * Successfully extracted **41 variants** with upstream Open Reading Frame consequences (e.g. `5_prime_UTR_premature_start_codon_gain_variant`, `uORF_frameshift`).
2. **Prioritization Score & Tiering**:
   * The UTR translation disruption rule in [`filter_variants.py`](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py) activated, assigning `NEW_IMPACT: HIGH`, `+20.0 pts` score bonus, and promoting severe uORF alterations to **Tier 2**.
3. **Dashboard Disclaimer Verification**:
   * Regenerated [`BAG3_clinical_prioritization_report.html`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_20260813_1028/reports/BAG3_clinical_prioritization_report.html).
   * `5UTR_consequence` warning disclaimer cleared; UTR tab now displays translation-modifying variants with complete Kozak and uORF context.
