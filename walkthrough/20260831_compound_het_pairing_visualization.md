# 🔗 Compound Heterozygous Variant Pairing & Visualization Walkthrough
**Date**: `2026-08-31 13:42:00 CEST`  
**Location**: `/data_lab_PGP/pipelines/annotation_pipeline_new/` (HPC) $\leftrightarrow$ `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/` (Local Sandbox)  

---

## 📌 Feature Overview

We introduced explicit **Compound Heterozygous Variant Pairing** to automatically identify, cross-link, and display complementary TRANS alleles within affected genes.

Clinicians can now instantly identify which specific paternal allele pairs with which maternal allele, alongside their genomic loci, HGVSc, HGVSp, and parental origin.

---

## 🛠️ Key Implementation Steps

1. **Variant Pairing Engine (`COMPOUND_HET_PAIR`)**:
   - Updated [`filter_variants.py`](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py#L1950) in `_mark_pairwise_trans()` to compute exact partner loci and annotations:
     - **Paternal Allele (`11:47360000 C>T`)**: `Paternal | Partner (Maternal): 11:47370000-T-C (c.1800+1T>C)`
     - **Maternal Allele (`11:47370000 T>C`)**: `Maternal | Partner (Paternal): 11:47360000-C-T (c.1500delC / p.Pro501fs)`
   - Added fallback formatting for unphased single-sample candidate variants (`Unphased Candidate | Partners (...)`).

2. **4-Tab Clinical HTML Report Enhancements**:
   - Updated [`generate_clinical_prioritization_report.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/python/generate_clinical_prioritization_report.py#L670):
     - Renders a prominent **🔗 Compound Heterozygous Allele Pairings Detected** callout box in Tab 1 (Executive Summary).
     - Includes `Compound Het Pair Details` column in the candidate variant table.

3. **Standalone Web Application (`app.py`)**:
   - Updated [`src/web_app/app.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/web_app/app.py#L305) with an interactive **🔗 Compound Heterozygous Allele Pairings Detected** expandable widget.

---

## 🧪 Verification Results

- **Unit Test Suite**: `tests/python/test_compound_het.py` passed with 100% success.
- **End-to-End Trio Pipeline**: `tests/python/test_multisample_trio.py` successfully generated [`MYBPC3_trio_clinical_prioritization_report.html`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/test_data/test_run/reports/MYBPC3_trio_clinical_prioritization_report.html) featuring the compound het pairings banner.
