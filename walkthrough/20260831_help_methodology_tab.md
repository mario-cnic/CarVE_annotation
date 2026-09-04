# ❓ Pipeline Help & Algorithm Methodology Tab Walkthrough
**Date**: `2026-08-31 13:45:00 CEST`  
**Location**: `/data_lab_PGP/pipelines/annotation_pipeline_new/` (HPC) $\leftrightarrow$ `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/` (Local Sandbox)  

---

## 📌 Executive Summary

We created a comprehensive **❓ Help & Algorithm Methodology** tab integrated into both:
1. **The HTML Clinical Prioritization Dashboard** ([`MYBPC3_trio_clinical_prioritization_report.html`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/test_data/test_run/reports/MYBPC3_trio_clinical_prioritization_report.html)) $\rightarrow$ 5th Navigation Tab (`tab-help`).
2. **The Standalone Streamlit Web Application** ([`src/web_app/app.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/web_app/app.py)) $\rightarrow$ 5th Navigation Tab (`tab_help`).

---

## 📖 Content Included in the Help Tab

### 1. 🏆 4-Tier Clinical Stratification Framework Table
- **Tier 1 (🔴 Critical Pathogenic Candidate)**: ClinVar Pathogenic / Likely Pathogenic OR Loss-of-Function OR Priority Score $\ge 60.0$ with gnomAD AF $\le 0.001$.
- **Tier 2 (🟠 Likely Deleterious / Strong Candidate)**: In silico Priority Score $\ge 35.0$ OR confirmed **De Novo**, **Autosomal Recessive (Hom)**, or **Compound Heterozygous (TRANS)** variants.
- **Tier 3 (🟡 VUS / Moderate Potential)**: Priority Score $\ge 15.0$ OR Splicing Delta $\ge 0.20$ floor (SpliceAI 20kb, Pangolin, SPiP).
- **Tier 4 (🟢 Benign / Tolerated)**: Common Population Variants (gnomAD v4.1 AF $> 0.01$) or ClinVar Benign.

### 2. 📊 Multi-Evidence Priority Scoring Formula
- **ClinVar Classification Weights**: Pathogenic (+40 pts), Likely Pathogenic (+30 pts), Benign (-30 pts).
- **Population Frequency**: Rare (AF $\le 0.0001$, +15 pts), Common (AF $> 0.01$, -25 pts).
- **In Silico Missense**: AlphaMissense ($0-1.0 \times 25$ pts), REVEL ($0-1.0 \times 25$ pts), CADD Phred ($\text{Score} / 2.0$).
- **Deep Splicing Predictors**: SpliceAI 20kb custom context ($\text{Max } \Delta \times 30$ pts), Pangolin Heart LV ($\text{Score} \times 25$ pts), SPiP ($\text{Score} \times 25$ pts), SpliceVault (+10 pts).
- **Pedigree Inheritance Bonuses**: De Novo (+20 pts), Recessive / Compound Het TRANS (+15 pts), Shared Affected (+10 pts).

### 3. 🧬 Pedigree Inheritance & TRANS Phasing Definitions
- Formal definitions of **De Novo**, **Autosomal Recessive**, **Compound Heterozygous (TRANS vs CIS)**, and **Unphased Candidates**.

---

## 🧪 Verification & Output

- **Report Dashboard**: Regenerated [`MYBPC3_trio_clinical_prioritization_report.html`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/test_data/test_run/reports/MYBPC3_trio_clinical_prioritization_report.html) with 5 fully interactive tabs.
- **Streamlit Web App**: Accessible live under tab 5 (`❓ Pipeline Help & Methodology`).
