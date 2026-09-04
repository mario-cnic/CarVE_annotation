# 📊 Population Allele Frequency Scale & Help Tab LaTeX Formatting Walkthrough
**Date**: `2026-08-31 13:55:00 CEST`  
**Location**: `/data_lab_PGP/pipelines/annotation_pipeline_new/` (HPC) $\leftrightarrow$ `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/` (Local Sandbox)  

---

## 📌 Summary of Changes

### 1. Population Allele Frequency Scoring Scale Update
In [`filter_variants.py`](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py#L1115), population allele frequency (gnomAD v4.1 PopMax AF) scoring was updated to strictly enforce 3 distinct frequency brackets:
- **Ultra-Rare ($\text{AF} < 1 \times 10^{-4}$ or $0.0001$)**: **Positive Score Bonus (+15.0 pts)**.
- **Moderate Rarity ($1 \times 10^{-4} \le \text{AF} \le 0.01$ or $1\%$)**: **Neutral Score (0.0 pts)**.
- **Common Variant ($\text{AF} > 0.01$ or $1\%$)**: **Heavy Penalty (-30.0 pts)**.

### 2. Help Tab LaTeX & Math Formatting Fixes
Updated KaTeX and HTML math formatting across both [`generate_clinical_prioritization_report.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/python/generate_clinical_prioritization_report.py#L1485) and [`src/web_app/app.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/web_app/app.py#L480):
- Fixed backslash escaping for KaTeX symbols: `\\ge`, `\\le`, `\\Delta`, `\\text{AF}`.
- Replaced ambiguous text with explicit mathematical formulas: $\text{AF} < 1 \times 10^{-4}$.

---

## 🧪 Verification & Results

- **Entire Test Suite**: `test_compound_het.py`, `test_multisample_trio.py`, and `test_multigene_plotting.py` executed with **100% Pass Rate**.
- **Trio Variant Re-scoring Result**:
  - `11:47328406 C>A` (De Novo, AF = 0.0): $+15\text{ pts (AF)} + 20\text{ pts (De Novo)} = 35.0\text{ pts}$ (Tier 3).
  - `11:47350000 G>A` (Recessive Hom): $+15\text{ pts (AF)} + 35\text{ pts (LoF)} + 15\text{ pts (Recessive)} = 65.0\text{ pts}$ (Tier 2).
