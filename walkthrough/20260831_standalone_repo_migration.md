# 🚀 Standalone Web App Repository Migration & Fix Walkthrough
**Date**: `2026-08-31 14:06:00 CEST`  
**Target Repository**: `/data_lab_PGP/pipelines/clinical_variant_prioritization/` (HPC) $\leftrightarrow$ `/home/mruizp/data_lab_PGP/pipelines/clinical_variant_prioritization/` (Local Sandbox)  

---

## 📌 Summary of Bug Fix & Active Launch

1. **Demo Dataset File Stream Bug Fix**:
   - **Root Cause**: `parse_uploaded_file()` in `app.py` originally assumed all inputs were Streamlit `UploadedFile` objects equipped with `.getvalue()`. When opening local demo files via standard Python file handles (`open("demo_data/MYBPC3_trio_test.vcf", "rb")`), `.getvalue()` threw an AttributeError, causing `df_raw` to remain `None`.
   - **Fix Implemented**: Updated [`app.py`](file:///home/mruizp/data_lab_PGP/pipelines/clinical_variant_prioritization/app.py#L60) to check `hasattr(uploaded_file, "getvalue")` and fallback to `.read()` for standard file streams:
     ```python
     file_bytes = uploaded_file.getvalue() if hasattr(uploaded_file, "getvalue") else uploaded_file.read()
     ```

2. **Standalone App Launch**:
   - Launched the new application from `/home/mruizp/data_lab_PGP/pipelines/clinical_variant_prioritization/`.
   - **Active URL**: `http://localhost:8502`

---

## 📁 Repository Structure (`clinical_variant_prioritization`)

```text
/home/mruizp/data_lab_PGP/pipelines/clinical_variant_prioritization/
├── app.py                         # Standalone Streamlit Web Application
├── run_app.sh                     # Executable launcher script
├── requirements.txt               # Lightweight Python dependencies
├── environment.yml                # Conda/Mamba environment setup
├── README.md                      # Comprehensive standalone documentation
├── .gitignore                     # Git ignore rules
├── src/
│   ├── filter_variants.py         # Standalone 4-Tier classification & Pedigree engine
│   └── report_generator.py        # Standalone HTML dashboard & Plotly charts engine
└── demo_data/
    ├── MYBPC3_trio_test.vcf       # Demo multi-sample trio VCF
    └── MYBPC3_trio.ped            # Demo pedigree file
```
