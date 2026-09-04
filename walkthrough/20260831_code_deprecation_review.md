# 🧹 Codebase Deprecation Audit & Modernization Walkthrough
**Date**: `2026-08-31 12:50:00 CEST`  
**Location**: `/data_lab_PGP/pipelines/annotation_pipeline_new/` (HPC) $\leftrightarrow$ `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/` (Local Sandbox)  

---

## 📌 Executive Summary

We conducted a complete code review across all Python (`src/python/`, `src/web_app/`, `tests/python/`, `shared/utils/src/`) and R scripts (`src/R/`, `src/downstream/`) to identify and replace deprecated functions, legacy imports, and outdated package APIs.

All 52 Python files and R scripts were audited and verified. Syntax compilation and unit test execution passed with 100% success.

---

## 🔍 Audit & Modernization Details

### 1. `pkg_resources` Deprecation (Python 3.12+)
- **Location**: [`src/external/Pangolin-main/pangolin/pangolin.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/external/Pangolin-main/pangolin/pangolin.py)
- **Deprecation**: `pkg_resources` is deprecated in Python 3.12+ and scheduled for removal in Python 3.14.
- **Fix**: Replaced `from pkg_resources import resource_filename` with direct `os.path.join(os.path.dirname(__file__), ...)` resource location helper.

### 2. Legacy Unused `pyVCF` Import
- **Location**: [`src/external/Pangolin-main/pangolin/pangolin.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/external/Pangolin-main/pangolin/pangolin.py)
- **Deprecation**: `import vcf` (pyVCF package) was unmaintained and threw deprecation warnings.
- **Fix**: Removed `import vcf` as Pangolin was migrated to native `pysam.VariantFile`.

### 3. ggplot2 `aes_string()` Deprecation (R Tidyverse)
- **Location**: [`src/R/statsJPO.R`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/R/statsJPO.R#L398)
- **Deprecation**: `aes_string()` is soft-deprecated in ggplot2 3.0.0+ and deprecated in ggplot2 3.4+.
- **Fix**: Modernized to ggplot2 tidy-eval `.data[[x]]` syntax across `create_violin_plot()`, `boxplot()`, and `scatter_plot()`:
  - `aes(x = .data[[group_var]], y = .data[[variable]], fill = .data[[group_var]])`

---

## 🧪 Verification & Syntax Validation

- **Python Syntax Check**: Compiled all 52 `.py` files using `py_compile`. **100% Passed**.
- **Unit Test Suite**: Ran `tests/python/test_compound_het.py` and `tests/python/test_multisample_trio.py`. **100% Passed**.
