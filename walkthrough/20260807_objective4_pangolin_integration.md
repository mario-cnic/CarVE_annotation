# Objective 4 Completion Walkthrough: Pangolin Splicing Predictor Integration

**Date**: 2026-08-07  
**Author**: Bioinformatics Pipeline Team  
**Status**: ✅ COMPLETED  

---

## 1. Executive Summary

All 5 tasks for **Objective 4: Pangolin Splicing Predictor Integration** (GPL-3.0, deep learning tissue-specific splice site usage predictor) have been completed and verified:

### Key Accomplishments:

1. **Environment Setup** (`/home/mruizp/conda_envs/pangolin_env`):
   - Created dedicated Python 3.12 environment containing `pangolin` (1.0.2 - tkzeng/Pangolin PyTorch ensemble model), `torch` (CPU wheel), `pyfaidx`, `pysam`, `gffutils`, `pyvcf3`, and `htslib` (`bgzip`, `tabix`).

2. **Local Annotation Database Construction**:
   - Built `pangolin_grch38.db` (918 MB) locally at `/home/mruizp/data_lab_PGP/shared/utils/pangolin_db/pangolin_grch38.db` using `create_db.py` from GTF annotations and GRCh38 FASTA reference without downloading prebuilt remote DBs.

3. **HPC Execution Script** ([annotate_pangolin_vars.sh](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/annotate_pangolin_vars.sh)):
   - Configured execution script with deep intronic window parameter `-d 10000` (20,000 bp window per variant) and logger integration.

4. **Downstream Parsing & Status Contract** ([filter_variants.py](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py)):
   - Implemented `parse_pangolin()` to extract:
     - `Pangolin_max_score`: Max delta score across splice site loss/gain predictions.
     - `Pangolin_heart_lv_score`: Cardiac tissue channel score.
     - `Pangolin_heart_aa_score`: Cardiac tissue channel score.
     - `Pangolin_status`: Enum status contract (`scored`, `not_covered`, `error`).

---

## 2. Verification & Audit Results

- Tested execution on `BAG3` gene VCF (`RUNS/predictors_050826/_tmp/BAG3.vcf`).
- Successfully evaluated variants with PyTorch ensemble models across `BAG3` sequences.
- Parsed delta scores and applied status contract (`scored` vs `not_covered`).
