# Walkthrough: gnomAD v4.1 Root Cause Resolution & 7-Gene Cohort Reprocessing

**Date & Time**: 2026-08-20 11:00:00 CEST

## Objective
Investigate why `gnomADv4_AF_grpmax_joint` appeared missing or 0% populated in the 7-gene cohort (`BAG3`, `DSP`, `FLNC`, `LMNA`, `MYBPC3`, `PKP2`, `TTN`) despite previous VEP chromosome synonym fixes, resolve the root cause, re-parse clean datasets, and regenerate the clinical prioritization reports with the new interactive transcript visualization maps.

---

## 🔍 Root Cause Analysis & Empirical Discovery

1. **VCF Annotation Verification**:
   - Inspected `RUNS/run_20260813_1028/annotation/MYBPC3.annVEP.vcf.gz` using `cyvcf2` and `pysam`.
   - **Result**: VEP custom annotation **did** successfully annotate `gnomADv4` fields into the VCF `CSQ` tag (`124,370` non-empty entries for `gnomADv4_AF_grpmax_joint` in `MYBPC3.annVEP.vcf.gz`).

2. **Parsing Bottleneck Identification**:
   - In `src/hpc/vcf2parsed.sh`, Stage 1 invokes `vcf_parser_pysam.py` with `--vep_columns resources/all_but_old_gnomad_vep_cols.txt`.
   - Inspection of [resources/all_but_old_gnomad_vep_cols.txt](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/all_but_old_gnomad_vep_cols.txt) revealed that **`gnomADv4` field names were missing from the user-specified column list**.
   - Because `vcf_parser_pysam.py` strictly filters CSQ fields against `selected_vep_columns`, it omitted all `gnomADv4_*` fields during VCF-to-TSV conversion.
   - Stage 2 (`filter_variants.py`) subsequently filled missing `gnomADv4_AF_grpmax_joint` with NaN, causing 0% non-null values in the `.parsed.clean.pq` files.

---

## 🛠️ Fix Implemented

1. **Updated Column Registry**:
   - Added all `gnomADv4_*` field names (`gnomADv4_AF_joint`, `gnomADv4_AF_grpmax_joint`, `gnomADv4_grpmax_joint`, `gnomADv4_AF_exomes`, `gnomADv4_AF_genomes`, etc.) to [resources/all_but_old_gnomad_vep_cols.txt](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/all_but_old_gnomad_vep_cols.txt).

2. **Reprocessed 7-Gene Cohort Datasets**:
   - Re-parsed annotated VCFs to `.parsed.clean.pq` using `vcf2parsed.sh` for all 7 genes in `RUNS/run_20260813_1028/`.

---

## 📊 Empirical Verification Results

| Gene | Total Variants | Populated `gnomADv4_AF_grpmax_joint` | Populated % | Status |
| :--- | :--- | :--- | :--- | :--- |
| **BAG3** | 8,089 | 3,778 | **46.70%** | ✅ Verified |
| **DSP** | 22,168 | 10,427 | **47.04%** | ✅ Verified |
| **FLNC** | 20,437 | 9,467 | **46.32%** | ✅ Verified |
| **LMNA** | 93,450 | 43,231 | **46.26%** | ✅ Verified |
| **MYBPC3** | 28,207 | 13,585 | **48.16%** | ✅ Verified |
| **PKP2** | 47,101 | 22,765 | **48.33%** | ✅ Verified |
| **TTN** | 485,122 | 228,020 | **47.00%** | ✅ Verified |

---

## 📈 Clinical Prioritization HTML Reports

Regenerated interactive 4-tab clinical prioritization HTML reports for all 7 genes incorporating:
1. **Multi-Track SVG Transcript Visualizer & Variant Position Map** (General, Splicing, and Missense tracks).
2. **gnomAD v4.1 Joint PopMax Allele Frequencies & Population Ancestry Tags** (`gnomADv4_AF_grpmax_joint`).

Reports generated:
- [BAG3_clinical_prioritization_report.html](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_20260813_1028/reports/BAG3_clinical_prioritization_report.html)
- [DSP_clinical_prioritization_report.html](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_20260813_1028/reports/DSP_clinical_prioritization_report.html)
- [FLNC_clinical_prioritization_report.html](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_20260813_1028/reports/FLNC_clinical_prioritization_report.html)
- [LMNA_clinical_prioritization_report.html](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_20260813_1028/reports/LMNA_clinical_prioritization_report.html)
- [MYBPC3_clinical_prioritization_report.html](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_20260813_1028/reports/MYBPC3_clinical_prioritization_report.html)
- [PKP2_clinical_prioritization_report.html](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_20260813_1028/reports/PKP2_clinical_prioritization_report.html)
- [TTN_clinical_prioritization_report.html](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_20260813_1028/reports/TTN_clinical_prioritization_report.html)

---

## 🦎 6. Pangolin Engine Audit & Multi-Level Fixes

1. **Chromosome Name Resolution (`pangolin.py`)**:
   - `pangolin_grch38.db` stores chromosome names without `chr` prefix (`11`, `1`, `X`).
   - GRCh38 VCF inputs use `chr11`. `gffutils.region(('chr11', ...))` returned empty sets, causing Pangolin to skip 100% of variants.
   - Updated `get_genes()` in [pangolin.py](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/external/Pangolin-main/pangolin/pangolin.py#L70-L82) to try both `chr` and non-`chr` prefixes.

2. **Column Registry Update**:
   - Added raw INFO tag `Pangolin` to [resources/all_but_old_gnomad_vep_cols.txt](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/all_but_old_gnomad_vep_cols.txt#L463) to allow `vcf_parser_pysam.py` to retain raw Pangolin predictions for downstream extraction.

3. **Environment & PATH Resolution**:
   - Updated [annotate_pangolin_vars.sh](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/annotate_pangolin_vars.sh) to resolve `PANGOLIN_ENV=/home/mruizp/conda_envs/pangolin_env` and export `$PANGOLIN_ENV/bin` into `$PATH`.

---

## 🧪 7. Self-Contained Test Suite & `--test` Flag

1. **Test Dataset**:
   - Created lightweight 50-variant test VCF: `test_data/raw_vcfs/MYBPC3_test.vcf.gz` (3.2 MB total folder size).
2. **`--test` Flag Integration**:
   - Updated `main.sh` to accept `--test`. Routes input to `test_data/raw_vcfs/`, outputs to `test_data/test_run/`, and forces `--overwrite-all`.

---

## ⚡ 8. Quality Auditor Subshell Pipe Optimization (`audit_run_results.py`)

- **Subshell Deadlock Resolution**: Replaced slow, subshell-piped `zcat file | grep | wc` calls in `audit_run_results.py` with fast indexed variant counting via `bcftools index -n`.
- **Performance Impact**: Reduced audit query time per VCF from 30+ seconds down to **0.002 seconds**, completing full run audits in **0.40 seconds**.

