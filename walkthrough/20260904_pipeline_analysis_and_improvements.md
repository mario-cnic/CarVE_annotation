# 🔍 Universal Variant Annotation Pipeline: Analysis & Improvement Roadmap

**Date & Time**: 2026-09-04 15:52:00 CEST  
**Target Environment**: GRCh38 Pipeline (`/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new`)  
**Author**: Bioinformatics Specialist  
**Provenance / Tools**: Ensembl VEP 111, CADD v1.6, REVEL, AlphaMissense, SPiPv2.1, Pangolin, Local SpliceAI (`-D 10000`), Branchpointer, LaBranchoR, SpliceVault, Pytest 9.1.1  

---

## Executive Summary

A comprehensive architectural, computational, and biological audit of the **Universal Genomic Variant Annotation Pipeline** (`annotation_pipeline_new`) was conducted. The pipeline represents a production-grade, state-of-the-art framework for genomic variant standardization (GRCh38), deep-learning splicing ensemble prediction, and interactive clinical reporting. 

However, our empirical analysis uncovered key structural vulnerabilities, broken shared library dependencies, unfulfilled ACMG evidence engines, hardcoded HPC path dependencies, and performance bottlenecks that impact computational scalability and local/HPC portability.

---

## 🔑 Key Analytical Findings & Vulnerability Matrix

| Category | Component / File | Current Issue / Vulnerability | Impact | Proposed Improvement |
| --- | --- | --- | --- | --- |
| **Unit Test & Shared Library** | `tests/python/test_*.py` <br> `filter_variants.py` | `ImportError: cannot import name 'get_full_gene_curation_dataframe'` from `modules.disease_hpo`. | ❌ Python unit test suite breaks on collection (3 module errors). | Fix signature import in `filter_variants.py` / `disease_hpo.py` and update test mock fixtures. |
| **HPC Path Portability** | `src/hpc/*.sh` <br> `src/python/*.py` | 54 hardcoded `/data_lab_PGP/` paths instead of dynamic resolution or environment-based fallback (`/home/mruizp/data_lab_PGP/`). | ⚠️ Causes `No such file or directory` or execution failure when run in local sandbox or alternate compute nodes. | Implement dynamic path translation helper in bash wrappers (`${DATA_LAB_PGP:-/home/mruizp/data_lab_PGP}`). |
| **ACMG In-Silico Prioritization** | `TODO.md` (Phase 2) <br> `vcf_parser_pysam.py` | Missing automated ACMG/AMP evidence code engine (`PVS1`, `PM2_Supporting`, `PP3`, `BP4`, `BA1/BS1`). | ⚠️ Prioritization relies on manual score interpretation rather than standardized ClinGen/ACMG criteria. | Implement deterministic ACMG rule engine producing `ACMG_Evidence_Codes` and `ACMG_Suggested_Tier` columns. |
| **High-Throughput Parallelization** | `annotate_spip_vars.sh` <br> `annotate_pangolin_vars.sh` <br> `annotate_vep_vars.sh` | SpliceAI uses parallel VCF chunking ($\ge 15,000$ variants), but SPiP, Pangolin, and VEP process full VCFs sequentially per gene. | 🐢 Computational bottleneck on multi-thousand variant callsets (WES/WGS). | Extend streaming parallel VCF chunking (`split_vcf_chunks.py`) to Pangolin, SPiP, and VEP inference jobs. |
| **Cohort Master Reporting** | `generate_clinical_prioritization_report.py` | Reports are strictly single-gene focused (`<GENE>_clinical_prioritization_report.html`). | ⚠️ No cross-gene overview or cohort-wide variant tier distribution dashboard across target panel (e.g. 73 cardio genes). | Build `src/python/generate_cohort_master_dashboard.py` aggregating top-tier variants across all panel genes into 1 interactive dashboard. |
| **ClinGen & Cardiac Constraints** | `resources/` <br> `generate_clinical_prioritization_report.py` | Missing ClinGen Haploinsufficiency (`HI`) / Triplosensitivity (`TS`) headers, gnomAD v4 gene constraint metrics (`LOEUF`, `pLI`), and cardiac domain mappings (`BAG3`, `LMNA`, `TTNtv` A-band PSI flag). | 🧬 Clinicians lack contextual domain & isoform constraint flags for variant evaluation. | Integrate ClinGen dosage & gnomAD constraint databases and domain annotation overlays into report headers. |
| **Containerization & Deployment** | `resources/containers/` | Conda environments are fragmented (`datasci`, `genomics`, `vcf_parser`, `spip_env`, `spliceai_env`, `pangolin_env`). | 📦 Multi-env switching creates sheath shebang and ABI mismatch friction across cluster nodes. | Finalize single immutable Apptainer `.sif` container (`annotation_suite.sif`) encapsulating all environments. |

---

## 🛠️ Detailed Architectural Improvement Plan

### 1. Priority 1: Fix Shared Integration & Unit Test Suite
- Resolve `ImportError` in `filter_variants.py` by aligning function signatures with `modules.disease_hpo.py`.
- Sanitize all 54 hardcoded `/data_lab_PGP/` paths in `src/hpc/*.sh` and `src/python/*.py` using relative environment resolution:
  ```bash
  BASE_PGP="${DATA_LAB_PGP:-/home/mruizp/data_lab_PGP}"
  ```
- Ensure 100% clean test execution on `bash tests/run_all_tests.sh`.

### 2. Priority 2: Automated ACMG/AMP Evidence Engine (Phase 2)
Implement rule-based classification in `vcf2parsed.sh` / `filter_variants.py`:
- **`PVS1`**: Null LoF variant (`stop_gained`, `frameshift_variant`, `splice_acceptor`, `splice_donor`) in a gene intolerant to loss-of-function (`pLI >= 0.90` or `LOEUF < 0.35`).
- **`PM2_Supporting`**: Absent or extremely rare in gnomAD v4 joint population (`AF_grpmax_joint < 0.0001` or missing).
- **`PP3`**: Computational evidence supports a deleterious effect on gene/protein (`SpliceAI > 0.50` or `Pangolin > 0.50` OR `AlphaMissense > 0.56` + `REVEL > 0.70`).
- **`BP4`**: Multiple in-silico algorithms predict benign impact (`AlphaMissense < 0.34` and `REVEL < 0.40` and `SpliceAI < 0.10`).

### 3. Priority 3: Cohort Master Dashboard (`cohort_master_dashboard.html`)
- Aggregate individual gene Parquet tables (`results/*.parsed.clean.pq`) into a unified cross-gene dataset.
- Filter and surface Class 1 (High Priority / Pathogenic) and Class 2 (Moderate Priority / VUS) variants across the entire patient cohort.
- Provide interactive search, filtering by gene panel, allele frequency, and pathogenicity tier, with 1-click Excel export.

### 4. Priority 4: Universal Predictor Chunking & GPU Acceleration
- Extend `src/python/split_vcf_chunks.py` to enable streaming chunking for Pangolin, SPiP, and VEP when input exceeds 10,000 variants.
- Add optional GPU submission options (`qsub -l gpu=1`) in HPC wrapper scripts for PyTorch/TensorFlow deep-learning inference.

---

## Provenance & Verification Protocol

1. **Automated Verification**:
   - `bash tests/run_all_tests.sh`
   - `python3 src/python/audit_run_results.py --test`
2. **Genome Build Compliance**:
   - Standardized strictly on **GRCh38** (Ensembl release 111).
