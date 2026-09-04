# 📋 Master Annotation Pipeline & Clinical Dashboard Roadmap

---

## 🚀 Phase 1: High-Performance Compute & Orchestration (Inference & Checkpointing)
- [x] **State-Aware Pipeline Checkpointing & Granular Flags**
  - [x] Implement `-s "$file"` file integrity checks across all pipeline stages (Steps 1–5)
  - [x] Add dynamic SGE dependency hold chaining (`active_predictor_jobs`)
  - [x] Add granular CLI override flags: `--overwrite-all`, `--force-spliceai`, `--force-vep`, `--force-pangolin`, `--force-spip`, `--force-branchpoint`, `--force-merge`, `--force-vcf2parsed`, `--force-reports`, `--skip-genes`
- [x] **Parallel SpliceAI VCF Chunking Engine**
  - [x] Implement streaming header-preserving VCF chunker (`src/python/split_vcf_chunks.py`)
  - [x] Integrate multi-worker parallel inference pool matching `$NSLOTS` in `src/hpc/annotate_spliceai_vars.sh`
  - [x] Add automated volume thresholding ($\ge 15,000$ variants) and individual chunk checkpointing
  - [x] Lossless `bcftools concat -a` merge and `tabix` indexing
  - [x] Unit test suite in `tests/python/test_vcf_chunking.py` (100% passed)
- [x] **Built-in Run Quality & Completeness Auditor**
  - [x] Implement `src/python/audit_run_results.py` with predictor completeness %, deep learning throughput, and schema diffing
  - [x] Integrate `--audit-run` and `--audit-only` CLI flags into `main.sh`
  - [x] Add automated Step 6 post-run SGE audit job (`audit_${RUN_NAME}`) chained to terminal tasks

---

## 🧬 Phase 2: Biological & Clinical Variant Prioritization
- [x] **gnomAD v4 Joint Population Frequency Integration**
  - [x] Fix UCSC `chr10` $\leftrightarrow$ Ensembl `10` contig mapping via `--synonyms chr_synonyms.txt` in `annotate_vep_vars.sh`
  - [x] Add `gnomADv4 AF grpmax` across all 4 dashboard tabs with fallback detection and warning alerts
- [x] **5' UTR Annotator & Multi-Tab Dashboard**
  - [x] Integrate 5' UTR consequence, uORF, and start/stop disruption annotations into `vcf_parser_pysam.py`
  - [x] Implement 4-tab Clinical Prioritization Report (`generate_clinical_prioritization_report.py`)
- [x] **Centralized Configuration & Shared Library Fallbacks**
  - [x] Implement `config/env.sh` and `src/python/config.py` path resolver
  - [x] In-memory fallbacks for deprecated shared utilities (`get_full_gene_curation_dataframe`, `filter_by_custom_gene_list`) without editing external shared files
  - [x] ACMG calculation externalized to downstream app pipeline; lightweight annotations retained in-pipeline
- [x] **Cardiomyopathy Domain & Cardiac Isoform / PSI Filter**
  - [x] Add cardiac ventricle Percent Spliced In (PSI > 85%) flag for `TTN` truncating variants (TTNtv in A-band)
  - [x] Map critical DCM structural domains (`BAG3` BAG domain, `LMNA` rod domain, `FLNC` Ig-like folds)

---

## 📊 Phase 3: Reporting & Clinical Dashboard Enhancements
- [x] **Multi-Gene Cohort Master Dashboard (`generate_cohort_master_dashboard.py`)**
  - [x] Aggregate top-tier prioritized variants across all single-gene Parquet tables into a unified master cohort dashboard
  - [x] Cohort-level summary metrics (distribution of High/Moderate candidates, gene breakdown, splicing vs missense burden)
  - [x] Interactive Plotly visualizations and 1-click CSV candidate exporter
  - [x] Comprehensive pytest coverage in `tests/python/test_generate_cohort_master_dashboard.py`

---

## 🔒 Phase 4: Containerization & Infrastructure Stability
- [x] **Unified Apptainer / Singularity SIF Container (`annotation_pipeline.sif`)**
  - [x] Create definition file `resources/containers/annotation_pipeline.def` packaging Python datasci, R, genomics binaries, and Plotly
  - [x] Automated build script `src/hpc/build_container.sh` with environment-resolved output paths
  - [x] Integrated `--use-container` and `--sif <PATH>` CLI flags in `main.sh` and execution wrapper in `config/env.sh`

