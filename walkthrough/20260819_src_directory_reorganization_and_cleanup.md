# 🗂️ Source Directory Reorganization & Obsolete Script Archival

**Date & Time**: `2026-08-19 10:47:00 CEST`  
**Target Environment**: Pipeline Source Space (`src/`) & Local Sandbox  
**Scope**: Codebase restructuring, modularization, and safe archival of legacy scripts.

---

## 🎯 1. Overview & Rationale

As the annotation pipeline grew with advanced functional predictors (VEP 111, CADD, REVEL, AlphaMissense, SPiP, Pangolin, SpliceAI `-D 10000`, Branchpointer, LaBranchoR, SpliceVault, 5' UTR Annotator), the `src/` directory accumulated one-off scripts, deprecated wrappers, and legacy TSV processors alongside core production tools.

To enhance maintainability, clarity, and prevent execution errors, we performed a structured reorganization of all scripts into clean, functional domains while safely moving obsolete files to `./trash/src_archive/` (strictly adhering to data preservation bans on `rm`).

---

## 📁 2. Reorganized Modular Directory Architecture

```
src/
├── python/                     # Core active production modules (10 scripts)
│   ├── universal_variant_converter.py        # Step 1: Standardization & cDNA conversion
│   ├── annotate_spliceai.py                  # Step 2: SpliceAI deep learning inference
│   ├── split_vcf_chunks.py                   # Step 2: Parallel SpliceAI VCF chunking engine
│   ├── annotate_branchpointer.py             # Step 2: Branchpointer + LaBranchoR runner
│   ├── filter_and_summarize.py               # Step 5: CLI variant filtering module
│   ├── query_new_transcripts.py              # Core: Ensembl/RefSeq transcript resolver
│   ├── generate_interactive_report.py        # Step 5: Single-gene Plotly interactive dashboard
│   ├── generate_clinical_prioritization_report.py # Step 5: 4-tab clinical tiering dashboard
│   ├── audit_run_results.py                  # Step 6: Automated single-run quality auditor
│   └── track_pipeline_progress.py            # Master: Multi-pipeline execution tracker
│
├── hpc/                        # Core SGE cluster drivers & wrappers (11 scripts)
│   ├── variant_converter.sh                  # Step 1 SGE driver
│   ├── annotate_vep_vars.sh                  # Step 2 VEP SGE driver
│   ├── annotate_spip_vars.sh                 # Step 2 SPiP SGE driver
│   ├── annotate_pangolin_vars.sh             # Step 2 Pangolin SGE driver
│   ├── annotate_spliceai_vars.sh             # Step 2 Parallel SpliceAI driver
│   ├── annotate_branchpointer_vars.sh        # Step 2 Branchpointer driver
│   ├── merge_vep_spip.sh                     # Step 3 Dynamic VCF merge driver
│   ├── vcf2parsed.sh                         # Step 4 VCF-to-Parquet conversion driver
│   ├── gene_coords.sh                        # Genomic interval resolver
│   ├── logger.sh                             # HPC standard execution logger
│   └── build_container.sh                    # SIF Apptainer container builder
│
├── R/                          # Core plotting & statistical libraries (2 scripts)
│   ├── plot_annotation_results.R             # Step 5: ggplot2 static publication figures
│   └── statsJPO.R                            # Master statistical methods library
│
├── downstream/                 # Cohort burden & association analysis (12 scripts)
│   ├── 01_build_master.R & .sh               # Master cohort dataset builder
│   ├── 02_cohort_enrich.R & .sh              # Case-control cohort enrichment
│   ├── 03_burden.R & .sh                     # Gene-level variant burden testing
│   ├── 04_forest_plot.R & .sh                # Odds ratio forest plotting
│   ├── 05_missense_predictors.R              # Missense pathogenicity comparison
│   ├── 06_missense_spatial_burden.R          # Linear domain & spatial clustering
│   ├── 07_domain_burden.R                    # Functional domain burden analysis
│   └── 08_describe_ttn_results.R             # TTN-specific exon/domain PSI audit
│
├── tools/                      # Reference builders, backfills, & utilities (14 scripts)
│   ├── backfill_variant_prioritization.py    # Backfill PRIORITY_TIER / NEW_IMPACT
│   ├── backfill_spliceai_custom_columns.py   # Backfill custom SpliceAI columns
│   ├── build_labranchor_grch38.py            # LaBranchoR GRCh38 bed generator
│   ├── build_labranchor_ism_grch38.py        # LaBranchoR in-silico mutagenesis db
│   ├── compare_spliceai_predictions.py       # SpliceAI local vs VEP validator
│   ├── combine_valencia_hic.py & .sh         # Multi-center cohort combiner
│   ├── parse_raw_hic.py                      # Hospital raw table parser
│   ├── eda_parsed_pq.py                      # Parquet exploratory metrics
│   ├── noonan_string_audit.py                # Noonan syndrome PPI network audit
│   ├── pre_run_readiness_audit.py            # Pre-run file/resource auditor
│   ├── split_input_by_gene.py                # Raw table gene splitter
│   ├── store_tsv_as_parquet.py               # TSV-to-Parquet fast converter
│   ├── update_mappings_to_grch38.py          # Gene transcript mapping updater
│   └── validate_output.py                    # Output integrity verification tool
│
└── external/                   # Third-party source dependencies
    └── Pangolin-main/                        # Extracted Pangolin model repository
```

---

## 📦 3. Safely Archived Files (`./trash/src_archive/`)

In adherence to data preservation rules, the following 10 obsolete/redundant files were safely moved to `./trash/src_archive/`:

| Archived File | File Type & Reason | Replacement / Superseded By |
| :--- | :--- | :--- |
| `src/hpc/validate_output.sh` | 0-byte empty file | `tests/python/test_validate_output.py` |
| `src/external/pangolin.zip` | 167 MB redundant archive | Unzipped `src/external/Pangolin-main/` |
| `src/hpc/vcf2tsv.sh` | Deprecated legacy TSV parser | `src/hpc/vcf2parsed.sh` (Direct Parquet generation) |
| `src/hpc/tsv2xlsx.sh` | Deprecated legacy converter | Multi-format output in `universal_variant_converter.py` |
| `src/hpc/merge_vcfs.sh` | Redundant merge wrapper | `src/hpc/merge_vep_spip.sh` |
| `src/hpc/subset_gnomAD.sh` | One-off legacy utility | VEP direct gnomAD v4 tabix query |
| `src/hpc/UTRannotation_single_thread.sh` | Deprecated R script | VEP `--plugin UTRannotator` in `annotate_vep_vars.sh` |
| `src/R/refseq_to_ensembl.R` | Redundant translation script | `src/python/query_new_transcripts.py` |
| `src/python/refseq_to_ensembl.py` | Redundant translation script | `src/python/query_new_transcripts.py` |
| `src/python/filter_gene_transcript_by_missing.py` | Obsolete one-off helper | `src/python/query_new_transcripts.py` |

---

## 🧪 4. Testing & Verification Results

All tests across Python, R, and Bash were executed and verified:

```
============================================================================
 🧪 Master Pipeline Unit Testing Suite Runner 
============================================================================
 [1/3] Python Unit Tests (pytest)     : 27 / 27 PASSED (100%)
 [2/3] R Unit Tests (Rscript)         : 2 / 2 PASSED (testthat)
 [3/3] Bash Unit Tests (sh)           : 2 / 2 PASSED
============================================================================
 🎉 ALL UNIT TEST SUITES PASSED CLEANLY! (100% SUCCESS)
============================================================================
```
