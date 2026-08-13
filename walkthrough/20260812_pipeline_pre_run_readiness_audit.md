# Pipeline Pre-Run Readiness & Verification Audit Report

**Date**: 2026-08-12 15:12:35  
**Build**: GRCh38 / Ensembl v104 / RefSeq MANE v1.3  
**Status**: ✅ ALL CHECKS PASSED (100% READY)

---

## 1. System & Architecture Inventory

| Component | Status | Size (MB) | Path |
| :--- | :--- | :--- | :--- |
| Pangolin GRCh38 Database | ✅ Ready | 917.36 | `/home/mruizp/data_lab_PGP/shared/utils/pangolin_db/pangolin_grch38.db` |
| LaBranchoR GRCh38 Top BP Database | ✅ Ready | 2.49 | `/home/mruizp/data_lab_PGP/resources/annotation/labranchor/labranchor_grch38_top.bed.gz` |
| LaBranchoR Top BP Index | ✅ Ready | 0.43 | `/home/mruizp/data_lab_PGP/resources/annotation/labranchor/labranchor_grch38_top.bed.gz.tbi` |
| LaBranchoR GRCh38 ISM Database | ✅ Ready | 366.33 | `/home/mruizp/data_lab_PGP/resources/annotation/labranchor/labranchor_grch38_ism.tsv.gz` |
| SpliceVault GRCh38 Database | ✅ Ready | 842.96 | `/home/mruizp/data_lab_PGP/resources/annotation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz` |
| SpliceVarDB Reference VCF | ✅ Ready | 0.69 | `/home/mruizp/data_lab_PGP/resources/annotation/splicevardb/splicevardb.to.annotate.tsv.vcf.gz` |
| GRCh38 Reference FASTA | ✅ Ready | 3099.36 | `/home/mruizp/data_references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta` |
| Master Orchestration Script | ✅ Ready | 0.01 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/main.sh` |
| Universal Variant Converter Wrapper | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/variant_converter.sh` |
| VEP Annotation Wrapper | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/annotate_vep_vars.sh` |
| SPiP v2.1 Annotation Wrapper | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/annotate_spip_vars.sh` |
| Pangolin Annotation Wrapper | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/annotate_pangolin_vars.sh` |
| SpliceAI (-D 10000) Wrapper | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/annotate_spliceai_vars.sh` |
| Branchpoint Predictor Wrapper | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/annotate_branchpointer_vars.sh` |
| Multi-Predictor Annotation Merger | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/merge_vep_spip.sh` |
| VCF to TSV Parsing Wrapper | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/vcf2tsv.sh` |
| TSV to Parquet/XLSX Clean Filter Wrapper | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/hpc/tsv2xlsx.sh` |
| Pangolin Conda Spec | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/conda/pangolin_env.yml` |
| SpliceAI Conda Spec | ✅ Ready | 0.0 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/conda/spliceai_env.yml` |
| VCF Parser Column Inclusions | ✅ Ready | 0.01 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/all_but_old_gnomad_vep_cols.txt` |
| Column Dictionary & Data Schema | ✅ Ready | 0.01 | `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/column_dictionary.md` |
| Filter Variants Logic & Schema | ✅ Ready | 0.08 | `/home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py` |

---

## 2. Integrated Splicing & Functional Predictor Matrix

| Predictor | Version / Spec | Target Genome Window / Metric | Output Columns | Status Contract |
| :--- | :--- | :--- | :--- | :--- |
| **SPiP** | v2.1 (Multi-threaded) | Donor / Acceptor / Exonic | `SPiP_interpretation`, `SPiP_prediction`, `SPiP_score`, `SPiP_mechanism` | `scored`, `not_covered`, `error` |
| **SpliceVault** | Empirical GTEx/SRA | RNA Aberrant Events & Cryptic Splice | `SpliceVault_top_events`, `SpliceVault_Predictions_Decoded`, `SpliceVault_site_sample_count`, `SpliceVault_SpliceAI_delta` | `aberrant_event_detected`, `no_events_found`, `not_covered` |
| **Intron Offset** | Refactored Signed | Distance to Nearest Splice Junction | `intron_offset_signed`, `splice_side` | Exact integer (neg = intron, pos = exon) |
| **Pangolin** | PyTorch Ensemble | 20kb Window (`-d 10000`) | `Pangolin_max_score`, `Pangolin_heart_lv_score`, `Pangolin_heart_aa_score` | `scored`, `not_covered`, `error` |
| **SpliceAI Local** | TensorFlow 2.15 | 20kb Window (`-D 10000`) | `spliceAI_MAX`, `SpliceAI_status` | `scored`, `not_covered`, `error` |
| **Branchpointer** | Signal-Feature Model | -18 to -44 bp Upstream of 3' SS | `Branchpointer_prob`, `Branchpointer_U2_energy`, `Branchpoint_disrupted` | `scored`, `not_covered`, `error` |
| **LaBranchoR** | Bi-LSTM Deep Learning | 206,249 Top BP Sites (GRCh38) | `LaBranchoR_score`, `LaBranchoR_acc_dist` | `scored`, `not_covered`, `error` |

---

## 3. Unit Test Verification

- **Python Suite (`pytest`)**: 24/24 unit tests passed cleanly (100%).
- **R Suite (`Rscript`)**: Publication plot generation and statistical tests passed cleanly.
- **Bash Suite (`sh`)**: Parameter validation and logging tests passed cleanly.
- **Result**: **100% SUCCESS across all test suites**.

---

## 4. Pipeline Execution Guidance

To launch the full pipeline run across your target raw directory on the SGE cluster:

```bash
bash main.sh --raw-dir /path/to/raw_variants/ --build GRCh38 --output-format pq
```

All jobs (`varconv`, `vep`, `spip`, `pangolin`, `spliceai`, `branchpoint`, `merge`, `vcf2tsv`, `tsv2xlsx`, `filter`, `plot`, `report`) will automatically execute in parallel DAG order.
