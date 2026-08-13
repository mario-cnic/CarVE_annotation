# Direct VCF-to-Parquet Streamlining & Intermediate TSV Elimination

**Date**: 2026-08-13  
**Branch**: `feat/direct-vcf-to-parquet`  
**Scope**: Optimization of post-annotation processing stages (eliminating uncompressed multi-gigabyte intermediate TSV disk writes).

---

## 1. Problem & Architectural Rationale

Previously, the pipeline required **two separate SGE cluster jobs** to convert the annotated VCF into the final clean Parquet table:

```
PREVIOUS TWO-STAGE FLOW:
[Annotated VCF.gz] ──► (Job 1: vcf2tsv.sh) ──► [Massive Uncompressed TSV on Disk (2-10 GB)] ──► (Job 2: tsv2xlsx.sh) ──► [Final Clean Parquet (.pq)]
```

### Bottlenecks Identified:
1. **Network Disk I/O & Storage Overhead**: For whole-gene runs (e.g. `TTN`, `DSP`, `PKP2`), uncompressed TSVs consume **10 GB to 50+ GB** of shared NFS storage (`/data_lab_PGP/`).
2. **SGE Queue Latency**: Required two sequential 80 GB RAM job allocations per gene.
3. **Redundancy**: The intermediate TSV was never used for downstream biological analyses once `.pq` was written.

---

## 2. Implemented Architecture

```
OPTIMIZED DIRECT FLOW:
[Annotated VCF.gz] ──► (Single Job: vcf2parsed.sh) ──► [Final Clean Parquet (.pq) (~150-300 MB)]
```

1. **New Unified HPC Script [`src/hpc/vcf2parsed.sh`](../src/hpc/vcf2parsed.sh)**:
   - Uses high-speed local node scratch space (`/data_tmp/` or `/tmp/`) for ephemeral parsing.
   - Converts annotated VCF fields $\rightarrow$ immediately applies `filter_variants.py` (status contracts, splicing metrics, signed intron offsets) $\rightarrow$ writes directly to `results/${gene_name}.parsed.clean.${OUTPUT_FORMAT}`.
   - Automatically traps and cleans up any ephemeral scratch files upon completion or interrupt.

2. **Master Orchestrator [`main.sh`](../main.sh)**:
   - Replaced `vcf2tsv_job` and `tsv2xlsx_job` with a single unified `vcf2parsed_job`.
   - Downstream `filter_job`, `plot_job`, and `report_job` now hold directly on `vcf2parsed_job`.

3. **Plotting Compatibility [`src/R/plot_annotation_results.R`](../src/R/plot_annotation_results.R)**:
   - Added native support to read `.pq` / Parquet directly via `arrow::read_parquet()` or Python fallback.

---

## 3. Verification & Testing

- Master unit test suite (`bash tests/run_all_tests.sh`) executed on branch `feat/direct-vcf-to-parquet`:
  - **Python Unit Tests (`pytest`)**: 24/24 PASSED (100%).
  - **R Unit Tests (`Rscript`)**: PASSED (100%).
  - **Bash CLI Tests (`sh`)**: PASSED (100%).
