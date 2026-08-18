# Pipeline Diagnostic & Fix: gnomAD v4 Chromosome Alias Resolution & Clinical Report Integration

**Date**: 2026-08-18 14:36:00 CEST  
**Genome Assembly**: GRCh38  
**Components Modified**:
* [`src/hpc/annotate_vep_vars.sh`](../src/hpc/annotate_vep_vars.sh)
* [`src/hpc/subset_gnomAD.sh`](../src/hpc/subset_gnomAD.sh)
* [`main.sh`](../main.sh)
* [`src/python/generate_clinical_prioritization_report.py`](../src/python/generate_clinical_prioritization_report.py)

---

## 1. Problem Diagnosis

### Root Cause
1. **Contig Notation Mismatch**:
   * The whole-genome reference `/references/genomes/Homo_sapiens/annotations/gnomAD/GRCh38/v4.1/gnomAD.v4.1.vcf.gz` (805 GB) uses **UCSC chromosome notation** (`chr1`, `chr2`, ..., `chr10`, `chrX`).
   * Input cohort VCFs fed to VEP use **Ensembl chromosome notation** (`1`, `2`, ..., `10`, `X` without `chr`).
2. **Silent Failure**:
   * When VEP executed `--custom file=.../gnomAD.v4.1.vcf.gz,short_name=gnomADv4,format=vcf,type=exact,coords=0`, tabix queried `10:...` against an index containing only `chr10:...`.
   * Tabix returned 0 matches, leaving all `gnomADv4_*` fields empty (100% `NA`) across every gene without raising an error code.
3. **Downstream Impact**:
   * R burden testing and cohort enrichment scripts (`02_cohort_enrich.R`, `03_burden.R`, `07_domain_burden.R`) evaluated `is.na(gnomADv4_AF_joint)` as `TRUE` for all variants, bypassing rarity thresholds.

---

## 2. Solutions Implemented

### A. VEP `--synonyms` Flag in [`annotate_vep_vars.sh`](../src/hpc/annotate_vep_vars.sh)
Added the native VEP synonym mapping dictionary:
```bash
--offline --cache --cache_version 111 --dir ${vep_files}/cache \
--synonyms ${vep_files}/cache/homo_sapiens/111_GRCh38/chr_synonyms.txt \
--fasta ${gatk_bundle}/v0/Homo_sapiens_assembly38.fasta --assembly GRCh38 \
```
* Enables bidirectional alias translation (`10` $\leftrightarrow$ `chr10`, `X` $\leftrightarrow$ `chrX`) across all `--custom` VCF lookups.

### B. Optional Coordinate-Buffered Subsetting in [`subset_gnomAD.sh`](../src/hpc/subset_gnomAD.sh)
* Documented as an optional standalone utility for projects extracting whole-gene gnomAD catalogs.
* Added a $\pm 20\text{kb}$ buffer around gene boundaries to capture regulatory and deep intronic variants:
```bash
buf_start=$(( start > 20000 ? start - 20000 : 0 ))
buf_end=$(( end + 20000 ))
```
* Enforced `--rename-chrs /data_lab_PGP/shared/utils/data/map_chr.txt` to strip `chr` prefixes.

### C. Master Orchestrator Reporting Integration in [`main.sh`](../main.sh)
* Added automatic job submission for the 4-tab clinical prioritization dashboard in Step 5:
```bash
clinical_report_job=$(qsub -N "clinrep_${gene_name}" -P BIGN -A PGP -l h_vmem=20G -pe smp 1 \
    -hold_jid "$vcf2parsed_job" \
    -o "$ERROR_LOG_DIR/${gene_name}/${gene_name}.clinrep.out" \
    -e "$ERROR_LOG_DIR/${gene_name}/${gene_name}.clinrep.err" \
    -b y $PYTHON_EXE src/python/generate_clinical_prioritization_report.py --input "$final_output_file" --output "$REPORTS_MASTER_DIR/${gene_name}_clinical_prioritization_report.html" | awk '{print $3}')
```

---

## 3. Verification & Results

1. **Locus Variant Match Rate Verification**:
   * Evaluated `BAG3` (chr10) against `gnomAD.v4.1.vcf.gz`:
     * Total `BAG3` input variants: **20,650**
     * Exact matches in gnomAD v4.1 with chromosome alias alignment: **13,831 (67.0%)**
2. **Unit Test Suite**:
   * Executed `python3 -m pytest tests/python/`: **25/25 tests passed (100%)**.
