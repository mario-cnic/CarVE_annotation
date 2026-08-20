# 📊 Comprehensive Multi-Pipeline Progress & Execution Audit
**Audit Execution Timestamp**: `2026-08-19 10:49:26 CEST`  
**Target Environment**: Local Sandbox & HPC SGE Cluster State  
**Genome Build**: `GRCh38` (Ensembl v111 / GATK Bundle v0)  

---
## 📌 Executive Cross-Pipeline Summary
| Pipeline Run | Total Genes | Completed | In Progress / Pending | Total Input Variants | Clean Output Records | Schema Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`run_20260813_1028`** | `7` | **`5`** (71.4%) | `2` | `395,858` | `1,939,179` | ⚠️ 512-551 cols |
| **`run_more_genes_20260817_1143`** | `118` | **`75`** (63.6%) | `43` | `3,939,837` | `972,678` | ⚠️ 534-551 cols |

## 🚀 Pipeline Run: `run_20260813_1028`
- **Path**: `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_20260813_1028`
- **Progress**: `5 / 7` genes complete (**71.4%**)
- **Total Input Variants Processed**: `395,858`
- **Total Output Consequence Records**: `1,939,179`
- **Stage Breakdown**:
  - `✅ Complete`: **5** genes
  - `3. Merge`: **2** genes

### Per-Gene Execution Status Matrix
| Gene   |   Input_Vars | Stage       | VarConv   | VEP   | SPiP   | Pangolin   | Branchpoint   | SpliceAI                           | Merged   | Final_PQ                   | Log_Issues                                              |
|:-------|-------------:|:------------|:----------|:------|:-------|:-----------|:--------------|:-----------------------------------|:---------|:---------------------------|:--------------------------------------------------------|
| BAG3   |        20650 | ✅ Complete | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.33 MB gz)               | ✅       | ✅ 51,872 rows (551 cols)  | None                                                    |
| DSP    |        29784 | 3. Merge    | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (2.00 MB uncompressed) | ❌       | ⏳ Pending                 | None                                                    |
| FLNC   |        25957 | ✅ Complete | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.40 MB gz)               | ✅       | ✅ 110,966 rows (551 cols) | None                                                    |
| LMNA   |        41098 | ✅ Complete | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.51 MB gz)               | ✅       | ✅ 982,298 rows (551 cols) | LMNA.spliceai.err, LMNA.spip.err                        |
| MYBPC3 |        28207 | 3. Merge    | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (1.50 MB uncompressed) | ❌       | ⏳ Pending                 | None                                                    |
| PKP2   |        63777 | ✅ Complete | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (1.13 MB gz)               | ✅       | ✅ 607,659 rows (551 cols) | PKP2.spliceai.err, PKP2.spip.err, PKP2.branchpoint.err  |
| TTN    |       186385 | ✅ Complete | ✅        | ✅    | ❌     | ❌         | ✅            | ✅ Done raw (2.54 MB uncompressed) | ✅       | ✅ 186,384 rows (512 cols) | TTN.branchpoint.err, TTN.pangolin.err, TTN.spliceai.err |

### Key Predictor Coverage in Completed Parquet Files
| Gene   |   Output_Records | gnomADv4   | REVEL   | AlphaMissense   | SPiP_cov   | Branch_cov   | SpliceAI_cov   | 5UTR_cov   | MaxEnt_cov   | Tier_cov   | Score_cov   |
|:-------|-----------------:|:-----------|:--------|:----------------|:-----------|:-------------|:---------------|:-----------|:-------------|:-----------|:------------|
| BAG3   |            51872 | 0.0%       | 3.0%    | 3.0%            | 89.3%      | 100.0%       | 100.0%         | 0.1%       | 0.3%         | 100.0%     | 100.0%      |
| FLNC   |           110966 | 0.0%       | 7.2%    | 7.4%            | 95.8%      | 100.0%       | 100.0%         | 0.0%       | 2.3%         | 100.0%     | 100.0%      |
| LMNA   |           982298 | 0.0%       | 0.7%    | 1.9%            | 94.8%      | 100.0%       | 100.0%         | 0.1%       | 0.9%         | 100.0%     | 100.0%      |
| PKP2   |           607659 | 0.0%       | 0.5%    | 1.0%            | 96.8%      | 100.0%       | 100.0%         | 0.0%       | 0.4%         | 100.0%     | 100.0%      |
| TTN    |           186384 | 0.0%       | 0.0%    | 28.4%           | ❌ Missing | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |

### Schema Consistency
- ⚠️ **Schema Variance**: Columns range from `512` to `551` across genes.
  - `TTN`: 512 cols

================================================================================

## 🚀 Pipeline Run: `run_more_genes_20260817_1143`
- **Path**: `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_more_genes_20260817_1143`
- **Progress**: `75 / 118` genes complete (**63.6%**)
- **Total Input Variants Processed**: `3,939,837`
- **Total Output Consequence Records**: `972,678`
- **Stage Breakdown**:
  - `✅ Complete`: **75** genes
  - `3. Merge`: **39** genes
  - `2. Predictors (SPiP)`: **2** genes
  - `2. Predictors (VEP)`: **1** genes
  - `2. Predictors (VEP/SPiP/Branch/SpliceAI)`: **1** genes

### Per-Gene Execution Status Matrix
| Gene     |   Input_Vars | Stage                                    | VarConv   | VEP   | SPiP   | Pangolin   | Branchpoint   | SpliceAI                            | Merged   | Final_PQ                  | Log_Issues                         |
|:---------|-------------:|:-----------------------------------------|:----------|:------|:-------|:-----------|:--------------|:------------------------------------|:---------|:--------------------------|:-----------------------------------|
| ACAD9    |        19913 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.35 MB gz)                | ✅       | ✅ 19,913 rows (534 cols) | None                               |
| ACADVL   |        15381 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.22 MB gz)                | ✅       | ✅ 15,167 rows (534 cols) | None                               |
| ACTA1    |         4465 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.05 MB gz)                | ✅       | ✅ 4,465 rows (534 cols)  | None                               |
| ACTA2    |        22705 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.25 MB gz)                | ✅       | ✅ 11,102 rows (534 cols) | None                               |
| ACTC1    |         8646 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.12 MB gz)                | ✅       | ✅ 8,646 rows (534 cols)  | None                               |
| AGK      |        32271 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (1.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| AGPAT2   |        14238 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.22 MB gz)                | ✅       | ✅ 14,238 rows (534 cols) | None                               |
| AKT1     |        26342 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.40 MB gz)                | ✅       | ✅ 26,342 rows (534 cols) | None                               |
| ALPK3    |        27020 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.45 MB gz)                | ✅       | ✅ 27,020 rows (534 cols) | None                               |
| ANKRD1   |         7694 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.11 MB gz)                | ✅       | ✅ 7,694 rows (534 cols)  | None                               |
| ATP5F1E  |         7928 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| ATPAF2   |        15753 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| BAG3     |        20650 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.33 MB gz)                | ✅       | ✅ 20,650 rows (534 cols) | None                               |
| BRAF     |        76742 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (5.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| BSCL2    |        16544 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.34 MB gz)                | ✅       | ✅ 16,544 rows (548 cols) | None                               |
| C10ORF71 |        14127 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.21 MB gz)                | ✅       | ✅ 14,127 rows (534 cols) | None                               |
| CACNA1C  |       228714 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (8.06 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| CALR     |         9688 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.13 MB gz)                | ✅       | ✅ 9,688 rows (534 cols)  | None                               |
| CAV3     |         7491 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.11 MB gz)                | ✅       | ✅ 7,491 rows (534 cols)  | None                               |
| CAVIN4   |         6664 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.10 MB gz)                | ✅       | ✅ 6,664 rows (534 cols)  | None                               |
| COA5     |         8373 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.13 MB gz)                | ✅       | ✅ 8,373 rows (551 cols)  | None                               |
| COA6     |         7540 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.11 MB gz)                | ✅       | ✅ 7,540 rows (534 cols)  | None                               |
| COQ2     |        11729 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.18 MB gz)                | ✅       | ✅ 11,729 rows (534 cols) | None                               |
| COX15    |        21332 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.31 MB gz)                | ✅       | ✅ 14,376 rows (534 cols) | None                               |
| COX6B1   |         6934 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.10 MB gz)                | ✅       | ✅ 6,934 rows (534 cols)  | None                               |
| CRYAB    |        12791 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| CSRP3    |        11844 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.16 MB gz)                | ✅       | ✅ 9,889 rows (534 cols)  | None                               |
| CTF1     |         9873 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.14 MB gz)                | ✅       | ✅ 9,873 rows (534 cols)  | None                               |
| CTNNA1   |        68175 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (5.28 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| CTNNB1   |        20669 | 2. Predictors (SPiP)                     | ✅        | ✅    | ❌     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | CTNNB1.spip.err                    |
| DES      |        10052 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.13 MB gz)                | ✅       | ✅ 10,052 rows (534 cols) | None                               |
| DLD      |        15997 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.25 MB gz)                | ✅       | ✅ 15,997 rows (534 cols) | None                               |
| DMD      |       533174 | 2. Predictors (VEP)                      | ✅        | ❌    | ✅     | ❌         | ✅            | ✅ Done raw (12.71 MB uncompressed) | ❌       | ⏳ Pending                | None                               |
| DNAJC19  |         7020 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| DOLK     |         5831 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.14 MB gz)                | ✅       | ✅ 5,831 rows (534 cols)  | None                               |
| DSC2     |        22753 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| DSG2     |        21674 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (1.00 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| EMD      |         4533 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.05 MB gz)                | ✅       | ✅ 4,533 rows (534 cols)  | None                               |
| FHL1     |        18565 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.19 MB gz)                | ✅       | ✅ 7,607 rows (534 cols)  | None                               |
| FHL2     |        31160 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| FHOD3    |       132380 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (4.13 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| FKRP     |        13909 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.22 MB gz)                | ✅       | ✅ 13,909 rows (534 cols) | None                               |
| FNIP1    |        50255 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (3.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| GATA5    |        16207 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.21 MB gz)                | ✅       | ✅ 16,207 rows (534 cols) | None                               |
| GATAD1   |        11198 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.16 MB gz)                | ✅       | ✅ 11,198 rows (534 cols) | None                               |
| GLA      |         8005 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Finished passes                  | ❌       | ⏳ Pending                | None                               |
| GYG1     |        18795 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.31 MB gz)                | ✅       | ✅ 18,795 rows (534 cols) | None                               |
| HAND1    |         4390 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.05 MB gz)                | ✅       | ✅ 4,390 rows (534 cols)  | None                               |
| HAND2    |         8929 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.11 MB gz)                | ✅       | ✅ 8,929 rows (551 cols)  | None                               |
| HCN4     |        23468 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.36 MB gz)                | ✅       | ✅ 23,468 rows (534 cols) | None                               |
| HFE      |         8591 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.13 MB gz)                | ✅       | ✅ 8,591 rows (534 cols)  | None                               |
| HRAS     |         8138 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.10 MB gz)                | ✅       | ✅ 8,138 rows (534 cols)  | None                               |
| IDH2     |        15385 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.23 MB gz)                | ✅       | ✅ 15,385 rows (534 cols) | None                               |
| ILK      |        12270 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.20 MB gz)                | ✅       | ✅ 12,269 rows (534 cols) | None                               |
| ISL1     |         9691 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.13 MB gz)                | ✅       | ✅ 9,691 rows (534 cols)  | None                               |
| JPH2     |        34343 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| JUP      |        21906 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.31 MB gz)                | ✅       | ✅ 21,906 rows (534 cols) | None                               |
| KCNJ2    |         7951 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.09 MB gz)                | ✅       | ✅ 7,951 rows (534 cols)  | None                               |
| KCNJ8    |         5851 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.08 MB gz)                | ✅       | ✅ 5,851 rows (534 cols)  | None                               |
| KLF10    |         8393 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.11 MB gz)                | ✅       | ✅ 8,393 rows (534 cols)  | None                               |
| KLHL24   |        21337 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Finished passes                  | ❌       | ⏳ Pending                | None                               |
| KRAS     |        21744 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.31 MB gz)                | ✅       | ✅ 21,743 rows (534 cols) | None                               |
| LAMA2    |       187093 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (5.00 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| LAMP2    |        14811 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.27 MB gz)                | ✅       | ✅ 14,811 rows (534 cols) | None                               |
| LEMD2    |        15856 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| LIAS     |        21521 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.30 MB gz)                | ✅       | ✅ 21,521 rows (534 cols) | None                               |
| LMOD2    |         6176 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.09 MB gz)                | ✅       | ✅ 6,176 rows (534 cols)  | None                               |
| LRRC10   |         4031 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.05 MB gz)                | ✅       | ✅ 4,031 rows (534 cols)  | None                               |
| MEF2C    |        62271 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (5.00 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| MLYCD    |        29980 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (4.00 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| MRAS     |        20954 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.33 MB gz)                | ✅       | ✅ 20,954 rows (534 cols) | None                               |
| MRPL3    |        16372 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.26 MB gz)                | ✅       | ✅ 16,372 rows (534 cols) | None                               |
| MRPL44   |         8592 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.11 MB gz)                | ✅       | ✅ 7,359 rows (534 cols)  | None                               |
| MRPS22   |        11364 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.18 MB gz)                | ✅       | ✅ 11,364 rows (534 cols) | None                               |
| MTO1     |        28908 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.46 MB gz)                | ✅       | ✅ 28,907 rows (548 cols) | None                               |
| MYBPHL   |         8657 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.13 MB gz)                | ✅       | ✅ 8,657 rows (534 cols)  | None                               |
| MYH6     |        24228 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.36 MB gz)                | ✅       | ✅ 24,228 rows (534 cols) | None                               |
| MYH7     |        18996 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Finished passes                  | ❌       | ⏳ Pending                | None                               |
| MYL2     |         7108 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| MYL3     |        12149 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.15 MB gz)                | ✅       | ✅ 5,796 rows (534 cols)  | None                               |
| MYLK2    |        11422 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.18 MB gz)                | ✅       | ✅ 11,422 rows (534 cols) | None                               |
| MYLK3    |        27273 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.42 MB gz)                | ✅       | ✅ 23,112 rows (534 cols) | None                               |
| MYOM1    |        65735 | 2. Predictors (SPiP)                     | ✅        | ✅    | ❌     | ❌         | ✅            | ✅ Done raw (5.50 MB uncompressed)  | ❌       | ⏳ Pending                | MYOM1.spip.err, MYOM1.pangolin.err |
| MYOT     |        12580 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.19 MB gz)                | ✅       | ✅ 12,580 rows (534 cols) | None                               |
| MYOZ2    |        15989 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.25 MB gz)                | ✅       | ✅ 15,989 rows (534 cols) | None                               |
| MYPN     |        44062 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| MYZAP    |        34561 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (3.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| NAA10    |         9920 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.14 MB gz)                | ✅       | ✅ 9,920 rows (534 cols)  | None                               |
| NDUFB11  |         5046 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.08 MB gz)                | ✅       | ✅ 5,046 rows (534 cols)  | None                               |
| NEBL     |       163579 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (5.00 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| NEXN     |        25776 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.38 MB gz)                | ✅       | ✅ 25,776 rows (534 cols) | None                               |
| NF1      |       130452 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (5.00 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| NKX2-5   |         6182 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.08 MB gz)                | ✅       | ✅ 6,182 rows (534 cols)  | None                               |
| NONO     |        11156 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.18 MB gz)                | ✅       | ✅ 11,156 rows (534 cols) | None                               |
| NRAS     |         9100 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.14 MB gz)                | ✅       | ✅ 9,100 rows (534 cols)  | None                               |
| NUBPL    |        89653 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| PCCA     |       130352 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (5.00 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| PCCB     |        35701 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (1.00 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| PDHA1    |        12909 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.22 MB gz)                | ✅       | ✅ 12,909 rows (534 cols) | None                               |
| PDLIM3   |        17415 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.30 MB gz)                | ✅       | ✅ 17,415 rows (534 cols) | None                               |
| PERP     |         8728 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| PHKB     |        74315 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (4.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| PLEKHM2  |        32142 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (2.00 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| PLN      |         6486 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.16 MB gz)                | ✅       | ✅ 6,486 rows (534 cols)  | None                               |
| PMM2     |        32589 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (1.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| PPA2     |        35080 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (1.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| PPCS     |         9459 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.10 MB gz)                | ✅       | ✅ 6,779 rows (534 cols)  | None                               |
| PPP1CB   |        21487 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.33 MB gz)                | ✅       | ✅ 21,487 rows (534 cols) | None                               |
| PPP1R13L |        22020 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.38 MB gz)                | ✅       | ✅ 22,020 rows (534 cols) | None                               |
| PRDM16   |       140750 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (7.18 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| PRKAG2   |       109263 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (5.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| PSEN2    |        16797 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.27 MB gz)                | ✅       | ✅ 16,797 rows (534 cols) | None                               |
| PTPN11   |        33214 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (2.00 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| QRSL1    |        18958 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.31 MB gz)                | ✅       | ✅ 18,958 rows (534 cols) | None                               |
| RAF1     |        37974 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (0.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| RBM20    |        60695 | 3. Merge                                 | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done raw (5.50 MB uncompressed)  | ❌       | ⏳ Pending                | None                               |
| RIT1     |        10069 | ✅ Complete                              | ✅        | ✅    | ✅     | ❌         | ✅            | ✅ Done (0.15 MB gz)                | ✅       | ✅ 10,069 rows (534 cols) | None                               |
| RPL3L    |        15750 | 2. Predictors (VEP/SPiP/Branch/SpliceAI) | ✅        | ❌    | ❌     | ❌         | ❌            | ⏳ Pending                          | ❌       | ⏳ Pending                | None                               |

### Key Predictor Coverage in Completed Parquet Files
| Gene     |   Output_Records | gnomADv4   | REVEL   | AlphaMissense   | SPiP_cov   | Branch_cov   | SpliceAI_cov   | 5UTR_cov   | MaxEnt_cov   | Tier_cov   | Score_cov   |
|:---------|-----------------:|:-----------|:--------|:----------------|:-----------|:-------------|:---------------|:-----------|:-------------|:-----------|:------------|
| ACAD9    |            19913 | 0.0%       | 4.7%    | 4.8%            | 95.5%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| ACADVL   |            15167 | 0.0%       | 7.3%    | 7.6%            | 99.1%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| ACTA1    |             4465 | 0.0%       | 6.2%    | 6.5%            | 57.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| ACTA2    |            11102 | 0.0%       | 2.7%    | 2.8%            | 96.8%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| ACTC1    |             8646 | 0.0%       | 2.9%    | 2.9%            | 99.9%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| AGPAT2   |            14238 | 0.0%       | 4.8%    | 5.1%            | 88.9%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| AKT1     |            26342 | 0.0%       | 0.0%    | 5.6%            | 92.6%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| ALPK3    |            27020 | 0.0%       | 11.6%   | 11.6%           | 93.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| ANKRD1   |             7694 | 0.0%       | 6.1%    | 6.1%            | 79.5%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| BAG3     |            20650 | 0.0%       | 4.9%    | 5.0%            | 89.6%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| BSCL2    |            16544 | 0.0%       | 4.1%    | 4.3%            | 97.7%      | 100.0%       | 100.0%         | 0.2%       | ❌ Missing   | 100.0%     | 100.0%      |
| C10ORF71 |            14127 | 0.0%       | 16.4%   | 16.4%           | 93.2%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| CALR     |             9688 | 0.0%       | 6.4%    | 6.4%            | 80.8%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| CAV3     |             7491 | 0.0%       | 3.0%    | 3.2%            | 76.1%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| CAVIN4   |             6664 | 0.0%       | 8.6%    | 8.6%            | 84.5%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| COA5     |             8373 | 0.0%       | 1.6%    | 1.8%            | 88.5%      | 100.0%       | 100.0%         | 0.1%       | 0.6%         | 100.0%     | 100.0%      |
| COA6     |             7540 | 0.0%       | 4.7%    | 4.9%            | 79.9%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| COQ2     |            11729 | 0.0%       | 0.0%    | 8.0%            | 84.7%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| COX15    |            14376 | 0.0%       | 3.6%    | 3.8%            | 100.0%     | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| COX6B1   |             6934 | 0.0%       | 0.0%    | 1.5%            | 69.5%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| CSRP3    |             9889 | 0.0%       | 2.9%    | 3.1%            | 92.3%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| CTF1     |             9873 | 0.0%       | 8.2%    | 8.6%            | 78.9%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| DES      |            10052 | 0.0%       | 10.0%   | 10.0%           | 68.9%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| DLD      |            15997 | 0.0%       | 4.3%    | 4.5%            | 86.8%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| DOLK     |             5831 | 0.0%       | 12.6%   | 12.6%           | 84.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| EMD      |             4533 | 0.0%       | 8.7%    | 8.9%            | 53.0%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| FHL1     |             7607 | 0.0%       | 3.8%    | 3.9%            | 93.5%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| FKRP     |            13909 | 0.0%       | 10.9%   | 11.4%           | 92.9%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| GATA5    |            16207 | 0.0%       | 8.9%    | 8.9%            | 67.3%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| GATAD1   |            11198 | 0.0%       | 5.0%    | 5.0%            | 82.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| GYG1     |            18795 | 0.0%       | 3.0%    | 3.2%            | 89.2%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| HAND1    |             4390 | 0.0%       | 9.8%    | 9.8%            | 61.1%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| HAND2    |             8929 | 0.0%       | 6.1%    | 6.4%            | 80.6%      | 100.0%       | 100.0%         | 1.6%       | 0.3%         | 100.0%     | 100.0%      |
| HCN4     |            23468 | 0.0%       | 12.2%   | 12.2%           | 93.2%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| HFE      |             8591 | 0.0%       | 5.3%    | 5.7%            | 57.2%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| HRAS     |             8138 | 0.0%       | 3.3%    | 3.5%            | 62.0%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| IDH2     |            15385 | 0.0%       | 5.2%    | 5.4%            | 89.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| ILK      |            12269 | 0.0%       | 5.2%    | 5.5%            | 92.2%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| ISL1     |             9691 | 0.0%       | 4.5%    | 4.6%            | 91.8%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| JUP      |            21906 | 0.0%       | 5.4%    | 5.5%            | 90.1%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| KCNJ2    |             7951 | 0.0%       | 5.4%    | 5.7%            | 87.8%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| KCNJ8    |             5851 | 0.0%       | 6.3%    | 2.1%            | 71.2%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| KLF10    |             8393 | 0.0%       | 8.1%    | 8.4%            | 73.3%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| KRAS     |            21743 | 0.0%       | 0.7%    | 0.8%            | 95.8%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| LAMP2    |            14811 | 0.0%       | 3.3%    | 3.4%            | 90.6%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| LIAS     |            21521 | 0.0%       | 0.0%    | 2.9%            | 88.1%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| LMOD2    |             6176 | 0.0%       | 15.6%   | 16.1%           | 77.1%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| LRRC10   |             4031 | 0.0%       | 9.4%    | 9.4%            | 59.1%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MRAS     |            20954 | 0.0%       | 1.2%    | 1.3%            | 92.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MRPL3    |            16372 | 0.0%       | 3.4%    | 3.6%            | 91.0%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MRPL44   |             7359 | 0.0%       | 7.1%    | 7.1%            | 65.5%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MRPS22   |            11364 | 0.0%       | 0.0%    | 5.0%            | 87.6%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MTO1     |            28907 | 0.0%       | 3.4%    | 3.6%            | 67.8%      | 100.0%       | 100.0%         | 0.1%       | ❌ Missing   | 100.0%     | 100.0%      |
| MYBPHL   |             8657 | 0.0%       | 6.4%    | 6.4%            | 81.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MYH6     |            24228 | 0.0%       | 11.2%   | 11.7%           | 91.6%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MYL3     |             5796 | 0.0%       | 4.6%    | 4.7%            | 56.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MYLK2    |            11422 | 0.0%       | 7.8%    | 8.1%            | 85.0%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MYLK3    |            23112 | 0.0%       | 6.0%    | 6.1%            | 96.9%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MYOT     |            12580 | 0.0%       | 4.9%    | 5.2%            | 85.5%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| MYOZ2    |            15989 | 0.0%       | 2.0%    | 2.0%            | 90.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| NAA10    |             9920 | 0.0%       | 2.0%    | 2.1%            | 74.7%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| NDUFB11  |             5046 | 0.0%       | 4.3%    | 4.5%            | 85.8%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| NEXN     |            25776 | 0.0%       | 3.9%    | 4.2%            | 85.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| NKX2-5   |             6182 | 0.0%       | 10.9%   | 10.9%           | 70.0%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| NONO     |            11156 | 0.0%       | 3.2%    | 3.4%            | 92.1%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| NRAS     |             9100 | 0.0%       | 2.2%    | 2.2%            | 92.5%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| PDHA1    |            12909 | 0.0%       | 2.9%    | 3.0%            | 93.9%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| PDLIM3   |            17415 | 0.0%       | 1.0%    | 3.6%            | 87.2%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| PLN      |             6486 | 0.0%       | 1.2%    | 1.2%            | 100.0%     | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| PPCS     |             6779 | 0.0%       | 7.7%    | 8.2%            | 100.0%     | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| PPP1CB   |            21487 | 0.0%       | 0.8%    | 0.9%            | 88.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| PPP1R13L |            22020 | 0.0%       | 10.9%   | 11.3%           | 94.8%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| PSEN2    |            16797 | 0.0%       | 3.9%    | 4.0%            | 84.0%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| QRSL1    |            18958 | 0.0%       | 3.4%    | 3.5%            | 95.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |
| RIT1     |            10069 | 0.0%       | 2.7%    | 2.8%            | 82.4%      | 100.0%       | 100.0%         | ❌ Missing | ❌ Missing   | 100.0%     | 100.0%      |

### Schema Consistency
- ⚠️ **Schema Variance**: Columns range from `534` to `551` across genes.
  - `ACAD9`: 534 cols
  - `ACADVL`: 534 cols
  - `ACTA1`: 534 cols
  - `ACTA2`: 534 cols
  - `ACTC1`: 534 cols
  - `AGPAT2`: 534 cols
  - `AKT1`: 534 cols
  - `ALPK3`: 534 cols
  - `ANKRD1`: 534 cols
  - `BAG3`: 534 cols
  - `BSCL2`: 548 cols
  - `C10ORF71`: 534 cols
  - `CALR`: 534 cols
  - `CAV3`: 534 cols
  - `CAVIN4`: 534 cols
  - `COA6`: 534 cols
  - `COQ2`: 534 cols
  - `COX15`: 534 cols
  - `COX6B1`: 534 cols
  - `CSRP3`: 534 cols
  - `CTF1`: 534 cols
  - `DES`: 534 cols
  - `DLD`: 534 cols
  - `DOLK`: 534 cols
  - `EMD`: 534 cols
  - `FHL1`: 534 cols
  - `FKRP`: 534 cols
  - `GATA5`: 534 cols
  - `GATAD1`: 534 cols
  - `GYG1`: 534 cols
  - `HAND1`: 534 cols
  - `HCN4`: 534 cols
  - `HFE`: 534 cols
  - `HRAS`: 534 cols
  - `IDH2`: 534 cols
  - `ILK`: 534 cols
  - `ISL1`: 534 cols
  - `JUP`: 534 cols
  - `KCNJ2`: 534 cols
  - `KCNJ8`: 534 cols
  - `KLF10`: 534 cols
  - `KRAS`: 534 cols
  - `LAMP2`: 534 cols
  - `LIAS`: 534 cols
  - `LMOD2`: 534 cols
  - `LRRC10`: 534 cols
  - `MRAS`: 534 cols
  - `MRPL3`: 534 cols
  - `MRPL44`: 534 cols
  - `MRPS22`: 534 cols
  - `MTO1`: 548 cols
  - `MYBPHL`: 534 cols
  - `MYH6`: 534 cols
  - `MYL3`: 534 cols
  - `MYLK2`: 534 cols
  - `MYLK3`: 534 cols
  - `MYOT`: 534 cols
  - `MYOZ2`: 534 cols
  - `NAA10`: 534 cols
  - `NDUFB11`: 534 cols
  - `NEXN`: 534 cols
  - `NKX2-5`: 534 cols
  - `NONO`: 534 cols
  - `NRAS`: 534 cols
  - `PDHA1`: 534 cols
  - `PDLIM3`: 534 cols
  - `PLN`: 534 cols
  - `PPCS`: 534 cols
  - `PPP1CB`: 534 cols
  - `PPP1R13L`: 534 cols
  - `PSEN2`: 534 cols
  - `QRSL1`: 534 cols
  - `RIT1`: 534 cols

================================================================================
