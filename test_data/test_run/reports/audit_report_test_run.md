# 🧪 Master Pipeline Execution & Quality Audit: `test_run`
**Audit Timestamp**: `2026-08-20 12:49:28`  
**Target Build**: `GRCh38`  
**Run Directory**: `test_data/test_run`  
**Results Directory**: `test_data/test_run/results`  

---

## 📊 1. Executive Summary
- **Total Target Genes**: `1`
- **Genes Fully Finished**: `1 / 1` (100.0%)
- **Genes In Progress / Pending**: `0`
- **Genes with Critical Errors**: `1`
- **Total Input Variants**: `50`
- **Total Clean Output Variants Generated**: `50`

---

## ⏱️ 2. Per-Gene Step Tracking & Running Times
| Gene   |   Variants | Current_Step   | SpliceAI_Progress   | Total_Elapsed   | Completed_At        |
|:-------|-----------:|:---------------|:--------------------|:----------------|:--------------------|
| MYBPC3 |         50 | Done           | ✅ Done (0.00 MB)   | ✅ 5m 01s       | 2026-08-20 12:18:41 |

---

## 📋 3. Stage-by-Stage Predictor Matrix
| Gene   | Status   | VarConv   | VEP     | SPiP    | Pangolin       | Branchpoint   | SpliceAI          | Final_PQ   | Log_Issues                              |
|:-------|:---------|:----------|:--------|:--------|:---------------|:--------------|:------------------|:-----------|:----------------------------------------|
| MYBPC3 | FINISHED | ✅        | ✅ Done | ✅ Done | ❌ Not created | ✅ Done       | ✅ Done (0.00 MB) | ✅ 50 vars | MYBPC3.pangolin.err, MYBPC3.clinrep.err |

---

## 🧬 4. Key Predictor Completeness (% Non-Null in Output Parquets)
| Gene   |   Columns | gnomADv4_AF_grpmax   | REVEL_score   | AlphaMissense   | SPiP   | Pangolin   | SpliceAI   | Branchpoint   | 5UTR_Annotator   | SpliceVault   | ClinVar   |
|:-------|----------:|:---------------------|:--------------|:----------------|:-------|:-----------|:-----------|:--------------|:-----------------|:--------------|:----------|
| MYBPC3 |       551 | 54.0%                | 0.0%          | 0.0%            | 100.0% | N/A        | 100.0%     | 100.0%        | 0.0%             | 100.0%        | 0.0%      |

---

## 📐 5. Schema & Column Integrity
- **Column Count Range**: `551` to `551` columns per parquet table.
- **Schema Consistency**: ✅ **100% Unified Schema** across all finished genes!

