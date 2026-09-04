# 🧬 Pedigree & Multi-Sample Inheritance Analysis Engine Walkthrough
**Date**: `2026-08-31 11:45:00 CEST`  
**Genome Assembly**: `GRCh38` (Ensembl v111 / GATK Bundle v0)  
**Location**: `/data_lab_PGP/pipelines/annotation_pipeline_new/` (HPC) $\leftrightarrow$ `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/` (Local Sandbox)  

---

## 📌 Executive Summary

We developed and integrated a comprehensive **Pedigree & Multi-Sample Inheritance Engine** ([`analyze_pedigree_inheritance()`](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py#L2054)) into the annotation pipeline and Standalone Web Application ([`src/web_app/app.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/web_app/app.py)).

Clinicians and bioinformaticians can now evaluate genetic inheritance models across multi-sample VCFs and variant tables, searching for **De Novo** variants in trios, segregating **Autosomal Dominant** or **Autosomal Recessive** variants, **Compound Heterozygous** variants, and **X-Linked** mutations.

If single-sample data is uploaded or no pedigree information is provided, the engine **gracefully falls back** to single-sample clinical prioritization without throwing errors or requiring configuration.

---

## 🛠️ Supported Inheritance Models

| Inheritance Model | Genotype Requirements | Priority Score Bonus |
| :--- | :--- | :---: |
| **`De Novo`** | Proband is `HET` or `HOMALT`; both Father & Mother are `HOMREF` (`0/0`). | **+20 pts** |
| **`Autosomal Recessive (Hom)`** | Proband is `HOMALT`; both Father & Mother are carrier `HET`. | **+15 pts** |
| **`Compound Heterozygous`** | Proband has $\ge 2$ complementary heterozygous mutations in the same gene. | **+15 pts** |
| **`Autosomal Dominant`** | Present (`HET` or `HOMALT`) in all affected members; absent (`HOMREF`) in unaffected members. | **+15 pts** |
| **`X-Linked`** | Variant on ChrX; male proband is hemizygous `HET`/`HOMALT`, mother is carrier `HET`. | **+15 pts** |
| **`Shared Affected`** | Present in all specified affected individuals. | **+10 pts** |
| **`Single Sample / Not Applicable`** | Default graceful fallback for single-sample VCFs or tables. | — |

---

## 💻 Web App & Dashboard Integration

1. **Auto-Detection**: The Web App automatically checks if uploaded files contain multiple sample columns (`GT_<SampleID>`).
2. **Pedigree Upload & Role Mapping**:
   - Upload `.ped` files (standard 6-column PED format) or `.json` files.
   - Interactive **Sample Mapping Interface** to select `Proband (Index)`, `Father`, `Mother`, `Affected Relatives`, and `Unaffected Relatives`.
3. **Multi-Column Filtering**: Filter data tables by `INHERITANCE_MODEL` in addition to Gene, Priority Tier, and Splicing scores.
4. **HTML Dashboard & Export**: Embedded HTML reports and exported `.xlsx` / `.pq` tables include sample genotype summaries and inheritance badges.
