# 🧬 Compound Heterozygosis Scientific Review & Refactoring Walkthrough
**Date**: `2026-08-31 12:25:00 CEST`  
**Genome Assembly**: `GRCh38` (Ensembl v111 / GATK Bundle v0)  
**Location**: `/data_lab_PGP/pipelines/annotation_pipeline_new/` (HPC) $\leftrightarrow$ `/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/` (Local Sandbox)  

---

## 📌 Executive Summary

We conducted a deep scientific literature review of **compound heterozygosity** according to ACMG/ClinGen standards and audited the implementation in [`filter_variants.py`](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py#L1883).

We corrected a critical logical flaw where **CIS** variants (inherited from the same parent) were previously being incorrectly flagged as compound heterozygous. We replaced this with a strict **pairwise TRANS matching engine** (`_mark_pairwise_trans()`), created a new unit test suite ([`tests/python/test_compound_het.py`](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/tests/python/test_compound_het.py)), and verified end-to-end multi-sample trio prioritization.

---

## 🔍 Key Literature Standards & Code Audit

1. **Strict TRANS Phasing Requirement**:
   - Compound heterozygosity occurs when an individual inherits two distinct mutant alleles within the same gene: one from the father (**paternal HET**, maternal HOMREF) and one from the mother (**maternal HET**, paternal HOMREF).
   - If both variants are inherited from the *same* parent (**CIS phasing**), the individual retains one intact, wild-type allele from the other parent and is **NOT compound heterozygous** for recessive disease.

2. **Fixed Code Logic Flaw**:
   - Previously, `_mark_with_dict()` used `or` logic across group variants, which marked CIS variants if a father had 2 HET variants in the same gene.
   - **Fix**: Implemented `_mark_pairwise_trans()`, which checks that a gene has at least **one paternal-only variant** AND at least **one maternal-only variant** for a `HET` proband.

3. **Unphased Single-Sample / Missing Parent Handling**:
   - For single-sample datasets or missing parents, genes with $\ge 2$ HET variants are flagged as `compound_het_candidate = True` (Unphased Candidate) requiring molecular/read-backed phasing.

4. **Support for Multiple Affected Family Members**:
   - `_mark_with_pandas()` now iterates across all affected individuals in a pedigree dataframe without throwing `ValueError` exceptions.

---

## 🧪 Unit Test & Verification Results

Executed `tests/python/test_compound_het.py` (3 tests passed in 0.029s):
- ✅ **Test 1 (TRANS Phasing)**: Variant A (Paternal) + Variant B (Maternal) $\rightarrow$ `compound_het = True` (Inheritance: `Compound Heterozygous`, Priority Score: 50.0).
- ✅ **Test 2 (CIS Phasing)**: Variant A (Paternal) + Variant B (Paternal) $\rightarrow$ `compound_het = False` (Inheritance: `Unclassified`).
- ✅ **Test 3 (Single Sample Unphased)**: 2 HET variants in same gene $\rightarrow$ `compound_het_candidate = True`.
