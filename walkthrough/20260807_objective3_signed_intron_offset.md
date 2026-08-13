# Objective 3 Completion Walkthrough: Signed `intron_offset` & Splice Side Refactor

**Date**: 2026-08-07  
**Author**: Bioinformatics Pipeline Team  
**Status**: ✅ COMPLETED  

---

## 1. Executive & Biological Rationale

Objective 3 is a **prerequisite structural refactor** of how intronic variant locations are parsed and represented throughout the pipeline.

### The Problem with Absolute Offsets
Previously, `intron_offset` was recorded as an **absolute positive integer**:
- `c.123+5C>T` (located 5 bp downstream of an exon, on the **Donor** side) $\rightarrow$ `intron_offset = 5`
- `c.123-5C>T` (located 5 bp upstream of an exon, on the **Acceptor** side) $\rightarrow$ `intron_offset = 5`

Because both variants produced `intron_offset = 5`, downstream pipeline steps could not distinguish donor-side intronic variants from acceptor-side intronic variants based on `intron_offset` alone.

### Biological Asymmetry Across Introns
Intronic splicing elements are asymmetric across human introns:
```
Exon 1 ===> [ Donor Motif (+1 to +6) ] ----------- Intron ----------- [ Branch Point (-18 to -44) | PPT (-10 to -25) | Acceptor Motif (-1 to -3) ] ===> Exon 2
```
- **Acceptor-Side-Only Regulatory Elements**: The **Branch Point (BP)** (located $-18$ to $-44\text{ nt}$ from exon) and the **Polypyrimidine Tract (PPT)** (located $-10$ to $-25\text{ nt}$ from exon) exist **ONLY on the acceptor side** (upstream of the exon).
- **Donor Side**: Contains the 5' donor splice motif ($+1$ to $+6\text{ nt}$).

Without a signed offset, acceptor-side-specific tools (like **Branchpointer** or **LaBranchoR**) would evaluate donor-side variants (`+18` to `+44`), producing invalid predictions and false positives, violating the `not_applicable` status contract.

---

## 2. Summary of Architectural Changes Made

1. **cDNA / HGVS Variant Converter** ([universal_variant_converter.py](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/src/python/universal_variant_converter.py)):
   - Implemented `compute_signed_intron_offset()` using regex parsing on HGVS string notation (`c.\d+([+-])(\d+)`).
   - Exports `INTRON_OFFSET_SIGNED` and `SPLICE_SIDE` into intermediate VCF INFO fields.

2. **Downstream Pipeline Post-Processing** ([filter_variants.py](file:///home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py)):
   - Implemented `build_signed_intron_offset(data)` and registered new standard columns:
     - `intron_offset_signed`: Positive for donor side (`+5`), negative for acceptor side (`-12`), `0` for exonic.
     - `splice_side`: Categorical Enum (`donor`, `acceptor`, `exonic`).
   - Positioned `intron_offset_signed` and `splice_side` immediately following `CDNA_NAME` in output Parquet schema (`first_cols`).

---

## 3. Verification & Unit Test Results

Ran unit tests on cDNA test cases:

| Input cDNA HGVS String | Derived `intron_offset_signed` | Derived `splice_side` | Status |
| --- | --- | --- | --- |
| `NM_000257.3:c.526+5C>T` | `+5` | `donor` | ✅ PASS |
| `ENST00000343260:c.100-12A>G` | `-12` | `acceptor` | ✅ PASS |
| `c.1504C>T` | `0` | `exonic` | ✅ PASS |
| `c.456-1G>A` | `-1` | `acceptor` | ✅ PASS |
| `c.789+2T>C` | `+2` | `donor` | ✅ PASS |

All unit tests and DataFrame transformation checks passed cleanly.
