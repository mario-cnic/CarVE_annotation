# Genotype tracking in the Nextflow output (2026-10-02)

Item: `TODO.md` "Genotype tracking" (Mario, priority). Plan approved 2026-10-02.

## Problem
`VCF_TO_TABLE` (`modules/local/vcf_to_table.nf`) calls the parser with `--add_info --add_vep` only, so no table ever had genotypes. The main table is variant x VEP transcript (panel7: 901 rows for 72 variants, up to 29 per variant), Module 3 reads wide `GT_<sample>` columns (`HET/HOMALT/HOMREF/MISSING`), and `carve-platform` plans a separate long genotype layer (`PLT-022`, `PLT-027`).

## Checked before designing
- Sample columns survive every stage that feeds the final VCF (chunker, VEP, `concat -a` + `sort`, `annotate -a`): all six VCFs of the real `S223` run carry `S223_S223` and `GT:AD:DP:GQ:PL`; a local run of the chunker, concat/sort and annotate on the 3-sample trio fixture left every `GT:GQ` unchanged. Not verified: VEP and the SPiP R script with more than one sample.
- The parser's existing wide mode (`--add_gt` etc.) also drops rows where no selected sample has a "1" allele, and `parse_genotype` leaves `1/2`, `./1` and haploid `1` unmapped.

## Decisions (Mario)
1. Separate long genotype table is canonical; an optional wide `GT_<sample>` view keeps Module 3 working for family-sized inputs.
2. Every variant x sample pair is written, with an explicit status. A pair that is absent from the table means no evidence, not homozygous reference.

## Implementation
- `src/python/extract_genotypes.py` + `EXTRACT_GENOTYPES` (`modules/local/genotypes.nf`): streams the merged VCF through `bcftools query` (the `datasci` env has pyarrow but not pysam, so pysam is not used), only FORMAT fields declared in the header. No file for inputs without samples or without GT. `gt_status`: `hom_ref`, `het`, `hom_alt`, `no_call`, `partial_no_call`, `haploid_ref`, `haploid_alt`, `multiallelic_other`, `alt_plus_other_allele` (unphased `1/0`). (The planned `gt_absent` status was dropped: `bcftools query` prints `.` for a missing GT, which is `no_call`, and inputs without a GT field write no table.) The input mode is derived from the VCF header by the script itself, not read from the check report.
- `src/python/genotypes_to_wide.py` + `ADD_WIDE_GENOTYPES` (`process_medium`): only with `--genotype_wide true`, Parquet only; streams the main table in batches and writes a NEW `<run_id>.parsed.clean.withGT.pq`; exit 1 above `genotype_wide_max_samples` (10), if the row count changes, or if any table row has no genotype row. Reason for batching: the real S223 table has 22,624,536 rows x 591 columns (5.0 GB, 22 row groups), far beyond what one pandas frame holds.
- `write_run_manifest.py completed` now records the input check (genotype mode, samples, normalisation counts) and the genotype file identity.
- Params in `nextflow.config`: `genotypes`, `genotype_wide`, `genotype_wide_max_samples`, script paths.

## Finding: unphased `1/0` is not a heterozygote (question 4)
In the S223 chr22 slice, 1,537 of 88,323 genotypes are `1/0`. All 1,537 have DP > AD ref + alt (reads for another allele) and 80.6% have zero reference reads; of the 63,069 `0/1` genotypes none has either property. Callers write heterozygotes `0/1`; `bcftools norm -m -any` writes `1/0` when the sample carries this ALT plus a different ALT (a spanning-deletion `*` ALT is one such other allele; 230 `*` records exist in the slice). These rows are therefore labelled `alt_plus_other_allele`, and map to `HET` in the wide view (one copy of this ALT). Whether Module 3 should treat them differently is for Mario to decide. Siblings at the same position are not a usable signal (none in the slice).

## Checks run (local)
- 64 tests pass (`test_extract_genotypes.py`, `test_check_vcf_assembly.py`, `test_write_run_manifest*.py`).
- chr22 of the real `S223.annotated.vcf.gz` (88,323 records, 1 sample): 88,323 rows; Locus, GT, GQ, DP and AD identical to `bcftools query`; 4.6 s. Statuses (before the `alt_plus_other_allele` split): 64,606 `het`, 23,717 `hom_alt` (a variant-only VCF has no `hom_ref` rows, as expected).
- Wide view on real data: the 385,038 chr22 rows (591 columns) of `S223.parsed.clean.pq` joined to the chr22 genotype table, exit 0 (no row without genotypes), 26 s in batches of 250,000 rows; extrapolated to 22.6 M rows about 25 min (not measured; memory not measured either).
- `Locus` format against the real table: all 88,093 distinct `chr22:` Locus values of `S223.parsed.clean.pq` are in the genotype table (0 missing); the genotype table has 230 more (ALT `*` records the main table does not contain).
- `nextflow run annotate_vcf.nf -profile local_dev -stub-run --genotype_wide true` on `panel7_refcorrected.vcf.gz`: all 19 tasks, including the new processes, completed (stub only: the stub of `EXTRACT_GENOTYPES` writes an empty file, the real script writes none for sites-only input).
- The stub run also exposed that bcftools 1.24 (the `local_dev` build) prints a longer `Lines` summary than 1.21 (cluster); the input check failed on it. Fixed in `parse_norm_summary`, tested against both formats.

## Not done / open
Cluster run (panel7, then the <family> family); Module 3 notification and a `DEC-` record; cohort-scale inputs; callset merging, coverage-based upgrade of no-calls, pedigree and person mapping (outside this pipeline by decision).
