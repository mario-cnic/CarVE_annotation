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
- `src/python/extract_genotypes.py` + `EXTRACT_GENOTYPES` (`modules/local/genotypes.nf`): streams the merged VCF through `bcftools query` (the `datasci` env has pyarrow but not pysam, so pysam is not used), only FORMAT fields declared in the header. No file for inputs without samples or without GT. `gt_status`: `hom_ref`, `het`, `hom_alt`, `no_call`, `partial_no_call`, `haploid_ref`, `haploid_alt`, `multiallelic_other`. (The planned `gt_absent` status was dropped: `bcftools query` prints `.` for a missing GT, which is `no_call`, and inputs without a GT field write no table.) The input mode is derived from the VCF header by the script itself, not read from the check report.
- `src/python/genotypes_to_wide.py` + `ADD_WIDE_GENOTYPES`: only with `--genotype_wide true`; writes a NEW `<run_id>.parsed.clean.withGT.<fmt>`; fails above `genotype_wide_max_samples` (10) or if the join changes the row count.
- `write_run_manifest.py completed` now records the input check (genotype mode, samples, normalisation counts) and the genotype file identity.
- Params in `nextflow.config`: `genotypes`, `genotype_wide`, `genotype_wide_max_samples`, script paths.

## Checks run (local)
- 64 tests pass (`test_extract_genotypes.py`, `test_check_vcf_assembly.py`, `test_write_run_manifest*.py`).
- chr22 of the real `S223.annotated.vcf.gz` (88,323 records, 1 sample): 88,323 rows; Locus, GT, GQ, DP and AD identical to `bcftools query`; 4.6 s. Statuses: 64,606 `het`, 23,717 `hom_alt` (a variant-only VCF has no `hom_ref` rows, as expected).
- `nextflow run annotate_vcf.nf -profile local_dev -stub-run --genotype_wide true` on `panel7_refcorrected.vcf.gz`: all 19 tasks, including the new processes, completed (stub only: the stub of `EXTRACT_GENOTYPES` writes an empty file, the real script writes none for sites-only input).
- The stub run also exposed that bcftools 1.24 (the `local_dev` build) prints a longer `Lines` summary than 1.21 (cluster); the input check failed on it. Fixed in `parse_norm_summary`, tested against both formats.

## Not done / open
Cluster run (panel7, then the <family> family); Module 3 notification and a `DEC-` record; cohort-scale inputs; callset merging, coverage-based upgrade of no-calls, pedigree and person mapping (outside this pipeline by decision).
