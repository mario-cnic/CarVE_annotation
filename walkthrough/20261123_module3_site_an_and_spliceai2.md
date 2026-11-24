# Module 3 handoff: gnomAD site AN for PM2 (SF-02) and SpliceAI2 columns (PR-03)

Written 2026-11-23T11:00:00+01:00 in `CarVE_annotation` (Module 2). Module 3 (`clinical_variant_prioritization`) was not read or edited for this note, so nothing here is checked against Module 3's code. Labels: FACT (code, table or run read now), INFERENCE, PROPOSAL.

Applies to tables built from `main` at `0562b5c` (SF-02) and `01e5881` (PR-03) onwards. Older tables are unchanged.

## A. Changed meaning: `absent_from_gnomad` and PM2 (SF-02)

| | Before | After |
| --- | --- | --- |
| Allele with no gnomAD v4 record | always `absent_from_gnomad`, PM2_Supporting, counts as rare | `absent_from_gnomad` only if `GNOMAD_SITE_AN` ≥ 74,894; otherwise `absent_low_an` |
| `absent_low_an` | did not exist | no PM2, no rarity bonus; treated like an unknown AF; kept by rarity filters for review |

- FACT: the threshold is −ln(0.05)/`pm2_absence_af`, i.e. 74,894 for the default 4e-5 (`nextflow.config` `pm2_absence_af`, parser `--pm2-absence-af`). With 0 observed alleles in that many, the AF is below 4e-5 with 95% confidence (Cardiomyopathy VCEP MYH7 v2.0 rule).
- FACT: `GNOMAD_SITE_AN` (new, numeric) is the gnomAD v4.1 *allele number at all callable sites*, exomes + genomes, minimum over the positions an indel spans, 0 where gnomAD has no callable genotype. Raw VEP columns `gnomADv4AN_exomes` / `gnomADv4AN_genomes` are also in the table (`&`-joined per position for indels).
- INFERENCE: rows that lose PM2 can lose priority score (15 points in the parser test `test_absence_needs_enough_site_an_for_pm2`) and tier. Expected mostly outside exome capture where only genomes count (genomes max AN about 152k, so a genomes-only site stays above the threshold unless coverage is poor), on chrY, and in poorly genotyped regions.
- Present alleles are not affected: `available_frequency` keeps its PM2 logic on the AF.
- Known limit: for a deletion, the BED overlap covers the deleted bases plus the anchor; for an insertion, its two flanking bases. Positions missing from the gnomAD file count as AN 0.

## B. New columns: SpliceAI2 (PR-03)

| Column | Type | Meaning |
| --- | --- | --- |
| `SpliceAI2_score` | float | max of the four delta scores (0-1) |
| `SpliceAI2_DS_AG`, `_DS_AL`, `_DS_DG`, `_DS_DL` | float | delta score of the strongest acceptor/donor gain/loss |
| `SpliceAI2_DP_AG`, `_DP_AL`, `_DP_DG`, `_DP_DL` | integer values (stored as float, null where not covered) | position of that site relative to the variant |
| `SpliceAI2_strand` | string | strand of the gene scored; with genes on both strands the higher-scoring row |
| `SpliceAI2_status` | string | `scored` / `not_covered` |

- FACT: exact-allele lookup in Illumina's GRCh38 precomputed scores (v2.0, Hugging Face revision `cd1b8b8`): every SNV and the population-observed indels inside GENCODE v48 protein-coding genes. Variant-level: every transcript row of a variant gets the same values.
- FACT: `not_covered` has no score, never 0. A real 0 is `scored` with `SpliceAI2_score` 0.
- FACT: informational only. No change to `SpliceAI*`, `SPLICE_MAX_UNIFIED`, `NEW_IMPACT`, ACMG criteria or priority score.
- PROPOSAL: Module 3 shows them for display and comparison; no use in tiering (no ClinGen calibration; author thresholds 0.1 / 0.25 / 0.5).
- Licence (`PR-02`, Mario's answer 2026-11-10): academic, non-commercial research. Showing the values in reports to families and in tables shared outside the lab is covered: YES. No redistribution of the score files without Illumina's consent.

## C. Run checks

- FACT: cluster fixture `an_spliceai2_fixture` (2026-11-23, commit `76478fd`, tree clean, 33/33 tasks completed), input `test_data/raw_vcfs/an_spliceai2_fixture.vcf.gz`, 9 synthetic variants, 174 table rows. Every variant matched `test_data/raw_vcfs/an_spliceai2_fixture.expected.tsv` in `GNOMAD_SITE_AN`, `POPULATION_AF_STATUS`, `SpliceAI2_status` / `_score` / `_strand`, and PM2 followed the status:
  - `absent_from_gnomad` (4, PM2_Supporting): site AN 1,613,700; 152,376; 152,332 (private insertion, flanks); 152,246 (intergenic).
  - `absent_low_an` (4, no PM2, `Unclassified`): AN 0 (chr21:16290571), 73,578 (SNV), 71,016 (deletion, minimum of three positions), 34,437 (chrY).
  - `available_frequency` control (AF 1.7e-6): PM2 from the AF, unchanged.
  - SpliceAI2: 3 `scored` (two at 0, one at 0.02 on a two-strand site, − strand kept), 6 `not_covered` (private indel, intergenic, outside protein-coding genes).
- Not yet measured: how many rows of a real cohort move from `absent_from_gnomad` to `absent_low_an`.

## Draft decision record (PROPOSAL, for `carve-platform`; number assigned there)

> **DEC-XXXX: PM2 absence needs gnomAD site allele number; SpliceAI2 added as informational (2026-11-23)**
>
> Context: an allele absent from gnomAD earned PM2_Supporting at any site, including sites gnomAD barely genotyped, where absence says little (SF-02). SpliceAI2 precomputed scores became available under an academic licence (PR-02/PR-03).
>
> Decision: Module 2 reports `absent_from_gnomad` only where gnomAD v4.1 exomes + genomes site AN ≥ −ln(0.05)/4e-5 = 74,894; below that the new status `absent_low_an` earns no PM2. New column `GNOMAD_SITE_AN`. New informational `SpliceAI2_*` columns, not used in impact, ACMG or score.
>
> Consequences for Module 3:
> - Accept the new `POPULATION_AF_STATUS` value `absent_low_an` and treat it like an unknown frequency.
> - If Module 3 derives PM2 itself from the gnomAD columns, apply the same site-AN rule or read Module 2's `ACMG_CRITERIA`, so both modules agree.
> - Show `SpliceAI2_*` as display-only; `not_covered` is not 0.
> - PM2, score and tier of absent variants differ across the change date; check the table's commit before comparing.

## Resolved 2026-11-24

- Filed as `carve-platform/decisions/DEC-0015-pm2-site-an-and-spliceai2.md` (PROPOSED, not committed there).
- FACT: Module 3 (`b725eaf`) computes PM2 itself, in both engines, from `POPULATION_AF_STATUS` only. Any status outside its four known values becomes `missing_frequency` (`modules/io_qc.py:68,80`, `engine2/modules/frequency.py:58-59`), so `absent_low_an` already gets no PM2 and no BA1/BS1: the modules agree. Shown as "frequency missing", and engine2's PM2 detail still says Module 2 does not emit allele number. The fixes are listed in DEC-0015 (Consequences) as a Module 3 task.
