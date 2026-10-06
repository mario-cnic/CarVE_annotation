# Working with the 6 GB family table (2 WES + 1 WGS): measured funnel and plan

Written 2026-10-05T16:00:00+02:00. Scope (Mario, 2026-10-05): the target is tables like `F2005_union` (3 samples, 6 GB), not 100 WGS. Labels: FACT (measured now), INFERENCE, PROPOSAL. Only aggregate counts are kept here; sample names stay in the family's own project folder.

## FACT: the table and DuckDB

`nf_work/annotation_out/F2005_union/F2005_union.parsed.clean.pq`: 23,686,438 rows (variant x VEP transcript), 5,244,357 distinct `Locus`, 614 columns, 6.0 GB; `F2005_union.genotypes.pq` 15.8 M rows (records x 3 samples). Machine: 12 cores, 62 GB RAM (17 GB in use).

DuckDB 1.5.5 over the Parquet in place (`memory_limit` 10 GB, 8 threads, nothing loaded into pandas), `src/tools/candidate_funnel_duckdb.py`:

| Query | Time | Peak RSS |
|---|---|---|
| row and distinct-Locus count | 1.2 s | 0.5 GB |
| the whole funnel below (8 counts) | 4.6 s | 0.7 GB |
| genotype statuses, carriers per sample, carriers per locus | 0.2-0.4 s each | 0.7 GB |
| reference row per (Locus, SYMBOL), window over 3 columns (6,643,027 rows) | 1.9 s | 1.1 GB |
| grid page, top 50 candidate rows, 17 columns, straight from the file | 0.4 s | 1.1 GB |
| all rows of one gene | 0.1 s | 1.1 GB |

So the file is not the problem for DuckDB. Module 3's upload path is (`variants.py:72` `await file.read()`, `:83` `pd.read_parquet`, scoring in pandas before the Parquet write). INFERENCE: Module 3's own benchmark (5.97 s of tiering and scoring per 10,000 rows, `TODO.md` `BUG-14`) extrapolates linearly to about 3.9 h for 23.7 M rows, before memory; it was not measured at this size.

## FACT: funnel (rows, distinct Locus); thresholds are illustrative, not validated for cardiomyopathy

| Step | Rows | Loci |
|---|---|---|
| all | 23,686,438 | 5,244,357 |
| transcript tier <= 2 | 2,183,065 | 2,086,558 |
| rare (gnomAD v4 joint grpmax AF <= 0.001 or missing) | 830,185 | 253,206 |
| rare and VEP IMPACT HIGH/MODERATE | 2,247 | 734 |
| rare and (SpliceAI custom or Pangolin) >= 0.2 | 1,371 | 140 |
| ClinVar pathogenic term in `CLIN_SIG`, any AF (old source, no stars) | 1,980 | 178 |
| **candidate** = rare and (impact or splice), or ClinVar pathogenic | **5,402** | **1,020** |
| candidate, tier <= 2 rows only | 906 | 871 |

With Module 2's `NEW_IMPACT` instead of `IMPACT` the candidate set is 6,072 rows. Scoring 5,402 rows is seconds with the existing pandas engine (INFERENCE from the 5.97 s per 10,000 rows figure).

Genotypes (per sample, status counts of 5,263,518 records): the two WES samples have about 4.25 M `no_call` each (81% of records, the WGS-only sites) and 583,243 / 512,987 `hom_ref`; the WGS sample has 177,541 `no_call` and no `hom_ref` (variant-only VCF). Candidate loci carried (het, hom_alt or alt_plus_other_allele): 511 / 513 / 540 per sample; 170 by all three, 204 by two, 646 by one.

## Pitfall found while measuring (FACT)

A first version of the script used `select *, gnomADv4_AF_grpmax_joint as af`. The table already has the 1000 Genomes column `AF`; DuckDB treats names case-insensitively, `a.af` resolved to that column, and the "rare" count came out 5.5x too high (4,551,824 rows instead of 830,185) without any error. Same class as `PLT-021`. Any SQL that Module 3 writes over these tables must name columns explicitly and never alias to a name that differs from an existing column only by case.

## PROPOSAL (Module 3 work; no Module 2 change needed at this scale)

1. **Ingest by path.** Register the run's Parquet and genotype table and query them in place with DuckDB; no upload, no pandas load, no rewrite.
2. **Two layers.** The full table stays queryable and unscored. A candidate set is selected in SQL (flags: rare, impact, splice, ClinVar, transcript tier; thresholds in one versioned, user-visible config), materialised to a small Parquet, and only that set goes through the existing scoring, ACMG and pedigree code. Result: seconds, not hours, and the engine is reused.
3. **Safety nets.** Record the candidate definition and its hash with the result; show how many rows each criterion removed; keep "include missing frequency" and "ClinVar pathogenic regardless of frequency" as explicit clauses (as above); an "expand" control re-runs the selection with looser thresholds, with a row-count warning above a cap (INFERENCE: a few hundred thousand rows is minutes).
4. **Family logic from the genotype table in SQL** (carrier per sample, segregation, de novo). `no_call` is unknown, not homozygous reference: in this table a WES `no_call` at a WGS-only site is missing coverage. Compound-het partners must come from the same candidate selection (a rare low-impact partner outside the filter is not seen).
5. **Not now:** the sites / variant_gene / variant_transcript split (`PLT-130`) and flag columns written by Module 2; they would slim the table but are not needed at 3 samples. Sparse genotype storage and 100-WGS benchmarks stay as design targets (D-6 still open).

## Not verified

Thresholds were not tuned on reviewed cases (recall of the candidate set on known variants is unmeasured); the Module 3 pipeline was not run on this table; timings are one machine, warm file cache not controlled.

## Addendum 2026-10-05T18:30:00+02:00: layout and the cost of fetching whole rows (answers a Module 3 report)

Module 3 reported that fetching the candidate rows with all 614 columns is slow (33-38 s warm, 106 s cold, 4-5 GB) and that row-group statistics prune almost nothing. Checked on `F2005_union` (DuckDB, in place, 8 threads, 12 GB limit):

- **FACT, layout:** the table is not in genomic order. Every row group spans `chr10` to `chrY` by its `Locus` statistics and the chromosome changes at nearly every sampled row. Cause: `src/python/parsing/filter_variants.py:429-433` sorts the table by `NEW_IMPACT`, then `SPiP_prediction`, then the frequency column (`io_qc.py:402-422`, pandas `sort_values`). Row groups hold 1,048,576 rows (23 groups).
- **FACT, the sort helps one predicate:** all 210,359 rows with `NEW_IMPACT` HIGH, MODERATE or LOW come first (HIGH rows 0-28,626, MODERATE to 115,652, LOW to 210,358), so they sit in row group 0 and the other 22 groups are pruned by the `NEW_IMPACT` statistics.
- **FACT, fetch cost by number of columns** (1,020 candidate loci, 9,799 rows over all transcripts, no `NEW_IMPACT` predicate): 17 columns 2.8 s; 60 columns 3.6 s; 150 columns 12.9 s; 300 columns 17.9 s; 612 columns 53 s (4.7 GB). Roughly linear in the number of columns.
- **FACT, full-column fetch of the 5,402-row candidate set (impact, splice, ClinVar definition):** 85.9 s without a `NEW_IMPACT` predicate; 7.3 s when `NEW_IMPACT IN (HIGH, MODERATE, LOW)` is added (4,091 rows). The missing 1,311 rows are `MODIFIER` rows (splice score >= 0.2 or ClinVar pathogenic via `CLIN_SIG`), which Module 2's `NEW_IMPACT` does not upgrade. Using `NEW_IMPACT` as the candidate rule would lose 24% of this set.
- **FACT, a permissive superset is tiny:** rare (AF <= 0.01 or missing) and (IMPACT HIGH/MODERATE/LOW or splice >= 0.1), or a ClinVar pathogenic term: 14,533 rows (3,087 loci), 0.06% of the table.
- **FACT, re-sorting the finished file is not practical:** an attempt to write a copy ordered by that superset flag with a single DuckDB `COPY ... ORDER BY` over all 614 columns was stopped after 40 min (`timeout`) with 84 GB of spill files; nothing was produced. The existing pandas sort in `filter_variants.py` already orders the whole frame, so changing its keys costs nothing extra at write time.
- **INFERENCE (not measured):** if Module 2 sorted by a permissive superset flag first, the 14.5 k superset rows would sit in one row group and every candidate definition inside the superset would read about one row group (the pruned fetch above took 7.3 s).

PROPOSAL, in order: (1) Module 3 fetches only the columns the engine and grid read (60-150 columns: 4-13 s instead of 53-86 s) and builds the candidate set once per definition, cached; no Module 2 change. (2) Optional, small, needs a `DEC-` record because it changes the file's row order: Module 2 sorts by a deliberate key (permissive superset flag, then `NEW_IMPACT`) instead of the current accidental one; applies to new runs only; genomic ordering inside a stratum waits for `chrom`/`pos` columns (`PLT-023`). (3) The sites / variant_gene / variant_transcript split (`PLT-130`, 4.46 rows per variant today) is the structural fix, not needed at 3 samples.
