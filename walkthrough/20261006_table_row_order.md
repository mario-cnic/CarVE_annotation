# Table row order: coordinate order instead of the legacy impact sort

Written 2026-10-06T09:20:00+02:00. Decision (Mario, 2026-10-06): the legacy sort by `NEW_IMPACT` does not make sense for a table that a query engine reads. Labels: FACT (measured or read now), INFERENCE, PROPOSAL.

## Why (FACT, from `20261005_family_table_candidate_funnel.md` and checks on 2026-10-06)

- `filter_variants.py:429-433` sorted the final table by `NEW_IMPACT`, `SPiP_prediction` (descending) and the frequency column, missing values first. In `F2005_union` every row group therefore spans `chr10` to `chrY`; only a predicate on `NEW_IMPACT` itself could prune row groups, and that misses 24% of a splice / ClinVar based candidate set.
- Nothing in Module 2 reads the order; Module 3's report tables sort by their own `VARIANT_PRIORITY_SCORE` before `head(250)`.
- The order was not a priority order either: inside the `HIGH` block the first 1,000 rows have no SPiP score (missing first), the last 1,000 all have one.
- The sort ran in pandas over the whole frame (23.7 M rows x 614 columns), a memory cost on top of the parse.

## What changed

- `src/python/parsing/filter_variants.py`: the sort is removed. Before `CHROM`/`POS` are dropped, `check_coordinate_order` (new, `modules/io_qc.py`) checks that each contig is one contiguous block and `POS` never decreases inside it; the result is logged (info when ordered, warning with the counts when not). Nothing fails on unordered input: the data is intact, only the order contract is broken. The contig order is the input VCF's.
- `tests/python/test_row_order.py` (6 tests; `io_qc.py` is loaded by path because other tests put a different `modules` package on `sys.path`).
- `src/python/parsing/PROVENANCE.md`: divergence recorded.

## Checks (local)

- Real fixture: `panel8_trio.annotated.vcf.gz` through the parser (stage 1) and `filter_variants.py` (stage 2), before and after. 1,015 rows x 594 columns both times, same columns, identical rows and values once both are sorted by `Locus`,`Feature`; the output is coordinate-sorted within every contig (7 contigs) and the log says "Rows are in coordinate order (7 contigs)".
- Suite: 141 passed (without `test_predictor_attribution.py`, which fails to collect in a full run because of the shared `modules` import; it passes alone, 25 tests, before and after).

## Not verified

- Not run on a WGS-size input or on the cluster; for a full VCF the order relies on `bcftools sort` earlier in the pipeline (`CONCAT_VEP`), which the check will report on.
- The gain for gene / region / locus fetches is expected, not measured. It needs numeric `chrom` / `pos` columns (`PLT-023`) before row-group statistics can prune by position; with only the string `Locus` they cannot.
- Existing tables (S223, F2005_union) keep the old order until re-run.
- The file's row order is a contract with Module 3: a `DEC-` record in `carve-platform` is needed (not written from here). Module 3 should not rely on `NEW_IMPACT` clustering (it no longer exists in new tables).

## Decision-record follow-up — 2026-10-08T10:26:41+02:00

The missing cross-module record is now registered as [DEC-0010](../../carve-platform/decisions/DEC-0010-table-row-order.md). The approved contract preserves input order rather than guaranteeing a new sort: coordinate-sorted input remains coordinate-sorted; unordered input is warned about. Cluster validation and Module 3 verification remain open.
