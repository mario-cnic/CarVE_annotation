# Module 3 handoff: table contract changes of the 2026-11-04/05 review

Written 2026-11-05T08:00:00+01:00 in `CarVE_annotation` (Module 2). Module 3 (`clinical_variant_prioritization`) and `carve-platform` are not in this container: neither was read nor edited, so nothing here is checked against Module 3's code. Labels: FACT (code or test read now), INFERENCE, PROPOSAL.

**2026-11-05T17:30+01:00:** the three branches are merged and pushed (`main` `1dc8c74`) and cluster-checked (`TODO.md` branch checklist). The draft below became `carve-platform/decisions/DEC-0013-labranchor-informational.md` (PROPOSED, not committed there). Module 3 was read for it at `e87b3a1`: its own `NEW_IMPACT` still has the position-only rule (`acmg.py:233-235`, CLI only), and its default-on "Branchpointer" predictor counts `Branchpoint_disrupted == "YES"` (`predictor_registry.py:61-62`).

## Scope

- Three branches, none merged into `main`, none run on the cluster:
  - `ccr-8536b091-fabr1g`: the fixes (48 commits over `main` when this was written).
  - `pred18-labranchor-positional`: `ccr` plus PRED-18 (1 commit).
  - `labranchor-ism`: `pred18` plus PRED-20 and the ISM columns (2 commits).
- Existing tables are unchanged. Everything below applies to tables built from these branches.
- Not covered: the cluster checks for these branches (`TODO.md`, item "Cluster checks for branch `ccr-8536b091-fabr1g`").

## A. Columns removed or added

| Change | Columns | ID | Branch | Module 3 action (PROPOSAL) |
| --- | --- | --- | --- | --- |
| Removed | `Branchpointer_prob`, `Branchpointer_U2_energy`, `Branchpoint_disrupted`, `Branchpoint_status` | PRED-18 | `pred18-labranchor-positional` | Stop reading them. FACT: Branchpointer was never run. `_prob` was a copy of `LaBranchoR_score`, `_U2_energy` a constant (-6.2), and `Branchpoint_status` was never `disrupted`. No information is lost. |
| Added | `LaBranchoR_ISM_d_best`, `LaBranchoR_ISM_d_max` (float), `LaBranchoR_ISM_dist`, `LaBranchoR_ISM_n_3pss` (integer), `LaBranchoR_ISM_status` (`scored` / `not_covered` / `not_applicable`) | PRED-18 follow-up | `labranchor-ism` | Display only. Do not use them in tiering until a threshold is chosen (none validated). See the notes below the table. |

Notes on the ISM columns:
- They are variant-level: every transcript row of a variant gets the same values, from the worst 3' splice site. They are not matched to the row's transcript.
- Only SNVs get values. Indels are `not_applicable`.
- The windows come from GENCODE v19 introns lifted to GRCh38.

Kept, with corrected meaning (PRED-18):
- `LaBranchoR_score` is the probability of LaBranchoR's top branch point within 2 bp of the variant. It describes the site; the alleles are not scored.
- `LaBranchoR_acc_dist` is the distance from the variant to that branch point's 3' splice site. The old dictionary said "from the branch point".
  2026-11-05T15:45+01:00: corrected by PRED-21. It counted one base too few on + and one too many on −; it now counts like `LaBranchoR_ISM_dist` and the HGVS intron offset (1 = last intron base). Values in tables built before this change are off by one.
- `LaBranchoR_status` stays `scored` / `not_covered`.

## B. Value changes in existing columns

| Column(s) | Before | After | ID | Branch |
| --- | --- | --- | --- | --- |
| `NEW_IMPACT`, then `PRIORITY_TIER` | `LaBranchoR_score >= 0.50` alone set `NEW_IMPACT` to HIGH | That rule is removed. Such rows keep their VEP or splicing-derived impact. `PRIORITY_TIER` can drop where that HIGH was the row's only functional signal (`scoring.py` `has_functional_signal` reads `NEW_IMPACT` MODERATE/HIGH). `VARIANT_PRIORITY_SCORE` does not read `NEW_IMPACT`. | PRED-18 | `pred18` |
| `HGVSc` | Text `"nan"` on rows without HGVSc | Null | MISC-33 | `ccr` |
| `intron_offset_signed`, `splice_side` | Only `c.<n>±k` was parsed. The rest fell back to `HGVS_OFFSET`, which is VEP's 3' shift of the notation, not an intron offset. Example: `c.-45+3` came out as acceptor -3, and 3'-shifted exonic indels as donor/acceptor. | 5'UTR, 3'UTR, non-coding (`n.`) and range HGVS are parsed. Without an intronic position the row is `0` / `exonic`. | PRED-8 | `ccr` |
| `SpliceVault_Predictions_Decoded` | `pos:` counted from the variant `POS` (0 when missing). Offset 0 read "at mutation site". | `pos:` = `SpliceVault_site_pos` ± offset by strand, omitted when the site is missing. Offset 0 reads "the annotated splice site". `SpliceVault_status` is unchanged. | PRED-7 | `ccr` |
| `SCORE_EXPLANATION`, `TIER_RATIONALE` | A PVS1-vetoed LoF (MYH7) was described as `LoF (+35.0)` / "Loss-of-Function" | Only active (non-vetoed) LoF is described as LoF, matching the score, which was already right | ACMG-7 | `ccr` |
| Disease evidence (`DISEASE_EVIDENCE_SCORE` and dependent score/tier) for `AGL` | `AGGL` (not a gene) in the storage/infiltrative core panel, so it never matched | `AGL`. FACT: rows change only if `AGL` is in `Complete_gene_list_V4.csv` (cluster file, not checked) | ACMG-9 | `ccr` |
| Categorical columns (`NEW_IMPACT`, `PRIORITY_TIER`, ...) in chunked tables | If the first chunk had the column all null, its categories and their order came from the values seen | Same categories and order as an unchunked table. Compare category values, not codes. | MISC-37 | `ccr` |
| `SpliceAI*`, `Pangolin*`, `SPiP*` on symbolic, breakend and `*` records | Pangolin scored `<CNV>`, breakends and similar as a literal-text insertion; SpliceAI did not skip breakends | Not sent to these predictors. The row has no score from them. | PRED-17 | `ccr` |

## C. Rows

- **MISC-28** (`ccr`): a VCF record without CSQ used to be missing from the table. It now gets one row with empty VEP fields.
  - INFERENCE: such a row has no `Feature` and no `MANE_SELECT`, so `TRANSCRIPT_PRIORITY_TIER` = 3 (`tag_transcript_priority.py:84-91`), and it is outside a tier 1-2 default view.
  - How many records VEP leaves without CSQ has not been measured. The run QC report's `completeness` row (MISC-2) now compares the table's distinct `Locus` with the merged VCF.
- **MISC-29, MISC-36** (`ccr`): a parse error on a record, or in the SpliceVault/MaxEntScan step, now fails the task. Before, the record or those columns went missing and the run still exited 0.
  - INFERENCE: tables built before the fix could be missing records or those columns, with no sign of it. MISC-2's completeness check covers the records on new runs.

## D. Run, input and provenance contract

- **Input refusals** (`ccr`). Neither case reaches Module 3; both now fail the run before any predictor starts.
  - MISC-34: the input already declares a predictor INFO tag (`SpliceAI`, `Pangolin`, `SPiP`, `LaBranchoR_*`, ...). The old value would have survived the merge.
  - MISC-39: a `<...>` ALT without an `INFO/END` header line.
  - INFERENCE: Module 1 (Sarek) VCFs carry neither.
- **`<table>.provenance.json`** (MISC-38, `ccr`): when a later launch with the same run id republishes a table, the old sidecar is moved to `<table>.provenance.superseded.<timestamp>.json` and a new one is written. Before, the old one stayed and described the earlier run.
  - PROPOSAL: Module 3 reads `<table>.provenance.json` only and ignores `*.superseded.*`.
- **`RUN_MANIFEST.json`** (`ccr`):
  - MISC-42: values now come from the params of the actual command line. Before, a `--fasta`, `--spliceai_distance` or other override on the command line was recorded as the config value, and the score-cache version id came from it too.
  - ACMG-6: `external_provenance.resources.disease_gene_list` (path, size, sha256) is new.
- **Gene-tier coverage** (MISC-4, `ccr`, rebuild pending on the cluster): once `resources/v5_genes_loc.bed` is rebuilt, SpliceAI/Pangolin/SPiP cover genes that were missing before, such as `TAZ`, mitochondrial genes, and `SLC37A4`/`SHANK3` (MANE has them on alt/fix contigs only). A `not_covered` status in older tables can mean "outside the BED" for those genes.
- **SpliceAI errors** (PRED-16, `ccr`): an unexpected SpliceAI exception now fails the chunk. Before, the record was written unscored and read as `not_covered`.
  - FACT: there is no evidence it ever happened.
  - The score cache's scored sites are migrated to the new code version (step in `TODO.md`); unscored ones are recomputed.
- **`-resume`** (MISC-41, `ccr`): a resumed run now reruns every step whose code changed.
  - The first resume on this branch reruns every step once, because each script block gained a `# code` line.
  - Score-cache version ids are unchanged: FACT for SpliceAI and Pangolin, checked with `annotation_cache.py`.

## E. No change to the table

These change no table column, row or value:
- PRED-2 and MISC-27: CLI options removed; the pipeline passed neither.
- MISC-30/31/32: parser genotype columns that VCF_TO_TSV does not request.
- ACMG-4 and ACMG-8: one definition instead of two, output identical.
- MISC-26: pandas metadata in the `withGT` table; no visible change on read.
- MISC-6: identity checks.
- MISC-44: resource identity in the `-resume` hash.
- MISC-40: local profile.
- LIFT-1/2: standalone liftover tool.
- PRED-19, MISC-22.
- MISC-43: `parse_missense` no longer shrinks `MISSENSE_COLUMNS`; one TSV per process, so pipeline tables were not affected (INFERENCE).

## Draft decision record (PROPOSAL, for `carve-platform`; number assigned there)

> **DEC-XXXX: Branch-point and parser fixes change the Module 2 table contract (2026-11-05)**
>
> Context: the Module 2 review of 2026-11-04/05 found that the "Branchpointer" columns were a positional LaBranchoR lookup that could raise `NEW_IMPACT` to HIGH from position alone (PRED-18), and several parser defects that dropped records or wrote wrong values (MISC-28/29/33/36, PRED-7/8).
>
> Decision: Module 2 removes `Branchpointer_prob`, `Branchpointer_U2_energy`, `Branchpoint_disrupted` and `Branchpoint_status`. `LaBranchoR_score` / `_acc_dist` / `_status` stay as informational, site-level columns with no effect on `NEW_IMPACT`, score or tier. New informational `LaBranchoR_ISM_*` columns give an allele-aware SNV effect from the LaBranchoR authors' in silico mutagenesis, not used in tiering. The value fixes in sections B-C of `CarVE_annotation/walkthrough/20261105_module3_table_contract_changes.md` apply to every table built from the merge commit onwards.
>
> Consequences for Module 3:
> - Stop reading the four removed columns.
> - Treat `HGVSc == "nan"` as missing in older tables.
> - Expect tier-3 rows without VEP fields.
> - Expect some rows to lose a LaBranchoR-only HIGH `NEW_IMPACT`, and their `PRIORITY_TIER` to drop with it.
> - Read only `<table>.provenance.json`, never `*.superseded.*`.
> - Old and new tables differ in these columns; do not compare them across the merge date without checking the table's commit (`carve.provenance` footer or the provenance sidecar).

## Open for Mario

- Merge order: `ccr` first, then `pred18`, then `labranchor-ism`. They stack, so `labranchor-ism` contains all three.
- Whether this needs a `DEC-` (CLAUDE.md: only if another project must know).
  - INFERENCE: yes for PRED-18, since columns are removed and `NEW_IMPACT` changes.
  - The value fixes alone could go in `STATUS.md`.
- A threshold for the ISM columns, if they are ever to affect tiering.

## Not done

- Nothing was changed in Module 3 or `carve-platform`.
- Module 3's reads of the removed columns were not checked (repo not available here).
  2026-11-05T17:30+01:00: checked; listed in DEC-0013 (Decision item 2, Evidence). Nothing changed in Module 3.
- No cluster run of any branch.
