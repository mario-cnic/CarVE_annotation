# Predictor score attribution: wrong-gene / wrong-transcript scores on overlapping loci

**2026-10-01.** Trigger: the Pangolin re-test on `panel7` returned `MADD` scores for variants the fixture labels `MYBPC3`.
Scope: read-only audit of the real `S223` table + alias-resolution feasibility. No pipeline behaviour changed.

## FACTS

**Why the label said MYBPC3.** `panel7_test.vcf.gz` is synthetic; its `gene=` tag was set by the generator from the padded window
(`v5_genes_loc.bed` = gene body ±10 kb). GENCODE 45: MYBPC3 chr11:47,331,406-47,352,702 (− strand); the 8 variants at
47,328,406-47,330,000 are 1.4-3.0 kb **downstream** (3') of it and inside MADD (47,269,161-47,330,031, + strand). Pangolin correctly scored MADD.

**How scores are attached (shared parser, vendored copy `src/python/parsing/modules/splicing.py`):**
| Predictor | Unit | Parser behaviour |
|---|---|---|
| SpliceAI (`parse_spliceai_custom`) | gene | takes only the **first** entry; never matches row gene |
| Pangolin (`parse_pangolin`) | gene (ENSG) | max \|score\| over **all** genes; gene id never read |
| SPiP (`parse_spip`) | transcript (`NM_`) | first entry whose gene symbol == row SYMBOL, else **first entry of any gene**; transcript never matched |
Rows are variant x VEP transcript; predictor strings are variant-level, so every row of a variant gets the same string.

**Audit** (`src/tools/audit_gene_attribution.py`, S223, 8,993,048 populated rows; rows, not distinct variants; not filtered to curated transcripts):
- SpliceAI (6,233,537 rows): first entry = row gene 85.2%; wrong entry used though the right one was present 81,564 (270 cross the 0.20 floor); row gene not scored, neighbour score attached 576,599 (3,553 with score >= 0.20); row has no symbol 261,471.
- SPiP (7,756,696 rows): row gene absent -> other gene's first entry used 1,237,426 (4,042 >= 0.20); gene has several SPiP transcripts 4,878,242; of 470,073 multi-transcript rows with a MANE row, 301,855 take a non-MANE transcript (136 cross 0.20); MANE not among SPiP entries 34,281.
- Pangolin: S223 predates the db fix (all 2,791,464 rows `ENSMUSG`), so gene matching is not assessable; 128,206 rows multi-gene, 582 of them with the 0.20 outcome depending on which gene.

**Alias vs true neighbour** (`src/tools/audit_alias_resolution.py`; 576,599 SpliceAI "row gene not scored" rows). Resolve each SpliceAI symbol to ENSG
by exon overlap with GENCODE 45 (same chrom+strand), no alias list: alias 187,734 rows (32.6%, 372 distinct pairs), true neighbour 388,434 (67.4%), unresolved 431 (3 symbols).
Spot check of the alias class: PARK2->PRKN, SEPT9->SEPTIN9, IKBKAP->ELP1, FAM134B->RETREG1. Neighbours are mostly overlapping antisense genes (HRG-AS1/KNG1) and include PMS2 rows carrying AIMP2's score.
Cross-check with HGNC complete set (release file dated 2026-09-29, sha256 `91d0ad20...9d889e`; current/prev/alias symbol -> ENSG): coordinate and HGNC classes agree on 10,488 of 11,043 distinct pairs (95.0%; 94.1% of rows);
HGNC has no record for 535 pairs (5.7% of rows; coordinates decide); they **disagree on 19 pairs (1,178 rows, 0.2%)**; coordinates leave 3 pairs unresolved, HGNC resolves 1 of them. Combined: 2 pairs (3 rows) unresolved.

## INFERENCE / limits
- Part of the "neighbour" rows may be alias cases the two methods both miss (HGNC-less symbols judged only by coordinates).
- The 0.20 flips are the rows that can create a false report flag; how many reach tier-1 transcripts (Module 3's input) was not measured.
- One sample (S223). SPiP alias resolution not tested (no SPiP annotation locally; `NM_`->ENSG via HGNC `refseq_accession` + MANE GTF is the proposed route).

## Module 3 relationship (checked)
Module 3 has its **own** copies (`clinical_variant_prioritization/src/filter_variants.py`, `src/modules/splicing.py`), already diverged from the shared ones (shared has `parse_maxentscan`).
Earlier statement in this session that the shared module is "used by Module 3" was wrong. `DEC-0003`: Module 2 owns parsing to the table, Module 3 owns interpretation.

## Actions taken
- Vendored the parsing code into `src/python/parsing/` (commit d5a915c, verbatim, sha256-verified, provenance in its PROVENANCE.md); Nextflow config repointed (fa87d38). Legacy bash (`config/env.sh`, `vcf2parsed.sh`) not repointed.
- HGNC complete set downloaded locally to `~/hgnc_build/` (not deployed anywhere).
- Next: `TODO.md` item "Predictor score attribution".

---

## Addendum 2026-10-01 (later): corrections and implementation

### Corrections to the audit above (found by checking the parser's real behaviour before coding)
- **SpliceAI multi-gene values were not "first entry used".** The old `parse_spliceai_custom` split on `,` but the table joins gene entries with `&`; a multi-gene value therefore failed `float()` and returned NaN for every field (verified: 23,770 of 23,770 multi-gene rows in row group 16 had `spliceai_custom_MAX` = NaN). Those rows fell back to VEP's per-transcript plugin score (`SpliceAI_pred_*`, ~50 bp; verified gene-matched: `SpliceAI_pred_SYMBOL == SYMBOL` on 766,859 of 766,859 rows). So the audit rows "wrong entry used though the right one was present" (81,564; 270 flips) were rows that **lost** their custom score, not rows with a wrong-gene score. The real SpliceAI wrong-gene exposure is single-entry values of another gene (and gene-less rows) - see numbers below.
- The old fallback used `custom > 0`, so a genuine 0.00 custom score was replaced by VEP's value (26 rows, none >= 0.20 in row group 16) - fixed with presence logic.
- The table is biallelic (`Locus` ALT == SpliceAI/SPiP entry ALT on 100% of rows; the `ref`/`alt` columns are empty, ALT is taken from `Locus`), and contains no RefSeq-keyed rows (`Feature` is always `ENS...`, `SOURCE` empty).
- `shared/utils/data/ensembl_to_refseq.tsv.gz` only maps MANE Select / MANE Plus Clinical, so it adds nothing beyond the row's own `MANE_SELECT` / `MANE_PLUS_CLINICAL` columns; not used.

### What was implemented (`ae5a38c`)
Matching rules, statuses, resources and tests are as in the plan; see `src/python/parsing/PROVENANCE.md` for the exact delta. New table columns: `spliceai_custom_match`, `spliceai_custom_anygene_MAX`, `Pangolin_match`, `Pangolin_anygene_max`, `SPiP_match`, `SPiP_anygene_max_prediction`. Existing column names unchanged.

### Results on S223 (8,993,048 populated rows; new parsers applied to the table's raw strings, nothing overwritten; `src/tools/validate_reparse_sample.py`, JSONs in `RUNS/gene_attribution_audit/`)
| | Old parser | New parser |
|---|---|---|
| SpliceAI custom score present | 6,062,036 | 5,585,352 |
| ... removed (other gene's entry / no row gene / ambiguous) | | 622,000 (3,831 had score >= 0.20) |
| ... recovered (multi-gene values the old parser dropped) | | 145,316 |
| ... value changed where both present | | 0 |
| SpliceAI match status | | matched 5,392,563 · alias_resolved 192,789 · no_entry_for_row_gene 470,247 · no_row_gene 171,111 · ambiguous 6,822 · unresolved 5 |
| SPiP present | 7,756,696 | 659,542 matched |
| SPiP status | | not_applicable_transcript 6,912,433 (89%) · no_entry_for_row_gene 184,721 |
| SPiP rows >= 0.20 removed | | 19,453; 142 more crossed 0.20 by choosing the row's own transcript (5,141 values changed) |
| Pangolin | 2,791,464 (all `ENSMUSG`) | 0: the S223 string is from the mouse db, so no entry can match a human `Gene` (expected; Pangolin must be re-run) |

Alias rows cross-checked against the coordinate table: 428 of 192,789 disagree, all from the HGNC-only route (symbol absent from the coordinate table, e.g. `MICALCL`; checked in row group 16: 44/44).
Small test (`test_data/test_run/annotation/MYBPC3.annotated.vcf.gz`, both parser stages run locally with the cluster envs): 238 rows in both outputs, only the matched-predictor columns and columns derived from them change; the 4 MYBPC3 transcript rows that carried `MADD`'s SpliceAI score of 0.30 moved from Tier 3 to Tier 4. Without the identity arguments `filter_variants.py` exits with an error (verified).

### Costs and limits
- SPiP is blank on rows whose transcript has no RefSeq equivalent (MANE Select / Plus Clinical / curated clinical transcript): 89% of S223 rows; tiers 1-2 are covered. Module 3 will see more empty SPiP cells.
- Pangolin rows can only be validated after the re-run with the human db; the logic is tested on synthetic strings only.
- Not validated through Nextflow on the cluster; `HGNC` file must be deployed first (below).

### Note for Module 3 (not written to `carve-platform`; for the user to decide)
Same column names. Values are now only the row's own gene/transcript. Rows without a matching entry are empty with a `*_match` reason, and the variant-level maximum is in `*_anygene_*`. Expect fewer populated SpliceAI/SPiP cells and some priority-tier shifts in this pipeline's own `PRIORITY_TIER`. Module 3's app reads the `.pq` as-is (`pd.read_parquet`) and does not call its own `parse_*`, so it is unaffected mechanically.

### Manual deployment (user)
```bash
mkdir -p /data_lab_PGP/resources/annotation/hgnc
cp ~/hgnc_build/hgnc_complete_set.txt /data_lab_PGP/resources/annotation/hgnc/hgnc_complete_set.2026-09-29.txt
sha256sum /data_lab_PGP/resources/annotation/hgnc/hgnc_complete_set.2026-09-29.txt   # expect 91d0ad20c34f26fb12c8cce3f80a65a6fc694f3a508040b19f6bfa29ae0d889e
chmod a-w /data_lab_PGP/resources/annotation/hgnc/hgnc_complete_set.2026-09-29.txt
```

### Addendum 2 (2026-10-01): verification checks run after the implementation (`src/tools/validate_reparse_checks.py`, S223, 9 populated row groups)
- **SPiP coverage by tier** (tier 1 = curated clinical ENST, 2 = row has `MANE_SELECT`, 3 = other). Tiers 1-2 are never `not_applicable_transcript`, but SPiP has **no entry for the row's own NM** on many of them: tier 1 36,776 rows -> 90.6% matched / 9.4% no entry; tier 2 895,997 rows -> 68.4% matched / 20.0% no entry / 11.6% no prediction; tier 3 8,060,275 rows -> 85.8% not applicable (no RefSeq equivalent in the pipeline's resources), 0.2% matched. Old SPiP >= 0.20 values removed: tier 1: 10, tier 2: 531, tier 3: 18,912 rows. The tier-3 blanking is a limit of this pipeline's transcript map (MANE + curated only), not evidence that SPiP has nothing for those transcripts.
- **Allele filter**: recounting the 184,721 SPiP `no_entry_for_row_gene` rows with the allele filter off gives 0 extra matches, so it is not dropping valid SPiP entries.
- **Independent check of accepted SpliceAI rows**: `alias_resolved` 192,789 rows - 0 disagree with the coordinate table; 428 are HGNC-only (symbol absent from the coordinate table). `matched` (symbol == row `SYMBOL`) 5,392,563 rows - 9,544 (0.18%) have a coordinate-table gene different from the row's `Gene`; in the single-entry subset (24 distinct symbol/gene pairs, 8,829 rows, e.g. TBCE, ASNS, MATR3, MUTYH) HGNC agrees with the **row's** gene every time, i.e. the coordinate table picked an overlapping duplicate gene ID, the symbol match was right. Consequence: the coordinate route can be wrong for overlapping/duplicated gene models; it can only turn a correct alias into `ambiguous` (dropped), not attach a wrong gene's score, because both routes must agree.
- **Memory**: reading the 9 columns of one 1,048,576-row group peaks at 4.09 GB RSS; with both parsers 4.17 GB (+0.08 GB per 1M rows; strings are shared, ~9 table rows per distinct SPiP string). Parser memory is not the constraint; the existing single-process load of the whole table is.
- **dbNSFP-derived scores (open, same bug class)**: in S223 only 37,128 rows carry dbNSFP values (REVEL_score, AlphaMissense_score, MetaRNN_score, SIFT_score, ... none multi-valued), but only 5,507 of them (14.8%) sit on the row of the dbNSFP-reported transcript (`Ensembl_transcriptid == Feature`) and 7,529 (20.3%) on the row of its gene (`Ensembl_geneid == Gene`). Not fixed here; needs its own look at how dbNSFP is joined and what Module 3 does with it.
- **Branchpointer / LaBranchoR (open)**: variant-level catalog intersection, not gene-keyed; not checked.
- **Open decision for the user**: for the 9.4% (tier 1) / 20.0% (tier 2) rows where SPiP has no entry for the row's NM, and for tier-3 rows, do we want to keep them empty (current) or show the same gene's other-transcript SPiP result in clearly labelled separate columns?
