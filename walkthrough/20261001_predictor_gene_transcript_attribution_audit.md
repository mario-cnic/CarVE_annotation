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

### Addendum 3 (2026-10-01): dbNSFP per-transcript selection + labelled same-gene SPiP columns (`f805585`)
**Correction to Addendum 2.** The dbNSFP finding there was wrong: I looked for `;` (the table uses `&`) and compared a whole `&`-joined list with a single ID. Measured properly on the 37,128 S223 rows that carry dbNSFP values: the row's gene IS in dbNSFP's gene list on 37,089 (99.9%); 29,573 rows (80%) have multi-valued (`&`) fields; the row's own transcript is in dbNSFP's transcript list on 27,872 (75.1%) and absent on 9,256 (24.9%, almost all `missense_variant` rows of transcripts dbNSFP did not score); the row's protein is never in the protein list when its transcript is not. So the real problem is **per-transcript lists attached whole to every transcript row** (up to 25 values per field), not wrong-gene scores.

**Fix.** 61 dbNSFP columns are transcript-aligned (item count == `Ensembl_transcriptid` item count on >= 99.9% of multi-valued rows; list derived from data by `src/tools/derive_dbnsfp_alignment.py` -> `resources/dbnsfp_transcript_aligned_columns.txt`, e.g. SIFT/Polyphen/AlphaMissense/EVE/ESM1b/VEST4/FATHMM/... scores and the id lists themselves). `parse_dbnsfp_by_row_transcript` (`predictors.py`) keeps for each row the element at the position of its own `Feature`; a transcript not in dbNSFP's list leaves them empty (`dbNSFP_match` = `no_entry_for_row_transcript`); `dbNSFP_transcripts_all` keeps the original list; single-valued columns are untouched. The old `parse_dbnsfp` stays in the call chain but is a no-op in whole-VCF mode (needs a single-transcript table). New CLI arg `--dbnsfp-aligned-columns` is required when dbNSFP columns are present (no silent fallback).
Validation on S223 row group 0: 27,872 matched / 9,256 no entry / rest no dbNSFP; no aligned column is multi-valued afterwards (AlphaMissense_score: 29,573 -> 0); independent check on matched rows: selected `Ensembl_geneid` equals the row's `Gene` on 27,827 of 27,872 (99.84%; the other 45 are transcripts that dbNSFP lists under another gene), selected `Ensembl_proteinid` equals the row's `ENSP` on 27,871 of 27,872. 5 new tests. Both parser stages re-run on the MYBPC3 test VCF: rc 0, and the identity/aligned-column arguments are refused when omitted.

**Not changed (listed so they are not assumed fixed):**
- `MetaRNN_score`, `MetaLR_score`, `REVEL_score` (the three columns the pipeline's `MISSENSE_COLUMNS` consumes) are single-valued in the table, even where the transcript list has up to 25 entries. They are variant-level as delivered; it cannot be told from the table whether dbNSFP/the plugin collapsed a per-transcript field to one value (max, first, ...). Open.
- 26 multi-valued dbNSFP columns whose lists are NOT aligned with the transcript list (Aloft_*, MutPred_*, GTEx_*, clinvar_*, Geuvadis, Interpro_domain, MutationTaster_*): left raw.

**SPiP labelled results (user decision).** Strict `SPiP_*` columns still hold only the row's own-NM entry. New labelled columns describe the GENE: `SPiP_samegene_prediction` (max over the row's gene's SPiP transcripts; symbol matched to the row's gene via the same symbol or HGNC), `SPiP_samegene_transcript` (the `NM_` it comes from), `SPiP_samegene_n_transcripts`. On the MYBPC3 test they fill 148 of 227 rows whose strict SPiP is empty.

### Addendum 4 (2026-10-01): dbNSFP bug found in review, regression checks, updated Module 3 note
- **Bug found and fixed before any run:** in the first version of `parse_dbnsfp_by_row_transcript`, when dbNSFP listed a single transcript that was not the row's own (e.g. COL6A5: list `ENST00000265379`, rows `ENST00000312481`, ...), the single-value branch ran before the position check and the row kept the foreign transcript's values while its status said `no_entry_for_row_transcript` (about 7,555 S223 dbNSFP rows). Now a transcript not in the list leaves every aligned column empty. Re-validated on S223 row group 0: no-entry rows with any aligned column non-null = 0 of 9,256; matched rows with `Ensembl_transcriptid == Feature` = 27,872 of 27,872; no aligned column holds a single value on a multi-transcript list. Two regression tests added (22 tests total).
- **SPiP regression check** after the rewrite of the loop/cache key: row group 16 strict counts identical to the earlier run (matched 100,994 / no entry 14,391 / not applicable 933,191; SpliceAI and Pangolin results identical too).
- **Module 3 note (replaces the earlier one):** column names unchanged. (1) SpliceAI/Pangolin/SPiP: values are only the row's own gene/transcript; empty rows carry a `*_match` reason; variant-level maxima in `*_anygene_*`; SPiP also has labelled gene-level `SPiP_samegene_*` columns (Module 3 does not reference any new column yet). (2) dbNSFP: the 61 transcript-aligned columns used to be `&`-joined lists of up to 25 values on every transcript row and are now one value (or empty, 24.9% of dbNSFP rows have no entry for their transcript). Module 3 references 28 of these columns by name in its own `src/modules/predictors.py` (its copy of the dbNSFP parsing) and `src/report_generator.py` (`SIFT_pred`, `Polyphen2_HVAR_pred`, `aapos`), so its report will show single values and blanks instead of lists. Module 3's app reads the `.pq` as-is and does not call its own `parse_*`.

### Addendum 5 (2026-10-01, evening): SPiP gaps, transcript-mapping audit, single-valued dbNSFP scores (`05a470a`)

**1. Why SPiP had no entry for the row's own transcript.** SPiP scores the RefSeq transcripts of its own bundled database (`RefFiles/dataRefSeqhg38.RData`: 65,484 transcripts, an older RefSeq snapshot; only 9 accessions in the `NM_00137x` range). Current MANE accessions (e.g. `NM_001371727.1`) are often absent while an older accession with the same exon structure is present (`NM_000813`). Among S223 row group 16 MANE rows that carry a SPiP value: MANE accession present 99,144; absent but SPiP has other transcripts of the same symbol 9,222; absent and no entry for the symbol 4,979. Accession matching alone cannot bridge this, exon structure can.
**Fix.** `src/python/build_enst_spip_map.py` maps each GENCODE 45 transcript (252,930) to SPiP-database transcripts on the same chromosome/strand with an identical intron chain (single-exon: identical bounds) -> `resources/enst_to_spip_nm.grch38.gencode45.tsv` (50,280 pairs: 4,770 identical exons, 42,880 same intron chain, 2,630 single exon; 50,139 Ensembl transcripts). `parse_spip` uses, in order: the row's own accession (`matched`), then the structural counterpart (`matched_by_exon_structure`; same splice junctions, UTR extent may differ), else `no_entry_for_row_gene`. Allele filter off recovers 0 rows (not a cause).
**Result on S223 (8,993,048 rows):** strict SPiP present 659,542 (7.3%) -> 1,875,156 (20.9%). Tier 1 (curated ENST) 89.8% (was 90.6%: see 2), tier 2 (MANE) 68.8%, tier 3 15.2% (was 0.2%). What remains empty is a transcript with no counterpart in SPiP's database (tier 2: 19.7% have no entry; tier 3: 68.4% not applicable, 2.4% no entry). Closing that needs a newer SPiP RefSeq database (`RefFiles/getRefSeqDatabase.r`) - not done. Old SPiP >= 0.20 values removed: tier 1: 10, tier 2: 522, tier 3: 16,348.
**Limit:** a same-intron-chain counterpart can differ in start codon/UTR; the status `matched_by_exon_structure` lets you filter these out.

**2. Audit of the curated `resources/gene_transcript_mapping.txt` (217 rows; `resources/gene_transcript_mapping.audit.2026-10-01.tsv`).** All 217 ENSTs exist in GENCODE 45; all NM accessions except OBSCN's exist in SPiP's database. Findings:
- **28 genes where the curated NM and ENST are different transcripts** (the ENST's own MANE/RefSeq counterpart is another accession): LAMP2(NM_013995->NM_002294), TNNT2(NM_001001430->NM_001276345), PKP2(NM_004572->NM_001005242), TBX5(NM_000192->NM_181486), FHL1(NM_001159702->NM_001159699), JUP(NM_021991->NM_002230), ABCC9(NM_005691->NM_020297), AKT1(NM_005163->NM_001382430), ALMS1(NM_015120->NM_001378454), BSCL2(NM_032667->NM_001122955), COQ2(NM_015697->NM_001358921), COX15(NM_004376->NM_078470), CRYAB(NM_001885->NM_001289808), DTNA(NM_001390->NM_001386795), DYSF(NM_003494->NM_001130987), FHL2(NM_201557->NM_001318895), GATA4(NM_002052->NM_001308093), MRAS(NM_012219->NM_001085049), NDUFB11(NM_019056->NM_001135998), NNT(NM_012343->NM_182977), NONO(NM_001145408->NM_007363), NRAP(NM_001261463->NM_198060), PDHA1(NM_001173454->NM_000284), PPP1R13L(NM_001142502->NM_006663), SLC25A3(NM_005888->NM_002635), SYNE2(NM_015180->NM_182914), TJP1(NM_003257->NM_001330239), VARS2(NM_001167734->NM_020442). The previous code paired the curated NM with the ENST to find SPiP entries, which would have scored the wrong isoform for these; that pairing is no longer used. **Not edited** - which of the two the clinicians intend is a curation decision.
- Version-only differences versus MANE v1.4 (e.g. `NM_005159.4` vs `.5`): 155 rows, harmless (versions are stripped everywhere).
- Symbol problems: `TAZ` (VEP: TAFAZZIN) listed twice with two ENSTs (`ENST00000369776` not MANE, `ENST00000601016` MANE) both with `NM_000116.4`; `ASNA1` (VEP: GET3); `NKX25` (typo) duplicates `NKX2-5`. Tier tagging looked transcripts up by gene symbol, so these genes could never reach tier 1. **Fixed:** `tag_transcript_priority.py` now defines tier 1 as "Feature is a curated ENST" (216 distinct transcripts); MYBPC3 test unchanged (tiers 1/2/3 = 17/9/212).
- `OBSCN` (`ENST00000680850`, `NM_001386125.1`): no SPiP transcript with that structure.

**3. The three "single-valued" scores - correction.** In the raw dbNSFP, `REVEL_score` and `MetaRNN_score` are per-transcript lists (checked with tabix on chr8:17646137 T>C: 7 transcripts, REVEL `0.661;...;.`); `MetaLR_score` is a single value. The table showed single values for all three because our own `parse_missense` collapses `MetaRNN_score`, `MetaLR_score`, `REVEL_score` to the maximum over transcripts. My earlier classification was made on that already-collapsed table. `src/tools/derive_dbnsfp_dependency.py` now classifies all 456 dbNSFP columns from the readme plus raw dbNSFP lines (2,511 raw lines, no column between 50% and 100% aligned): **67 transcript-aligned columns** (now including REVEL_score, MetaRNN_score, APPRIS, TSL, LIST-S2_*), 14 multi-valued but not positionally aligned (Aloft_*, MutPred_*, MutationTaster_*; left raw), 372 single-valued (variant-level, incl. MetaLR_score, CADD_phred), `Interpro_domain` excluded (domain lists use the same separator as the transcript list: broke alignment on the MYBPC3 test). `resources/dbnsfp_column_dependency.tsv` has the full table. Row selection happens before `parse_missense`, so REVEL/MetaRNN now are the row's own transcript's value; the old variant-level maximum is kept as `REVEL_score_anytranscript_max` / `MetaRNN_score_anytranscript_max`.

**4. What depends on what (the pipeline's predictors).**
| Level | Predictors / fields | How a row gets its value |
|---|---|---|
| Gene | SpliceAI (one entry per gene), Pangolin (one entry per gene), gene-level VEP fields (LoFtool, pLI, ...) | entry of the row's gene |
| Transcript | SPiP (RefSeq transcript), dbNSFP transcript-aligned columns (SIFT, PolyPhen, REVEL, MetaRNN, AlphaMissense, EVE, ESM1b, VEST4, FATHMM, ..., 67), VEP CSQ fields (HGVS, consequence, ...) | entry of the row's transcript (accession, else exon structure for SPiP) |
| Variant | CADD, conservation, gnomAD, MetaLR, Branchpointer/LaBranchoR (genome-wide catalog intersection), dbNSFP single-valued columns | same value on every row |

### Addendum 6 (2026-10-01, late): validation of the exon-structure SPiP route (`src/tools/validate_enst_spip_map.py`) - supersedes the SPiP numbers of Addendum 5
- **Agreement with known pairs.** Of 19,404 MANE (ENST, NM) pairs, 18,208 have the NM in SPiP's database and can be tested: the structural map gave the same NM for 17,024 (93.5%), a **different NM only: 0**, no counterpart for 1,184 (SPiP's older version of that transcript has another structure: a recall limit, not a mis-attribution).
- **Cross-gene pairs.** Of the first 50,280 structural pairs, 296 (0.6%) linked an Ensembl transcript to a SPiP transcript of a DIFFERENT gene (readthrough/bicistronic loci, e.g. PYURF/PIGY, LRRC51/LRTOMT, GDF1/CERS1, ABCF2-H2BK1/ABCF2) and 1,530 involved a SPiP symbol HGNC does not know. These are the bug class this work removes, so `build_enst_spip_map.py` now keeps only pairs verified as the same gene (same symbol, or resolved through HGNC): **48,454 pairs** (46,996 same symbol, 1,458 via HGNC); the manifest records the dropped counts. Re-validated on the rebuilt map: 0 cross-gene pairs, 0 "different NM only", 17,012 of the testable MANE pairs agree.
- **Reproducibility.** The map build now exports the SPiP database from `dataRefSeqhg38.RData` itself (`--spip-db-tsv/--spip-rdata/--rscript`); input hashes (RData sha256 `86de4151...19c2f`, GTF, MANE GTF, HGNC) are in `resources/enst_to_spip_nm.grch38.gencode45.tsv.manifest.json`. Rebuild whenever SPiP's database or the HGNC file changes.
- **Final S223 numbers (8,993,048 populated rows, rebuilt map):** strict SPiP present 1,865,109 (20.7%; 7.3% before), of which 1,206,082 by exon structure (13.4% of rows). Tier 1 89.8% (89.2% accession + 0.6% structure), tier 2 68.8%, tier 3 15.1%. Old SPiP >= 0.20 values removed: tier 1: 10, tier 2: 523, tier 3: 16,374.
- **Tier-1 rule check on S223:** symbol-keyed rule 36,792 rows, ENST-keyed rule 36,812; the 20 gained rows are exactly TAFAZZIN (18) and GET3 (2). Behaviour change to know about: both TAZ transcripts in the curated table are now tier 1, including the non-MANE `ENST00000369776` (previously, the later duplicate won in the dict and, because of the symbol mismatch, neither matched).
- **The 28 curated genes** (NM and ENST are different transcripts; list in Addendum 5) also affect tier-1 rows beyond SPiP: the tier-1 row and its HGVSc follow the curated ENST, not the NM a lab may report. Until the intended transcript is chosen per gene, tier-1 rows of those genes should not be used clinically.


### Addendum 7 (2026-10-01, late): what `resources/gene_transcript_mapping.txt` is used for (217 rows = 215 genes, 216 transcripts) and a correction
- **Whole-VCF Nextflow** (`annotate_vcf.nf`): only the ENST column, only in `TAG_TRANSCRIPT_PRIORITY` (tier 1 = curated ENST); all transcript rows stay in the table. The NM column is not used. The run manifest hashes the file. SPiP no longer uses it.
- **Legacy per-gene paths** (`main.nf`, `main.sh`): ENST becomes `--gene_set <ENST> --filter_by Feature` (the table keeps ONLY that transcript's rows); NM drives `src/hpc/gene_coords.sh` -> `shared/utils/src/gene_coords.py` (Ensembl REST) to build the gene's region in the BED; `query_new_transcripts.py --auto-append` and two `src/tools` utilities read/write it.
- **Module 3** does not reference the file or `TRANSCRIPT_PRIORITY_TIER`.
- **Correction:** earlier I said that for the 28 genes with NM != ENST the tier-1 row "and its HGVSc" follow the ENST. In whole-VCF mode only the tier flag follows it (every transcript's HGVSc is in the table); in the legacy paths the ENST decides which single transcript is kept and the NM the gene region.
- **Open decision (user):** 28 genes where the curated NM and ENST are different transcripts (ENST = MANE Select in all 28; the NM is not a MANE transcript in 25, MANE Plus Clinical in FHL1, DYSF, SLC25A3). Options: A keep ENST and fix the NM column (recommended default), B keep NM and change the ENST (candidates can be computed by exon structure), C per gene. Plus: TAZ listed twice, `NKX25` typo duplicate, TAZ/ASNA1 renamed (TAFAZZIN/GET3). The file has not been edited.

### Addendum 8 (2026-10-01, late): curated transcript table corrected (user decision: keep the ENST, fix the NM, remove duplicates)
`resources/gene_transcript_mapping.txt` (217 -> 215 rows; one row per gene, one transcript per gene):
- **28 NM values corrected** to the RefSeq accession (with the MANE v1.4 version) of the row's own ENST: LAMP2, TNNT2, PKP2, TBX5, FHL1, JUP, ABCC9, AKT1, ALMS1, BSCL2, COQ2, COX15, CRYAB, DTNA, DYSF, FHL2, GATA4, MRAS, NDUFB11, NNT, NONO, NRAP, PDHA1, PPP1R13L, SLC25A3, SYNE2, TJP1, VARS2. The old and new value of each is in `git diff`. The previous NM of FHL1, DYSF and SLC25A3 was the gene's MANE Plus Clinical transcript; if that is wanted as a second clinical transcript it can be added as its own row.
- **Duplicates removed:** `TAZ,NM_000116.4,ENST00000369776` (non-MANE; the MANE Select `ENST00000601016` row is kept) and `NKX25,NM_004387.4,ENST00000329198` (typo; the `NKX2-5` row is kept). Effect on tier tagging: `ENST00000369776` is no longer tier 1.
- **Not changed:** NM versions that differ from MANE v1.4 only by version (170 rows; same transcript); gene symbols `TAZ` and `ASNA1` (VEP now reports TAFAZZIN and GET3; tier 1 no longer depends on the symbol, but the legacy `main.sh` still looks genes up by the name in the file).
- **Re-audit** (`build_enst_spip_map.py` audit mode): 215 rows, 215 genes, 215 transcripts; no remaining NM/ENST disagreement and no duplicates. MYBPC3 test tiers unchanged (17/9/212). The ENST->SPiP map is unaffected (identical hash).
- **Cost of keeping the MANE Select ENST:** for AKT1, DTNA and OBSCN SPiP's database has no counterpart of that transcript (it only holds the older isoform), so their strict `SPiP_*` columns stay empty; `SPiP_samegene_*` still gives the gene-level result. (Before the fix, AKT1 and DTNA were "matched" to the older NM, which is a different isoform from the tier-1 ENST.)
- The run manifest hashes the file, so runs before and after this commit are distinguishable.
