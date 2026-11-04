# UTRAnnotator: high-impact grading of 5'UTR uORF changes (PRED-13)

2026-11-04T16:02:52+01:00. Code base `c0e5bfe` plus the uncommitted change described here. Mario approved: grade UTRAnnotator instead of "any value = HIGH", do not integrate the utr.annotation R package, and note the gnomAD v4.1 constraint item for later (`TODO.md`).

## What was wrong (FACT)

- `build_newImpact` set `NEW_IMPACT=HIGH` for any non-empty `5UTR_consequence` (`acmg.py:102-104` at `c0e5bfe`), and `build_priority_tier` gave the same rows +20 points and a direct Tier 2 route (`scoring.py:80-82`, `:223`).
- `F2005_union` (2026-11-03 run, 23,686,438 rows): 2,741 rows (1,707 distinct variants, 359 with a MANE Select row) carry a UTRAnnotator class out of 15,086 distinct 5'UTR variants. All 2,741 went from a non-HIGH VEP `IMPACT` to `NEW_IMPACT=HIGH`; 65 were M2 Tier 2. Classes: uAUG gain 868, uAUG loss 766, uORF frameshift 665, uORF stop loss 251, uORF stop gain 191 rows.
- The `acmg.py`/`scoring.py` branch for utr.annotation columns (`utr_num_kozak_gainedOrLost`, `mrl_gainedOrLost`, …) read columns that only the utr.annotation R package produces; this pipeline does not run it.

## Rule (from the literature, not tuned on our data)

According to PubMed, Whiffin et al. 2020, *Nat Commun* 11:2523 ([doi:10.1038/s41467-019-10717-9](https://doi.org/10.1038/s41467-019-10717-9)), the paper the UTRAnnotator plugin cites, reports:
- uAUG-creating variants that form an overlapping ORF (oORF) show selection "equivalent to missense variants". So MODERATE, not HIGH.
- Methods, "Sub-classifying genes": a uAUG gain is high-impact when it forms an oORF or elongates the CDS (in-frame oORF) from a strong or moderate Kozak start. A stop-removing variant is high-impact when the uORF start has a strong or moderate Kozak and/or translation evidence (sorfs.org), and the variant gives an oORF or CDS elongation.
- Discussion: frameshifts that turn an existing uORF into an oORF are also predicted high-impact. The authors also write that even high-impact variants may only cause partial LoF and need functional or segregation evidence.

New column `UTRAnnotator_status` (`predictors.py`, `parse_utrannotator`), per row:

| Value | Meaning |
| --- | --- |
| `high_impact` | uAUG gain with `type` in-frame/out-of-frame oORF and `KozakStrength` Strong/Moderate; uORF stop loss with `AltStop=False` (no other stop in the UTR, reads into the CDS) and Kozak Strong/Moderate or `Evidence=True`; uORF frameshift `ref_type=uORF` → `alt_type` oORF with the same Kozak/evidence condition (**INFERENCE**: our extension of the stop-removing rule, matching the paper's Discussion and its NF2 case) |
| `other` | any other UTRAnnotator class; uAUG loss and uORF stop gain are never high-impact |
| `unclassified` | the deciding field is missing (e.g. `KozakContext=-` when the Kozak context is unreadable near the cap or next to the CDS, `Evidence=NA`) |
| `no_uorf_effect` | 5'UTR row the plugin evaluated (`Existing_uORFs` set) without a uORF change |
| `not_evaluated` | 5'UTR row with no plugin output |
| `not_applicable` / `not_annotated` | not a 5'UTR row / `5UTR_consequence` column absent |

`NEW_IMPACT`: `high_impact` raises LOW/MODIFIER to MODERATE and never lowers a row that is already MODERATE or HIGH. In `build_priority_tier`, only `high_impact` counts as UTR disruption. Gene constraint (LOEUF) is not part of the column: the paper uses it to weigh genes, and that belongs to interpretation (Module 3). The dead utr.annotation branches are removed.

## Checks (local)

- `tests/python/test_utrannotator_status.py`: 5 passed (18 synthetic annotation cases, the no-effect/not-evaluated states, MODERATE without downgrading a Pangolin-HIGH row, status derived when the column is absent). Full runner: only `test_provenance_stamp` and `test_tag_transcript_priority` fail, the known local parquet-metadata failures (neither imports the parser).
- End to end on `concat_fixture_1` (the annotated VCF rebuilt to TSV with `vcf_parser_pysam.py`, then `filter_variants.py` with the pipeline's arguments, outputs in the session scratchpad): 199 rows, same keys. Only 2 rows differ: ITGB1 `chr10:32932877-A-G` on ENST00000474568 and ENST00000677363, a uAUG loss of a moderate-Kozak uORF: `NEW_IMPACT` HIGH → MODIFIER, M2 tier 2 → 3, score 45 → 25; `UTRAnnotator_status=other`. `SCORE_EXPLANATION` and `TIER_RATIONALE` change on the same 2 rows; no other column changed except the new one (`TRANSCRIPT_PRIORITY_TIER` is added later by `TAG_TRANSCRIPT_PRIORITY`). The 197 other rows are `not_applicable`, with `Existing_uORFs` null as read from the TSV.
- `F2005_union` UTR rows, recomputed in memory: 193 `high_impact` rows (137 variants, 29 with a MANE Select row; uAUG gain 119, frameshift 39, stop loss 35), 2,539 `other`, 9 `unclassified`. `NEW_IMPACT` now HIGH 31 (other evidence), MODERATE 236, LOW 39, MODIFIER 2,435. M2 Tier 2: 65 → 7 (16 to Tier 3, 42 to Tier 4). The 27,650 F2005 5'UTR rows without a UTRAnnotator class give 27,648 `no_uorf_effect` and 2 `not_evaluated`. The 9 `unclassified` rows have `KozakContext=-` (uAUG near the cap, next to the CDS, or `CapDistanceToStart=-1`). Many high-impact rows are common (e.g. gnomAD AF 0.41, 0.95); frequency still decides through the tier caps.
- Not run: the cluster fixture with the new code (pending).

## Plugin limitation found (PRED-15)

UTRAnnotator reports one class per row: `UTRAnnotator.pm:262` loops over `keys %five_prime_flag` and keeps the last flagged class, so with two effects the reported one depends on Perl's hash order. FACT: between `F2005_union` and `F1830_union_avi`, 47 of 1,780 shared (`Locus`, `Feature`) rows (34 variants) report a different class, and 6 change `high_impact` status. The NF2 variant in the paper (uAUG gain with a short uORF plus a frameshift into an oORF) is the kind of case that can be reported as the lower class. Not patched (shared plugin).

## Pending

- Commit, then the cluster fixture from that commit, with `concat_fixture_1`'s flags: `bash run_annotate_vcf.sh -profile standard --input_vcf test_data/raw_vcfs/clinvar_af_fixture.vcf.gz --annotation_cache_dir /mnt/scratch/mruizp/annotation_cache_test --chunk_size 5 --gene_restricted_chunk_size 2 --run_id utr5_fixture_1`. Pass: the table equals `concat_fixture_1` except `UTRAnnotator_status` and the 2 ITGB1 rows (`NEW_IMPACT`, `PRIORITY_TIER`, `VARIANT_PRIORITY_SCORE`, `TIER_RATIONALE`, `SCORE_EXPLANATION`).
- Module 3 still treats any non-empty `5UTR_consequence` as UTR disruption (`clinical_variant_prioritization/src/modules/features.py:116-119`) and expects utr.annotation columns. It should read `UTRAnnotator_status` (Module 3 session). Its README (`README.md:10`) still lists Kozak/PolyA/MRL.
- M2's own Tier 2 shortcut still admits `high_impact` rows directly (`scoring.py`, `tier2_mask`); tiers belong to Module 3 (`DEC-0008`, proposed).
