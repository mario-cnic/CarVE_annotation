# NMD judged at the premature stop codon (SF-01)

Updated: 2026-11-24T12:20:00+01:00

## What changed

`NMD_STATUS` is now computed by the parser from the transcript's exon structure, at the **new stop codon**, instead of being read from the VEP NMD plugin flag (which judged the variant position). Two new columns say why:

| Column | Values |
| --- | --- |
| `NMD_STATUS` | unchanged value set: `escape_predicted`, `nmd_predicted`, `not_evaluated`, `not_annotated`, `not_applicable` |
| `NMD_RULE` | escape: `last_exon`, `penultimate_last_50nt`, `intronless`; NMD: `none`, `start_proximal`; not evaluated: `splice_outcome_unknown`, `ptc_position_unknown`, `no_stop_in_new_frame`, `transcript_not_in_model`, `transcript_version_mismatch`, `hgvs_missing`, `incomplete_3prime`, `position_unknown`, `model_mismatch`; `no_transcript_model` with `not_annotated` |
| `NMD_PTC_BASIS` | `ptc` (judged at the new stop), `variant` (`fsTer?` with the variant already in the last exon), empty otherwise |

Rules (`src/python/parsing/modules/predictors.py`, `parse_nmd`; parameters at the top of the NMD block):

- New stop codon: from HGVSp. Frameshift `p.Xaa<p>Yaa fsTer<n>` -> stop at amino acid `p + n - 1` (HGVS counts the first changed residue as 1); nonsense `p.Xaa<p>Ter`; `delins…Ter` / `ins…Ter` at the position of `Ter`. Its first nucleotide is moved from variant-transcript to reference coordinates by the HGVSc length change (alt - ref), then to cDNA with the 5'UTR length from the exon table.
- Escape when the stop's first nucleotide is in the last exon, in the last 50 nt of the penultimate exon (`NMD_PENULTIMATE_WINDOW_NT`, ClinGen PVS1 tree, Abou Tayoun et al. 2018), or the transcript has one exon.
- `start_proximal`: stop within the first 100 coding nucleotides (`NMD_START_PROXIMAL_NT`). Reported only; `NMD_STATUS` stays `nmd_predicted`. The plugin counted this as escape; in the ClinGen tree it is not an NMD-escape branch. Not judged on transcripts whose CDS start is incomplete.
- Checks before judging: transcript in the exon table, same version as in HGVSc, no `cds_end_NF`/`mRNA_end_NF` tag (3' structure unknown), and `cDNA_position - CDS_position` equal to the table's 5'UTR length (catches `cds_start_NF` transcripts where VEP counts the CDS differently).
- Canonical splice rows stay `not_evaluated` (`splice_outcome_unknown`). The VEP `NMD` column is kept as it was, for comparison.

Exon table: `resources/transcript_exons.grch38.gencode45.tsv` (111,134 coding transcripts of GENCODE 45 = VEP 111 gene models), built by `src/python/build_transcript_exon_table.py` from `gencode.v45.annotation.gtf.gz` (sha256 `56a166e4…1ae68`, same file as `enst_to_spip_nm`); manifest beside it. Passed to `TSV_TO_TABLE` as `params.transcript_exon_table`, recorded in the run manifest resources and in the `tsv_to_table` cache key.

## Checks

**Coordinates, against the plugin** (`src/tools/validate_nmd_rules.py`): reproducing the plugin's own rules (variant position, 52-nt window, first 101 coding bases) from the exon table on the F2005_union table (run of 2026-11-03, code 814ab97 + dirty patch) gives the plugin's flag on 1,675 of 1,677 stop_gained/frameshift rows. The 2 others are insertions at edge cases (one on a `cds_start_NF` transcript, one at an exon-intron boundary).

**What moves** on the same table (all 3,327 truncating and canonical splice rows; the 1,650 splice rows stay `not_evaluated`):

| Rows | old -> new | All transcripts | MANE Select |
| --- | --- | ---: | ---: |
| frameshift | escape -> NMD | 109 | 14 |
| frameshift | NMD -> escape | 50 | 9 |
| stop_gained | escape -> NMD | 26 | 6 |
| frameshift + stop_gained | escape -> not evaluated | 174 | 13 |
| frameshift + stop_gained | NMD -> not evaluated | 96 | 3 |
| unchanged | | 1,222 | 404 |

- escape -> NMD: all 135 rows had the variant in the first 101 coding bases, the plugin's start rule. 67 have a start-proximal stop (`NMD_RULE` `start_proximal`, no longer escape); in 68 the new stop of the frameshift lies beyond the first 100 nt, or the transcript's CDS start is incomplete (3 stop_gained), so `none`.
- NMD -> escape: the new stop of the frameshift falls in the last exon (22) or the last 50 nt of the penultimate exon (28).
- -> not evaluated: mostly transcripts with an incomplete 3' end (128 rows) or an incomplete CDS start that VEP numbers differently (76), and rows with no computable stop (`p.?`, stop-lost frameshifts, intronic HGVSc, missing HGVSp: 53 + 9 + 4). On MANE Select rows only the last group occurs (16 rows).

Local tests: `bash tests/run_all_tests.sh` in `annotation_tests`, 462 passed, none skipped (`tests/python/test_nmd_status.py`: synthetic transcripts for every rule, both strands of the table builder; end-to-end table carries the new columns, chunked = whole).

## Effect on Module 3 (INFERENCE from `clinical_variant_prioritization/src/engine2/modules/pvs1.py`, read 2026-11-24)

Module 3 gives PVS1 (+8) only to `nmd_predicted` stop_gained/frameshift rows, and marks frameshifts as "judged at the variant position" (`frameshift_nmd_at_variant_position`). With a table built by this code:

- frameshift rows with `NMD_PTC_BASIS == "ptc"` are judged at the new stop: the approximation flag is no longer needed for them;
- some rows gain PVS1 (escape -> NMD, including start-proximal stops), some lose it (NMD -> escape or not evaluated). Whether start-proximal stops should keep full PVS1 is Module 3's call: `NMD_RULE == "start_proximal"` lets it decide.

## Pending

- Cluster fixture run with this code (command in the task note).
- Module 3: read `NMD_RULE`/`NMD_PTC_BASIS`; drop the variant-position approximation for tables that carry them.
