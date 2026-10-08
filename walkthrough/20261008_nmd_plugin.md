# NMD escape annotation (VEP NMD plugin)

Updated: 2026-10-08T16:24:42+02:00

## Scope and change

Mario approved adding the VEP NMD plugin (2026-10-08) after the review noted that PVS1 (`acmg.py`) is "LoF consequence and pLI ≥ 0.9" with no NMD-escape check, which the ClinGen PVS1 decision tree (Abou Tayoun et al. 2018) requires. This change is annotation only: `NEW_IMPACT`, scores, tiers and ACMG tags are unchanged.

- `modules/local/vep.nf`: `--plugin NMD` (no data file; `NMD.pm` is in `/data_lab_PGP/resources/annotation/VEP/VEP_plugins/`, the `--dir_plugins` already used).
- `resources/all_but_old_gnomad_vep_cols.txt`: `NMD` kept by the parser.
- `src/python/parsing/modules/predictors.py` `parse_nmd`, called in `filter_variants.py` after MaxEntScan: derived `NMD_STATUS`.

## What the plugin does (FACT, read from `NMD.pm`)

For rows with `stop_gained`, `frameshift_variant`, `splice_donor_variant` or `splice_acceptor_variant`, it writes `NMD_escaping_variant` when the variant lies in the last exon, in the last ~50 bp of the penultimate exon, in the first 100 coding bases, or in an intronless transcript (Coban-Akdemir 2018; Lindeboom 2016). Otherwise it writes nothing, including when VEP cannot place the variant in the CDS (canonical splice variants are intronic, so they are never evaluated).

## NMD_STATUS

| Value | Meaning |
| --- | --- |
| `escape_predicted` | plugin flag on an in-scope row |
| `nmd_predicted` | `stop_gained`/`frameshift_variant` with a complete `CDS_position` (no `?`) and no flag |
| `not_evaluated` | in-scope row the plugin could not place in the CDS (canonical splice, incomplete CDS position) |
| `not_annotated` | in-scope row in a table without the `NMD` column (runs before this change) |
| `not_applicable` | other consequences |

## Limitations

- Rules are judged at the variant position, not at the premature stop codon. A frameshift in the penultimate exon whose new stop falls in the last exon is not flagged.
- The first-100-coding-bases rule (reinitiation) is not an NMD-escape rule in the ClinGen PVS1 tree, which handles start-proximal variants separately.
- The penultimate-exon window is about 51–52 bp depending on strand in the plugin code.
- Predicted, not observed; per transcript row.

## Validation

Local only: `bash tests/run_all_tests.sh`, 214 passed (new `tests/python/test_nmd_status.py`; `test_vep_custom_columns.py` checks the flag and selector). Not run on the cluster: the plugin has not been loaded by the VEP 111 container yet.

## Open

Wire `NMD_STATUS` into PVS1 (ClinGen tree: escape → PVS1 downgraded) after the cluster check; that changes tiers and needs Mario's decision. See TODO Phase 5.
