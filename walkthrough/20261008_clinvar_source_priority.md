# ClinVar source priority

Updated: 2026-10-08T10:03:28+02:00

## Decision and implementation

Mario made custom ClinVar the TOP priority and requested the remaining sources as fallbacks. Active Module 2 parser/scoring only; outdated reporting and downstream analyses are outside scope. No existing result or reference was rewritten; sibling repositories remain unchanged.

FACT: `_extract_clinvar_status` previously omitted `ClinVarVCF_CLNSIG` and used the first recognized column present for the whole table. A synthetic custom-only Pathogenic row was returned as nonpathogenic.

New `src/python/parsing/modules/clinvar.py` selects per-row classification in this order:

1. `ClinVarVCF_CLNSIG` (dated custom VCF).
2. `clinvar_clnsig` (dbNSFP).
3. `CLIN_SIG` (VEP cache).
4. Legacy `CLNSIG`, `ClinicalSignificance`, `ClinVar`, `clinvar`, `CLINVAR_DISPLAY`.

Only empty/placeholder or not-provided classifications permit fallback. A custom benign, VUS, conflicting or other nonempty classification remains authoritative. Canonical pathogenic/likely pathogenic and benign/likely benign terms drive existing masks; mixed pathogenic/benign or uncertain/conflicting entries activate neither mask.

`acmg.py` now uses this selection for impact adjustment and existing ACMG tagging. `scoring.py` uses its review stars and exposes source-selected columns. Original annotations are retained; no annotation fields were removed. Existing score weights and ACMG rules were not recalibrated or scientifically revalidated by this source-selection change.

Review status is paired with its classification source: custom `CLNREVSTAT`, dbNSFP `clinvar_review`, no stars for cache `CLIN_SIG`. Missing/unrecognized review remains null rather than borrowing an older star value. Generic legacy aliases may use an explicitly supplied `ClinVar_Stars`. NCBI review levels are mapped from practice guideline (4), expert panel (3), multiple submitters without conflicts (2), criteria-provided single/conflicting submissions (1), and known no-assertion states (0). Comma and VEP ampersand encodings are accepted. Source: [NCBI review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).

Derived fields: `CLINVAR_SELECTED_CLNSIG`, `CLINVAR_SOURCE`, `CLINVAR_SELECTED_REVIEW`, `ClinVar_Stars`. Unknown stars cannot activate the existing expert-panel/guideline score bonus.

## Validation and remaining work

Local synthetic checks: source precedence, row-wise fallback, opposing older calls, custom VUS/conflict, review levels and VEP encoding, missing review without borrowed stars, mixed significance, legacy/missing inputs, impact and score integration.

- Focused ClinVar, population AF, trio, row-order and VEP schema tests: 46 passed.
- Predictor attribution regressions (separate invocation): 26 passed.
- `git diff --check`: passed.

No HPC annotation was run. Verify custom lookup and derived source/review columns on a fresh small cluster run. Module 3 has an independent copy and needs this contract before recomputing interpretation. Existing results retain their previous interpretation until deliberately regenerated. TODO items 2 and 2a track implementation and remaining validation.
