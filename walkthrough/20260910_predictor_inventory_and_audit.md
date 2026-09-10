# 🧬 Predictor & Annotation Inventory + Audit

**Date & Time**: 2026-09-10 (continued session)
**Scope**: Every functional/splicing predictor and annotation source actually wired into the pipeline, grounded in the real invocation commands (`src/hpc/annotate_*.sh`) and cross-checked against real production output (a 21-variant test run through the live `vcf_parser_pysam.py` → `filter_variants.py` path), not just the README's marketing description. Continues the [2026-09-10 full audit](./20260910_full_technical_methodological_audit.md) — see [BUG_TRACKER.md](../BUG_TRACKER.md) for live status of every finding below.
**Method**: Read every `annotate_*.sh` wrapper in full to determine what is *actually* run (not what the README claims), then traced each output column through `shared/utils/src/modules/*.py` to see whether it's genuinely consumed downstream or silently dead. Several findings below are new — discovered by cross-referencing the real VEP command against the parsing modules' actual column references and against real captured output, not by re-deriving the earlier fork audits.

---

## 1. Full Inventory — what's actually configured to run

| # | Predictor / Source | Tool & version | Invoked by | What it computes | Merged in as |
| --- | --- | --- | --- | --- | --- |
| 1 | **VEP core** | Ensembl VEP, `--cache_version 111`, `--everything` | `src/hpc/annotate_vep_vars.sh` | Consequence, IMPACT, HGVSc/p, exon/intron position, canonical/MANE, existing variation, native `SIFT`/`PolyPhen`, `--check_existing`, population AFs (1000G/ESP-style, deprecated in VEP but still emitted) | `annVEP.vcf.gz` |
| 2 | **CADD** | v1.7 (VEP plugin) — README/docs say v1.6, see DOC-1 | same, `--plugin CADD` | `CADD_PHRED`, `CADD_RAW` | VEP CSQ |
| 3 | **UTRAnnotator** | VEP plugin | same, `--plugin UTRAnnotator` | `5UTR_annotation`, `5UTR_consequence`, uORF gain/loss, Kozak context | VEP CSQ |
| 4 | **SpliceAI (VEP plugin copy)** | precomputed `spliceai_scores.raw.*.hg38.vcf.gz` | same, `--plugin SpliceAI` | `SpliceAI_pred_DS_AG/AL/DG/DL`, `DP_*` | VEP CSQ — **distinct from** the *local* SpliceAI (-D 10000) run separately (row 9) |
| 5 | **SpliceRegion** | VEP plugin | same, `--plugin SpliceRegion` | Extended splice-region flag | VEP CSQ |
| 6 | **REVEL (dedicated)** | VEP plugin, `new_tabbed_revel_grch38.tsv.gz` | same, `--plugin REVEL` | `REVEL` | VEP CSQ — **see PRED-9: never read downstream** |
| 7 | **Mastermind** | VEP plugin | same, `--plugin Mastermind` | Citation/variant-observation counts | VEP CSQ |
| 8 | **LoFtool / pLI** | VEP plugins | same | Gene-level LoF-intolerance metrics | VEP CSQ |
| 9 | **dbNSFP** | v4.8a, `ALL` fields | same, `--plugin dbNSFP` | ~250 columns: `REVEL_score`, `SIFT_score`/`SIFT4G_score`, `Polyphen2_HDIV/HVAR_score`, `MetaRNN_score`, `MetaLR_score`, `MutationTaster`, `PrimateAI`, `ClinPred`, GERP/phyloP/phastCons conservation, ClinVar fields, and more | VEP CSQ |
| 10 | **AlphaMissense** | VEP plugin, `AlphaMissense_hg38.tsv.gz` | same, `--plugin AlphaMissense` | `am_class`, `am_pathogenicity` | VEP CSQ — confirmed this (not the orphaned `annotate_alphamissense.py`) is the live source |
| 11 | **MaxEntScan** | VEP plugin | same, `--plugin MaxEntScan` | `MaxEntScan_ref`/`_alt`/`_diff` | VEP CSQ — **see PRED-11: never read downstream** |
| 12 | **SpliceVarDB** | VEP `--custom` (VCF), experimentally-validated splicing variants | same, `--custom ...short_name=splicevardb` | `splicevardb_gene/hgvs/method/classification/location/doi` | VEP CSQ — **see PRED-12: classification check is silently broken** |
| 13 | **SpliceVault** | VEP plugin | same, `--plugin SpliceVault` | Empirical GTEx/SRA aberrant-splicing events | VEP CSQ, decoded by `parse_splicevault`/`decode_splicevault_prediction` |
| 14 | **gnomAD v4.1** | VEP `--custom` (VCF), joint exomes+genomes | same, `--custom ...short_name=gnomADv4` | ~90 population-frequency fields (`AF_joint`, `AF_grpmax_joint`, per-population AFs, `AN`/`AC`/`nhomalt`, filter flags) | VEP CSQ |
| 15 | **SPiP v2.1** | standalone R/Perl tool | `src/hpc/annotate_spip_vars.sh` | Empirical + ML splicing probability, decoded mechanism (`Exon_skipping`, `Donor_disruption`, etc.) | Own annotated VCF, merged by `merge_vep_spip.sh` |
| 16 | **Pangolin** | standalone PyTorch model | `src/hpc/annotate_pangolin_vars.sh` | Splice gain/loss deltas, claimed cardiac `heart_lv`/`heart_aa` tissue scores | Own annotated VCF, merged |
| 17 | **Local SpliceAI (`-D 10000`)** | standalone PyTorch model, 20kb window | `src/hpc/annotate_spliceai_vars.sh` | Max delta scores across a wider window than the VEP-plugin copy (row 4) | Own annotated VCF, merged |
| 18 | **Branchpointer** | standalone R tool | `src/hpc/annotate_branchpointer_vars.sh` | Branch-point probability, U2 snRNA duplex energy, disruption flag | Own annotated VCF, merged |
| 19 | **LaBranchoR** | standalone BiLSTM model | same script as Branchpointer | Branch-point probability + acceptor distance | Own annotated VCF, merged |

Rows 1–14 all come from **one single `vep` invocation** (`annotate_vep_vars.sh`) — despite the README describing VEP, CADD, REVEL, AlphaMissense, UTRAnnotator, MaxEntScan, SpliceVarDB, SpliceVault, and gnomAD as if they were separately notable integrations, they are all VEP plugins/customs run in the same process and merged as CSQ subfields of a single `annVEP.vcf.gz`. Rows 15–19 are the five genuinely independent standalone tools, each with their own HPC wrapper script, each merged back in separately by `merge_vep_spip.sh`.

---

## 2. Findings — new in this pass (not covered by the earlier fork audits)

### PRED-9 — CRITICAL-adjacent: two REVEL sources, only the less-authoritative one is ever read
`annotate_vep_vars.sh` runs **both** a dedicated `--plugin REVEL,.../new_tabbed_revel_grch38.tsv.gz` (a specifically curated, versioned REVEL data file — someone deliberately installed this) **and** `--plugin dbNSFP,...,ALL`, which bundles its own copy of REVEL as `REVEL_score`/`REVEL_rankscore`. Confirmed twice: first via the real VEP-columns log line from an actual test run, then directly against the real Stage-2 `filter_variants.py` output dataframe (238-row `MYBPC3` parquet) — both `REVEL` (bare, from the dedicated plugin) and `REVEL_score` (from dbNSFP) are present as real, distinct, surviving columns at the exact point downstream code reads them. `shared/utils/src/modules/predictors.py`'s `MISSENSE_COLUMNS` list only contains `REVEL_score` — the dedicated plugin's `REVEL` column is computed on every single run, survives every intermediate transform, and is never read by anything downstream. Either the dedicated REVEL plugin is pointless dead weight (drop it, save the compute/IO), or dbNSFP's bundled REVEL is the wrong one to be trusting silently (dbNSFP bundles arbitrary third-party versions that may lag the dedicated REVEL release) — either way this is an undocumented, silent choice between two candidate values for the single most load-bearing missense score in the whole pipeline (it directly gates `PP3_Strong`/`BP4` in `acmg.py`).

### PRED-10 — Redundant missense-deleteriousness predictors computed, none of them ever used
`--everything` enables VEP's own native `SIFT`/`PolyPhen` scoring; `--plugin dbNSFP,...,ALL` separately bundles `SIFT_score`, `SIFT4G_score`, `Polyphen2_HDIV_score`, `Polyphen2_HVAR_score` — five distinct SIFT/PolyPhen-family values, computed on every variant, every run. Grepped the entire live codebase (`modules/predictors.py`, `modules/acmg.py`, `modules/scoring.py`, `modules/splicing.py`, `modules/io_qc.py`) for bare `SIFT`/`PolyPhen`: zero references. The only place these names appear at all is the orphaned legacy `vep_parser.py` (already flagged dead in DOC-2) and an old exploratory `.Rmd` figures script. Pure dead compute in the live pipeline.

### PRED-11 — MaxEntScan computed, never consumed
`MaxEntScan_ref`/`_alt`/`_diff` (a well-established, widely-used splice-site-strength delta score) is enabled via `--plugin MaxEntScan` and appears in real output, but a repo-wide grep of `shared/utils/src/` finds zero references to `MaxEntScan` anywhere outside the VEP invocation itself. It never feeds `build_spliceMAX`, never appears in any threshold filter, never appears in any ACMG/scoring logic. Dead predictor output.

### PRED-12 — CRITICAL: SpliceVarDB's experimental evidence is completely inert due to a column-name mismatch
This is the most consequential new finding in this pass. `build_newImpact` (`shared/utils/src/modules/acmg.py:68-71`) is supposed to escalate `NEW_IMPACT` to `HIGH` when SpliceVarDB (direct experimental — minigene/RT-PCR/RNA-seq — splicing-variant validation) reports a classification:
```python
if "classification" in data.columns:
    splicevardb_high_mask = data["classification"].notna() & (...)
else:
    splicevardb_high_mask = pd.Series(False, index=data.index)
```
The VEP custom-annotation that produces this data is configured with `short_name=splicevardb`, so every field it emits is prefixed: the real merged column is `splicevardb_classification`, not `classification`. Confirmed twice: against the raw `vcf_parser_pysam.py` output, and — to close the gap between "exists somewhere upstream" and "exists in the exact dataframe the broken check reads" — directly against the real Stage-2 `filter_variants.py` output (the same dataframe `build_newImpact` operates on, after `select_columns`/`build_cDNA_MANE`/`build_locus`/`fix_dup_columns`/`parse_genotype` have all already run): `splicevardb_classification` is present, bare `classification` is not. **This `if` branch has therefore never fired in production**: `splicevardb_high_mask` is unconditionally all-`False`, and SpliceVarDB — arguably the single most clinically trustworthy signal available to this pipeline, since it's direct experimental validation rather than a computational prediction — makes zero contribution to the impact/ACMG logic. This should be treated as at least as urgent as the C1 import fix, since unlike a crash, this fails silently and produces a plausible-looking report that simply never uses SpliceVarDB evidence. (The 21-variant `MYBPC3` test set used for verification has no SpliceVarDB hits at all — `splicevardb_classification` is `NaN` for all 238 output rows — so this doesn't demonstrate a *changed* result on this particular data; the defect is structural and confirmed by column-name analysis, not by observing a flipped classification in this sample.)

### PRED-13 — UTRAnnotator's string-based evidence check: NOT independently verified, flagged rather than confirmed
Two back-to-back checks in the same function (`build_newImpact`, `acmg.py:88-110`) handle UTR evidence very differently: the first (structured `gained`/`lost`/`true` columns) is targeted; the second, immediately below it, escalates to `HIGH` on **any** non-empty `5UTR_consequence` string with no differentiation by consequence severity. This is a real asymmetry in the code as written. However, whether it's an actual over-calling bug depends on `5UTR_consequence`'s real value distribution — populated narrowly (only for uAUG/uSTOP-class events) vs. broadly (for every 5'UTR-overlapping variant) — and the only data available this session (the 21-variant `MYBPC3` test set) has **zero** 5'UTR variants (`5UTR_consequence` is `NaN` for all 238 output rows), so this cannot be checked empirically right now. Per this document's own verification standard, this is listed as a **flagged hypothesis for follow-up with real UTR-variant data**, not a confirmed finding — do not action it without checking the value distribution on a dataset that actually contains 5'UTR variants first.

---

## 3. Findings already tracked (from the earlier fork audits) — grouped by predictor, for reference

| Predictor | Finding | Tracker ID |
| --- | --- | --- |
| Pangolin | `heart_lv`/`heart_aa` "cardiac tissue-specific" scores are a flat max across *all* tissues in the raw-string parse path — fabricated tissue specificity for the pipeline's headline cardiac feature | PRED-1 |
| dbNSFP (all fields) | Transcript-matched value resolution (`parse_dbnsfp`, `transcript=None`) is an unreachable no-op on every real multi-row run | PRED-2 |
| REVEL/MetaRNN/MetaLR | `parse_missense` takes the worst-case (max) score across *any* overlapping transcript when ambiguous, inflating pathogenicity vs. the clinically relevant transcript | PRED-3 |
| Local SpliceAI (-D 10000) | Custom parsing takes the first comma-separated gene entry with no symbol matching in a 20kb window that can span multiple genes | PRED-4 |
| SPiP v2.1 | Hardcoded positional field indexing (no schema-version check); silent wrong-gene fallback to `spip_entries[0]` | PRED-5 |
| Local SpliceAI (-D 10000) | Missing scores `fillna(0)` — "not scored" indistinguishable from "confidently benign" downstream | PRED-6 |
| SpliceVault | Decoded "direction" label ignores strand, contradicting its own (correct) position math | PRED-7 |
| (intron-offset, cross-predictor) | HGVS fallback misclassifies `c.-45+3`-style 5'UTR intronic notation as acceptor-side | PRED-8 |
| SPiP / Pangolin / VEP / Branchpointer | No SGE job-failure propagation — a crashed run silently looks completed | ORCH-1 |
| SPiP | A failed run can produce a "phantom completed" checkpoint (non-empty but content-empty output) | ORCH-2 |
| ACMG evidence (all predictors) | `PP3_Strong` assigned from a single tool clearing a threshold — ACMG/ClinGen SVI caps this at Supporting | ACMG-1 |
| ACMG evidence (frequency) | Two disagreeing BA1/BS1 threshold policies computed independently | ACMG-2 |

---

## 4. Cross-cutting observations

- **One VEP process, marketed as eight integrations.** Rows 1–14 of the inventory are all one `vep` call. This isn't a bug, but it does mean any failure of that single process (see ORCH-1: no `set -e` in `annotate_vep_vars.sh`) takes out CADD, REVEL, AlphaMissense, UTRAnnotator, MaxEntScan, SpliceVarDB, SpliceVault, and gnomAD frequencies simultaneously — there is no partial-failure granularity within this step, unlike the five genuinely independent standalone predictors.
- **Redundant computation is a real, recurring pattern, not a one-off.** REVEL (PRED-9), SIFT/PolyPhen (PRED-10) — the pipeline pays the compute/IO cost of multiple overlapping tools and then arbitrarily (or accidentally) uses only one, or none. This is worth a deliberate pass: either explicitly document which source is authoritative for each redundant pair and drop the others from the VEP command (faster, smaller VCFs, less ambiguity), or explicitly reconcile them (e.g. a documented "prefer VEP's own REVEL when present" rule) rather than have downstream code silently default to whichever field name it happens to reference.
- **The "unified status contract" claim (`scored`/`not_covered`/`error`) does not hold for every predictor.** SpliceVarDB has no `_status` column and no dedicated parser at all (unlike SpliceVault, which gets a full `parse_splicevault`/`decode_splicevault_prediction` treatment) — its raw VEP custom-annotation fields pass straight through. MaxEntScan and the dedicated REVEL/SIFT/PolyPhen plugins have no status contract either, but that's moot since nothing reads them (PRED-9 to PRED-11).
- **Recommended priority order for fixing this batch**: PRED-12 first (silent loss of the pipeline's most trustworthy evidence source, one-line fix), then PRED-9 (decide and document which REVEL source is authoritative), then decide whether to drop the dead MaxEntScan/SIFT/PolyPhen plugins entirely (PRED-10/PRED-11) or wire them into evidence scoring. PRED-13 needs a data check (run on any input with real 5'UTR variants and inspect `5UTR_consequence.value_counts()`) before it's worth prioritizing at all — it may turn out to be a non-issue.

---

## Resolution Log

### 2026-09-10 (continued) — PRED-12 fixed and verified

**Root cause, per user**: the column was `classification` before an earlier, unrelated column-renaming refactor renamed VEP custom-annotation fields to be `<short_name>_`-prefixed (`classification` → `splicevardb_classification`); `build_newImpact`'s check was never updated to track that rename. The new prefixed name is the more correct/appropriate one going forward (avoids future collision with any other column that might legitimately be called `classification`), so the fix updates the check rather than reverting the column name.

**Fix**: `shared/utils/src/modules/acmg.py:68-70` — both references to `"classification"` changed to `"splicevardb_classification"`.

**Verified**:
- Synthetic 3-row test: a real `"Pathogenic"` value now correctly escalates `NEW_IMPACT` to `HIGH`; `NaN` and empty-string values correctly leave `NEW_IMPACT` unaffected (confirms the fix is both effective and doesn't over-trigger on missing data).
- Full re-run of the real production path (`filter_variants.py` against the same 21-variant `MYBPC3` test data used throughout this session): exit 0, output shape unchanged (238×573), no crash. This test set has no real SpliceVarDB hits, so it doesn't demonstrate a changed classification on this particular data — the synthetic test is what confirms correctness; this run only confirms no regression.
- `pytest tests/python/ -q`: still 34/34. No existing test covered `build_newImpact`'s SpliceVarDB branch at all — logged as its own gap (this bug could recur silently again without a regression test; worth adding one).

**Committed**: `shared/utils@708a5c7`.

### 2026-09-10 (continued) — PRED-11 addressed: MaxEntScan now genuinely parsed and documented

**User decision**: verified the final parquet is not column-filtered in the default flow (`select_columns` only filters when `--columns` is passed, which `vcf2parsed.sh` never does — confirmed by reading `io_qc.py` directly; `MaxEntScan_ref`/`_alt`/`_diff` were already surviving as raw pass-through columns). "Include MaxEntScan" therefore meant genuinely integrating it, not recovering a dropped column. User chose **minimal scope**: add a real status contract and a documented disruption metric, but do not yet let it move `NEW_IMPACT`/`PRIORITY_TIER` — and specified **≥15% relative decrease from ref to alt** as the disruption threshold.

**Implemented**: `shared/utils/src/modules/splicing.py` gained `parse_maxentscan()` and `NEW_MAXENTSCAN_COLUMNS` (`MaxEntScan_pct_decrease`, `MaxEntScan_disrupted`, `MaxEntScan_status`), following the exact same row-wise/status-contract idiom as `parse_branchpointer`. `filter_variants.py` calls it alongside the other predictor parsers and adds the new (plus the previously-undocumented raw `MaxEntScan_ref`/`_alt`/`_diff`) columns to `first_cols` for visibility. `resources/column_dictionary.md` documents all six columns, explicitly noting the "informational only, does not move NEW_IMPACT/PRIORITY_TIER yet" status.

**Verified, including two corrections caught on review before committing**:
- Synthetic 6-case test: 20% decrease → `disrupted=YES`; 5% decrease → `NO`; site *strengthened* (alt > ref) → `NO` (by design — this metric only flags weakening, not cryptic-site creation); missing ref/alt → `not_covered`. Initial implementation tagged `ref<=0` as `error`, which was wrong — a non-positive reference score is a real, meaningful outcome (the reference sequence isn't a splice-site consensus at all, exactly the situation for a variant *creating* a cryptic site), not a predictor failure. Corrected to `not_covered`, matching the column dictionary's own description of this case as "undefined" rather than "failed"; re-verified all 6 cases after the fix.
- **Before trusting the parser, checked the raw intermediate TSV directly** for whether `MaxEntScan_ref`/`_alt`/`_diff` might contain non-numeric values (e.g. a delimited donor:acceptor pair) that `_to_float` could be silently swallowing to `NaN` — a failure mode that would look identical to "no data in range" but would actually mean the parser never produces a `scored` row on any real dataset. Confirmed via direct field extraction (`awk`) that all 238 rows are genuinely empty/`.` for all three raw columns in this test set, at fields 450–452 of the header — not a parsing bug, a genuine property of this small test set (MaxEntScan only scores variants near a canonical splice site, and none of these 21 `MYBPC3` variants happen to be there).
- Full production re-run (`filter_variants.py` against the real 21-variant `MYBPC3` test set, twice — before and after the `error`→`not_covered` correction): 238×576 (3 new columns), no crash, all rows correctly `not_covered`.
- `pytest tests/python/ -q`: still 34/34.
- **Not independently verified**: the direction of `(ref - alt) / ref` (positive = weakened = flagged) against a real splice-disrupting variant. This test dataset has zero non-missing MaxEntScan values, so there is no real case available this session to confirm the sign convention empirically — only the arithmetic was checked (a synthetic test can't distinguish "correct direction" from "consistently wrong direction," since the test's own expected values were written against the same formula). This predictor has not produced a single `scored` row against real data in this session. Treat the direction as implemented-per-standard-convention but not yet empirically confirmed, and check it against real data with actual MaxEntScan coverage before relying on it.

**Committed**: `shared/utils@19925ed`.
