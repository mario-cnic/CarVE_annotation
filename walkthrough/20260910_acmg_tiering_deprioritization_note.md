# 📌 ACMG Criteria & Priority Tier — Deprioritization Note

**Date & Time**: 2026-09-10
**Status**: Not yet deprecated, but no longer the primary classification path — documented here so this isn't lost/rediscovered from scratch later.

---

## What changed

Variant classification and dynamic filtering for real clinical/research use now happens in **`clinical_variant_prioritization`** (`/home/mruizp/data_lab_PGP/pipelines/clinical_variant_prioritization`), the standalone web app that was forked out of this pipeline's shared library on 2026-08-31 (see [walkthrough/20260831_standalone_repo_migration.md](./20260831_standalone_repo_migration.md) and [20260831_standalone_web_app.md](./20260831_standalone_web_app.md)). That app has its own independent, actively-developed copies of the ACMG/gene-curation/engine logic (confirmed during this session's audit — see [20260910_full_technical_methodological_audit.md](./20260910_full_technical_methodological_audit.md), Resolution Log for C1) and is where classification and filtering are actually being done dynamically today.

This means the ACMG-adjacent logic still living in **this** pipeline's shared dependency — `shared/utils/src/modules/acmg.py` (`build_newImpact`, `build_acmg_criteria`) and `modules/scoring.py` (`build_priority_tier`, `VARIANT_PRIORITY_SCORE`) — is **not the thing driving real classification decisions anymore**. It still runs on every `vcf2parsed.sh` invocation and still populates `NEW_IMPACT`, `ACMG_CRITERIA`, `PRIORITY_TIER`, `VARIANT_PRIORITY_SCORE` in `results/<GENE>.parsed.clean.pq` and the HTML reports — it just isn't the tool being used to actually classify variants going forward.

## What this changes about the earlier audit findings

The [2026-09-10 full audit](./20260910_full_technical_methodological_audit.md) and [predictor inventory](./20260910_predictor_inventory_and_audit.md) flagged several real issues in this exact logic:

- **ACMG-1**: `PP3_Strong` assigned from a single tool — ACMG/ClinGen SVI caps this at Supporting.
- **ACMG-2**: two disagreeing BA1/BS1 frequency thresholds.
- **ACMG-3**: `VARIANT_PRIORITY_SCORE` is an uncalibrated, uncited heuristic.
- **ACMG-4**: hardcoded MYH7 PVS1 veto duplicated in two files.
- **ACMG-5**: compound-het trans-phasing misses de novo + inherited pairs.
- **ACMG-6**: gene-disease curation outside version control.
- **PRED-12** (already fixed): SpliceVarDB classification silently inert due to a column-name mismatch.

**None of these findings are being withdrawn or downgraded in severity** — they're still real defects in the code as it exists, and PRED-12 in particular was worth fixing regardless (it's cheap, and the column still feeds the report even if it's not the primary classification tool). But this context changes **prioritization**: pouring more effort into ACMG-1/2/3/4/5/6 here is lower-value than it looked before this conversation, since the pipeline's own tiering output is heading toward being a secondary/reference signal rather than the actual clinical decision surface. If `clinical_variant_prioritization` is where classification really happens, its own ACMG/scoring modules (which are separate, independently-evolved code — not audited in this session) are the ones that matter more for correctness right now.

## What to do with this

- **Don't invest further effort fixing ACMG-1 through ACMG-6 here** unless something changes (e.g., a decision to keep using this pipeline's tiering as a cross-check, or a decision to formally deprecate it and these become moot anyway).
- **Do keep flagging new findings in this logic if found** — "not the primary path" isn't "dead code," and per the user, it is explicitly "not yet deprecated."
- **If/when this logic is formally deprecated**, that's the point to: update `README.md`'s feature list (item 6, "Automated ACMG/ClinGen priority tiering," currently describes what's now `clinical_variant_prioritization`'s job, not this pipeline's), update `resources/column_dictionary.md`'s `PRIORITY_TIER`/`ACMG_CRITERIA`/`VARIANT_PRIORITY_SCORE` entries with a deprecation notice, and reconsider whether `build_newImpact`/`build_acmg_criteria`/`build_priority_tier` should still run at all in `filter_variants.py`'s `main()` (right now they're not gated behind any flag — they always run).
- Tracked in [BUG_TRACKER.md](../BUG_TRACKER.md) with a note next to the ACMG-* rows pointing back here.
