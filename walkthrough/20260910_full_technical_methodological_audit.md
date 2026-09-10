# 🔬 Full Technical & Methodological Audit — Universal Variant Annotation Pipeline

**Date & Time**: 2026-09-10 (session audit)
**Scope**: Full-stack audit of ingestion → liftover → predictor annotation → merge → parsing → QC → filtering → prioritization/ACMG → reporting, covering both this repo (`annotation_pipeline_new`) and its runtime dependency, the shared library at `/home/mruizp/data_lab_PGP/shared/utils/src/` (`filter_variants.py` and `modules/*`), which is where most of the actual filtering/QC/scoring logic lives.
**Method**: Direct primary-source reading of every stage's code, cross-checked against actual execution (`python3 -c "import filter_variants"`, `pytest`, bash reproduction of shell-quoting behavior), plus targeted checks against the real 164-gene cardiac BED panel for claimed substring-collision bugs. Five parallel deep-dives covered: (A) core VCF-parsing/QC engine, (B) predictor score parsing, (C) ACMG/scoring/pedigree/gene-curation, (D) ingestion/liftover/predictor-merge, (E) downstream filtering/reporting/auditor.
**Relation to prior audits**: This builds on and sharpens [20260904_pipeline_analysis_and_improvements.md](./20260904_pipeline_analysis_and_improvements.md), which flagged the `ImportError` and hardcoded paths at a high level. This audit traces those issues to root cause and finds the claimed fix (TODO.md Phase 2) was never actually wired into production, plus a much larger set of methodological findings the prior audit didn't cover.

> **Update, same day (2026-09-10, later in session): C1/C2/C3 are now fixed.** See "Resolution Log" at the end of this document. Everything else in this document (C4 onward) is still open as of this update.

---

## Executive Summary — the four things that matter most

1. **The pipeline's core per-gene script (`filter_variants.py`) cannot currently be imported at all.** Every fresh run's `vcf2parsed.sh` Stage 2 will fail for every gene. Verified by direct execution, not inference.
2. **Even when it runs, the "Customizable Filtering Module" the README advertises never reaches any report a clinician opens.** All three report generators and the cohort dashboard read the *unfiltered* `results/<GENE>.parsed.clean.pq`, not `filtered/filtered_variants.pq`.
3. **Quality control is computed but never gates anything, anywhere, by default.** QUAL/GQ filtering is hardcoded off at parse time; the one opt-in QC gate downstream is never exposed via `main.sh`'s CLI.
4. **The unit test suite is green (34/34) while production is broken** — tests import a module (`config.py`) that patches around the exact defect that breaks the real invocation path, which never imports it. This is a false-confidence trap, not just a missed bug.

Everything below is organized by severity, with file:line references and a concrete failure scenario for each finding, so each is independently verifiable.

---

## 🔴 CRITICAL — Pipeline is not currently functional as documented

### C1. `filter_variants.py` fails to import — confirmed by direct execution
`filter_variants.py` imports `get_full_gene_curation_dataframe` and `filter_by_custom_gene_list` from `modules/disease_hpo.py` (lines 26–31/62–67 of `/home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py`). Neither function exists in `disease_hpo.py` — only `load_cardiac_disease_gene_curation`, `apply_disease_phenotype_curation`, `apply_hpo_phenotype_weights`, `get_dynamic_priority`, `prioritize_genes`, `validate_gene_categories`, `filter_by_gene_priority` do.

```
$ cd /home/mruizp/data_lab_PGP/shared/utils/src && python3 -c "import filter_variants"
ImportError: cannot import name 'get_full_gene_curation_dataframe' from 'modules.disease_hpo'
```

`filter_variants.py` is invoked directly as a standalone script by `src/hpc/vcf2parsed.sh:106` for **every gene, every run**. This is not a theoretical edge case — it's the main execution path.

`filter_variants.py` also imports `convert_to_parquet, DuckDBVariantEngine` from an `engine.py` (line 44) that does not exist anywhere under `shared/utils/` (searched recursively, including worktrees and conda site-packages) — a second, independent reason the import chain is broken.

**Failure scenario**: any user runs `bash main.sh --raw-dir ...` today. Every gene's `vcf2parsed` SGE job crashes at Stage 2. No `results/<GENE>.parsed.clean.pq` is produced. Nothing downstream (filtering, plots, reports, audit) has anything to work on.

### C2. The claimed fix for C1 (TODO.md Phase 2) is not actually wired into production
`TODO.md` marks as done: *"In-memory fallbacks for deprecated shared utilities (`get_full_gene_curation_dataframe`, `filter_by_custom_gene_list`) without editing external shared files"* (`src/python/config.py`).

Reading `config.py`:
```python
if not hasattr(_dh, "get_full_gene_curation_dataframe"):
    _dh.get_full_gene_curation_dataframe = None
if not hasattr(_dh, "filter_by_custom_gene_list"):
    _dh.filter_by_custom_gene_list = None
```
This monkey-patches the attribute onto the `modules.disease_hpo` module object **if and only if something imports `config.py` before `filter_variants.py`'s own top-level import statement runs**, and even then it sets the missing names to `None`, not to a working implementation — a later `get_full_gene_curation_dataframe(...)` call would raise `TypeError: 'NoneType' object is not callable` rather than `ImportError`.

Critically: **`vcf2parsed.sh` (line 106) invokes `filter_variants.py` directly as its own Python process — it never imports `config.py` first.** Grepped both `filter_variants.py` and `vcf2parsed.sh`: zero references to `config`. The fix exists in a file that is never on the execution path that actually breaks. The TODO checkbox is false.

### C3. The unit test suite is green because it accidentally routes around C1/C2, masking the production break
`tests/python/test_compound_het.py:18-21` (and `test_branchpoint_predictors.py`, `test_multisample_trio.py`) do:
```python
sys.path.insert(0, os.path.join(PIPELINE_ROOT, "src", "python"))
import config                                    # <-- patches modules.disease_hpo in-process
from filter_variants import mark_compound_het, analyze_pedigree_inheritance, build_priority_tier
```
Because the test imports `config` first, in the *same process*, the monkey-patch lands on the shared `modules.disease_hpo` module object before `filter_variants.py`'s import statement runs — so the import succeeds in tests only. `python3 -m pytest tests/python/ -q` → `34 passed`. Production (`vcf2parsed.sh` → `filter_variants.py` as a fresh, standalone process) never takes this path and is broken.

**This is worth calling out as its own category of problem**: the test suite is not exercising the actual invocation the pipeline uses in production. A 100%-passing test suite is currently providing false assurance about the one script every gene's output depends on.

**Reproduced directly, side by side** (both run from `annotation_pipeline_new/`, production interpreter `${PYTHON_DATASCI}` resolved via `config/env.sh`, which is `/home/mruizp/apps/miniforge3/envs/datasci/bin/python3`):
```
$ cd /home/mruizp/data_lab_PGP/shared/utils/src && python3 -c "import filter_variants"
ImportError: cannot import name 'get_full_gene_curation_dataframe' from 'modules.disease_hpo'

$ /home/mruizp/apps/miniforge3/envs/datasci/bin/python3 -c "import filter_variants"   # exact prod interpreter
ImportError: cannot import name 'get_full_gene_curation_dataframe' from 'modules.disease_hpo'   # identical — not an interpreter/PYTHONPATH artifact

$ python3 -c "
import sys; sys.path.insert(0,'src/python'); sys.path.insert(0,'/home/mruizp/data_lab_PGP/shared/utils/src')
import config          # <- the line only the test suite has
import filter_variants
print('OK with config pre-imported')"
OK with config pre-imported
```
Confirms C3 is a reproduced fact, not a plausible mechanism: the exact same broken module imports cleanly the moment `config.py` runs first, and the production script never does that.

### C3b. Root cause, dated: this isn't a bug in a release — it's an uncommitted, mid-refactor rewrite sitting directly on the production dependency path
The entire `modules/` package that every finding in this document lives inside (`io_qc.py`, `disease_hpo.py`, `acmg.py`, `pedigree.py`, `predictors.py`, `scoring.py`, `splicing.py`) is **untracked** in the shared-utils git repository:
```
$ cd /home/mruizp/data_lab_PGP/shared/utils && git status --short
 M src/filter_variants.py
 M src/vcf_parser_pysam.py
 M tests/test_filter_variants.py
?? src/modules/                     # <- entirely untracked
```
`git show HEAD:src/filter_variants.py | grep "modules\."` returns nothing — the last **committed** version of `filter_variants.py` is a 2,400+ line monolith that does not reference a `modules` package at all (`git diff --stat` shows −2,080/+339 lines on this file alone in the uncommitted working tree). File mtimes place the refactor's timeline precisely: `modules/disease_hpo.py`, `io_qc.py`, `splicing.py`, `predictors.py`, `acmg.py` were all written 2026‑09‑01 13:00; `modules/scoring.py` 2026‑09‑03; `modules/pedigree.py` and the rewritten `filter_variants.py` itself 2026‑09‑04 16:04–16:08 (`config.py`'s monkey-patch attempt in *this* repo was committed 5 minutes after that, 16:08, suggesting someone hit C1 immediately and tried to patch around it rather than finishing the rename).

**What this means**: the modular architecture this whole audit examines is not a released, versioned state of the shared library — it is live, in-progress, uncommitted development sitting directly in the path `vcf2parsed.sh` calls in production, with no committed checkpoint since 2026‑09‑01 and no rollback safety net. This is a more serious provenance problem than "there's a bug": a `git stash`, a disk issue, or simply someone else's unrelated `git checkout`/`git clean` in that shared directory would silently destroy over a week of unversioned refactor work, and there is currently no tagged/committed "last known good" state of `filter_variants.py` that matches the modular architecture the rest of this document describes.

**This also resolves an apparent contradiction between C1 and C4/C5 below**: if `filter_variants.py` cannot run, how do `test_data/test_run/results/MYBPC3.parsed.clean.pq` (mtime 2026‑08‑28) and its regenerated HTML reports (mtime 2026‑09‑04 16:49) both exist? Because `results/MYBPC3.parsed.clean.pq` **predates the refactor** (Aug 28, before the Sep 1–4 `modules/` rewrite even started) — it's a stale artifact from when the old monolithic `filter_variants.py` still worked. The Sep 4 16:49 report regeneration didn't need `vcf2parsed`/`filter_variants.py` to work at all: `main.sh`'s checkpointing (`is_valid_file()`, C-series finding below) sees the Aug 28 parquet already exists and is non-empty, skips Stage 2 entirely, and the report generators read that stale file directly. **Any run of this pipeline since 2026‑09‑01 that appears to "work" is either hitting pre-refactor cached outputs, or is not actually exercising the real per-gene parse-and-filter step.**

---

## 🔴 CRITICAL — Filtering & QC are structurally disconnected from every deliverable

### C4. Filtered output is orphaned — nothing downstream reads it
`main.sh` runs, per gene, after `vcf2parsed`:
- `filter_and_summarize.py` → writes `filtered/filtered_variants.pq` (main.sh:552-556)
- `plot_annotation_results.R --input "$final_output_file"` (main.sh:558-562)
- `generate_interactive_report.py --input "$final_output_file"` (main.sh:564-568)
- `generate_clinical_prioritization_report.py --input "$final_output_file"` (main.sh:570-574)

`$final_output_file` is `results/<GENE>.parsed.clean.pq` — the **unfiltered** table. All three reports, plus `generate_cohort_master_dashboard.py` (which globs `results/*.parsed.clean.pq`), read the unfiltered file. `filtered/filtered_variants.pq` and its `filtering_summary_metrics.tsv` are written to disk and then read by nothing.

**Failure scenario**: a user runs `--max-af 0.001 --min-revel 0.75 --consequences missense_variant`. Every clinical report and dashboard is built from the same unfiltered input the flags were supposed to narrow — the CLI filtering flags have zero effect on what a clinician opens. (Note: the report scripts do internally re-derive their own display thresholds and tier badges from `PRIORITY_TIER`/hardcoded cutoffs — see the "duplicated-logic drift" finding below — so a report is not literally an unfiltered dump of every row in every view; the precise claim is narrower and more important: *the user's own `--max-af`/`--min-revel`/etc. flags are silently ignored by every deliverable*, with no warning that they had no effect.) The advertised "Customizable Filtering Module" (README §4) does not reach the deliverable.

### C5. Quality filtering is hardcoded off at the only stage where it has real, wired thresholds
`src/hpc/vcf2parsed.sh:110` unconditionally passes `--skip_quality_filter` to `filter_variants.py`. This disables:
- `filter_variant_quality(data, qual_column="QUAL", variant_quality_score=20)`
- `filter_genotype_quality(data, columns=GENOTYPE_QUAL_COLUMNS, genotype_quality_score=30, ...)`

for **every** run through the documented orchestrator, with no CLI flag to turn it back on.

The one opt-in QC gate that does exist downstream, `filter_and_summarize.py --pass-qc-only` / `--min-qual` (`src/python/filter_and_summarize.py:39-40`), is **never exposed in `main.sh`'s argument parser or help text** (checked the full `case $1 in` block, lines 44-206 — no such flag).

Combined with C4: **running the pipeline exactly as the README's "Basic Run" example shows applies zero QUAL, zero genotype-quality, zero allele-frequency, and zero pathogenicity-score filtering anywhere in the path that reaches a clinical report.** `filtered_variants.pq` and `<GENE>.parsed.clean.pq` are named as if curated; by default, both are byte-for-byte the same variant set as whatever came out of annotation, QC-failing and common-benign variants included, with `filtering_summary_metrics.tsv` silently showing `Total_Input_Variants == Final_Retained_Variants` and no visible warning that nothing happened.

### C6. QC_STATUS is computed but exists only as a column nobody reads
`build_qc_status` (`shared/utils/src/modules/io_qc.py:149-226`) computes its own QC policy — `QUAL<30`, `QD<2.0`, `FS>60.0`, `MQ<40.0` → `LOW_QUAL` — independently of, and stricter than, `filter_variant_quality`'s `QUAL≥20` used elsewhere. Two different, uncited, disagreeing QC thresholds in one codebase. This is the *only* quality signal that reaches the output given C5, and given C4/the missing `--pass-qc-only` wiring, it never filters anything a user sees. This is functionally equivalent to "not flagging QC," even though a `QC_STATUS` column visibly exists in the table — exactly the failure mode the audit was asked to check for.

Also implemented via a row-wise `.apply(axis=1)` — will not scale to WES/WGS-sized inputs.

### C7. Genotype data is never extracted for "proper VCF" input — compound-het/inheritance features are dead in the main flow
`vcf2parsed.sh:95-102` calls `vcf_parser_pysam.py` with only `--add_info --add_vep` — never `--add_gt/--add_gq/--add_dp/--add_ad`. The intermediate TSV therefore never has `GT_*`/`GQ_*` columns. Downstream, `parse_genotype`, `filter_genotype_quality`, `merge_sample_cols`, and all of `modules/pedigree.py`'s compound-het/inheritance logic silently no-op with a warning log ("No GT_ columns found... skipping").

**Failure scenario**: a trio VCF (proband + parents) is submitted, exactly the "variant + genotype" input mode the pipeline claims to support. Genotypes are discarded before annotation begins. No compound-het pairing, no per-sample inheritance classification, no genotype-quality filtering ever happens, despite dedicated code (`mark_compound_het`, `analyze_pedigree_inheritance`, dedicated unit tests) existing for exactly this. (The unit tests pass because they call these functions directly with a hand-built DataFrame that already has `GT_*` columns — they never exercise the real `vcf2parsed.sh` invocation that would need to produce those columns first.)

### C8. `--consequences` is silently corrupted by shell-quoting when passed through `main.sh` — verified by reproduction
`main.sh:550`: `[ -n "$CONSEQUENCES" ] && filter_flags="$filter_flags --consequences \"$CONSEQUENCES\""`. The `\"` sequences store **literal double-quote characters** inside the `filter_flags` string variable. At the call site (`main.sh:556`), `$filter_flags` is expanded **unquoted**, which word-splits on whitespace but does *not* strip the embedded quote characters (those are just data at that point, not shell syntax). Reproduced directly:
```bash
$ CONSEQUENCES="missense_variant,stop_gained"
$ filter_flags="--consequences \"$CONSEQUENCES\""
$ py() { for a in "$@"; do echo "ARG=[$a]"; done; }; py $filter_flags
ARG=[--consequences]
ARG=["missense_variant,stop_gained"]     # <- literal quote characters land in argv
```
`filter_and_summarize.py`'s `consequences.split(',')` + `.strip()` does not remove quote characters, so the resulting match tokens are `"missense_variant` and `stop_gained"` — neither is ever a substring of any real (unquoted) `Consequence` value. **Every variant is filtered out** whenever `--consequences` is passed through `main.sh`, silently producing an empty `filtered_variants.pq`. Currently masked in practice by C4 (nothing reads that file), but this becomes a live silent-data-loss bug the moment C4 is fixed naively (e.g. by just pointing the reports at the filtered file).

### C9. Every gene's filtered output collides on the same run-level filename — an overwrite bug and an SGE data race
`main.sh:324`: `FILTERED_MASTER_DIR="${RUN_BASE_DIR}/filtered"` is a **run-level** directory, not per-gene (unlike `RESULTS_MASTER_DIR/${gene_name}.parsed.clean.pq`, which is correctly namespaced). `filter_and_summarize.py` (before this fix) always wrote to a fixed `<output-dir>/filtered_variants.<ext>` and `<output-dir>/filtering_summary_metrics.tsv`, with no gene component in the filename. `main.sh:552-568` submits one `filter_${gene_name}` SGE job per gene in the raw-dir loop, all writing into the same `$FILTERED_MASTER_DIR` with the same filenames.

**Failure scenario**: any run processing more than one gene (the normal case — `main.sh` loops over every file in `--raw-dir`) submits N independent, often-concurrent SGE jobs that all write `RUNS/<RUN>/filtered/filtered_variants.pq`. Whichever job's `to_parquet()`/`to_csv()` call lands last wins; every other gene's filtered output for that run is silently destroyed, with no error, no warning, and no indication in `filtering_summary_metrics.tsv` (also collided) that anything but one gene was ever filtered. This is on top of, and independent from, the fact that C4 means nothing was reading this file anyway — but it would have made C4's fix actively wrong if not caught first: pointing per-gene reports at `filtered/filtered_variants.pq` would have shown every gene's report the same single leftover gene's data.

*(Fixed alongside C8 below — see Resolution Log.)*

---

## 🟠 HIGH — Scientific-validity issues in ACMG/scoring/pedigree logic (`shared/utils/src/modules/{acmg,scoring,disease_hpo,pedigree}.py`)

- **PP3 mislabeled "Strong."** `build_acmg_criteria` (`acmg.py:147-168`) tags `PP3_Strong` from a single tool clearing a threshold (SpliceAI≥0.5 / SPiP≥0.5 / AlphaMissense≥0.564 / REVEL≥0.75). ACMG/ClinGen SVI guidance caps single-tool computational evidence at **Supporting**. This overstates evidence strength on every variant so tagged.
- **Two disagreeing BA1/BS1 frequency policies in one pipeline.** `acmg.py:151-152` uses flat AF>0.05 (BA1) / 0.01–0.05 (BS1). `scoring.py:125-131`'s `build_priority_tier` instead looks for `FAF_POPMAX_BA1`/`FAF_POPMAX_BS1` columns that are never populated anywhere in these modules, silently falling back to AF≥0.002 / 0.0001–0.002 — an order of magnitude stricter. The same variant can show "not BA1" in `ACMG_CRITERIA` while independently being pinned to Tier 4 by `PRIORITY_TIER`'s own gate — a self-contradictory report.
- **`VARIANT_PRIORITY_SCORE` is an uncalibrated black box presented as a clinical metric.** `build_priority_tier` (`scoring.py:87-234`) sums ~15 additive rules (LoF +35, splicing composite ×30, ClinVar +45/+50/−40, AF −100/−25/+15, QC −25, inheritance/compound-het +10–20, HPO +15…), clipped to [0,100], with no citation and no validation against a labeled truth set (e.g. ClinVar). This drives the tiers shown to reviewing clinicians. This is the single largest "trust the black box" exposure in the whole pipeline.
- **Hardcoded, duplicated, gene-specific PVS1 veto.** `symbol_series != "MYH7"` appears independently in both `acmg.py:134` and `scoring.py:44` — biologically plausible (MYH7 LoF isn't an established HCM mechanism) but uncited, can drift out of sync between the two copies, and doesn't generalize to other genes with known mechanism-specific caveats (e.g. TTN truncating variants only pathogenic in constitutively-expressed/high-PSI exons).
- **Compound-het logic misses de novo + inherited trans pairs.** `_mark_pairwise_trans` (`pedigree.py:80-139`) only recognizes a het variant as paternal-only/maternal-only when the *other* parent is exactly HOMREF — a genuine de novo variant (both parents HOMREF) paired with a truly inherited het variant in the same gene is a valid trans compound-het configuration that this logic never flags.
- **Gene-disease curation lives outside version control with zero provenance capture.** `load_cardiac_disease_gene_curation` (`disease_hpo.py:36-37`) reads hardcoded absolute paths (`/home/mruizp/data_lab_PGP/shared/utils/data/Complete_gene_list_V4.csv`, `full_cardiac_gene_list_V2.csv`), `lru_cache`d process-wide, with no file hash/mtime/row-count logged anywhere in the output. Combined with C1–C3, there is currently no way to trace which curation snapshot produced a given clinical tier — a direct violation of this project's own Rule of Provenance.
- Gene-panel constant lists (`DCM_CORE`, `HCM_CORE`, etc.) are copy-pasted verbatim in two separate functions in `disease_hpo.py` (lines 72-77 and 117-122) — an edit to one silently desyncs from the other.
- `apply_disease_phenotype_curation`'s own default argument string (`"🫀 Cardiomyopathies (General Panel - DCM, HCM, ARVC, RCM)"`, `disease_hpo.py:102`) contains the substring `"dcm"`, so its own default silently narrows to DCM-only genes despite explicitly claiming to be a general multi-phenotype panel. Not triggered on the live path today (`scoring.py` overrides with a different default that doesn't collide), but a live trap for any other caller.

## 🟠 HIGH — Predictor score parsing correctness (`shared/utils/src/modules/splicing.py`, `predictors.py`)

- **Pangolin "heart_lv"/"heart_aa" tissue-specific scores are fabricated in the raw-string parse path.** `splicing.py:251-299` (raw-string branch, 271-288): when parsing the raw Pangolin INFO string, the code takes `max(abs(score))` across **all** `tissue:score` pairs without checking which tissue, then assigns that single flat max to `Pangolin_max_score`, `Pangolin_heart_lv_score`, and `Pangolin_heart_aa_score` identically. The README specifically advertises these as cardiac-tissue-specific evidence (§2); this path fabricates that specificity. (The alternate path, reading precomputed per-tissue columns, is correct — the bug is specific to raw-string parsing.)
- **dbNSFP transcript-matching is an unreachable no-op on real data.** `parse_dbnsfp` is called as `parse_dbnsfp(data, transcript=None, ...)` (`filter_variants.py:339`); with `transcript=None`, `predictors.py:98-111` sets `transcript = data['Feature'].values` (the whole column) and only proceeds if `len(transcript) == 1` — false for any real multi-row dataframe. Every run hits the `else` branch and returns immediately, so all dbNSFP columns (MetaRNN, MetaLR, PolyPhen2, SIFT, MutationTaster) pass through as raw, ambiguous multi-transcript strings, silently undocumented as such.
- **`parse_missense` (REVEL/MetaRNN/MetaLR) takes the worst-case score across *any* overlapping transcript**, including non-canonical/non-MANE ones (`predictors.py:64-95`), inflating apparent pathogenicity relative to the clinically relevant transcript.
- **Local SpliceAI (`-D 10000`) custom parsing picks the first listed gene with no symbol matching.** `parse_spliceai_custom` (`splicing.py:395-425`) takes only `str(val).split(",")[0]` — in a 20kb window that can span multiple genes, whichever gene's prediction is listed first silently wins for every `spliceai_custom_*` column, unlike `parse_spip`, which at least attempts a symbol match.
- **SPiP parsing uses hardcoded positional field indexing with a silent wrong-gene fallback.** `splicing.py:52-136` pulls fields by numeric index into a pipe-delimited string with no schema-version check, and falls back to `spip_entries[0]` (line 74-75, potentially the wrong gene) when no entry's embedded gene tag matches the row's SYMBOL.
- **Missing SpliceAI scores collapse to numeric 0**, indistinguishable from "confidently benign" for any code that sorts/filters on the numeric column without also checking `SpliceAI_status` (`splicing.py:451-475`, `.fillna(0)` at line 469).

## 🟠 HIGH — Orchestration integrity (`main.sh`, `src/hpc/*.sh`)

- **No SGE job-failure propagation across the whole predictor layer.** `annotate_vep_vars.sh`, `annotate_spip_vars.sh`, `annotate_pangolin_vars.sh`, `annotate_branchpointer_vars.sh`, and `merge_vep_spip.sh` all capture `CMD_EXIT_CODE=$?` but never `set -e` and never `exit $CMD_EXIT_CODE` — each ends with an unconditional success `echo`. `qsub -hold_jid` only waits for job **completion**, not success, so a crashed predictor run lets every downstream job (merge → vcf2parsed → filter → reports) proceed on silently incomplete data. Only `annotate_spliceai_vars.sh` (the most recently added predictor) and `variant_converter.sh` get this right — the safe pattern was never backported to the others. This is a regression/inconsistency, not a considered design choice.
- **A failed SPiP run can produce a "phantom completed" checkpoint.** `annotate_spip_vars.sh`'s post-processing (`grep -v "^##SPiP output v2.1$" $OUTPUT > $TMP_FILE`) creates a non-empty `$TMP_FILE` via shell truncation even if `$OUTPUT` never existed — `main.sh`'s `is_valid_file()` (exists + size>0) and `merge_vep_spip.sh`'s `[ -s "$pred_file" ]` both then treat it as done, and the gene "looks finished" for SPiP on every future run unless `--force-spip` is manually passed.
- **`merge_vep_spip.sh` has no `set -e`, and its `CMD_EXIT_CODE=$?` (line 148) captures the exit status of the last `rm -f` in its cleanup loop, not of any actual `bcftools` merge operation** — the logger records a fake success no matter what happened during the merge itself.
- **Regression of a previously-fixed SGE freeze bug.** `gene_coords.sh` sources `~/.bashrc` in a headless SGE job, while sibling scripts explicitly carry the comment "Do not source ~/.bashrc in headless SGE cluster jobs" (a lesson from `walkthrough/2026-08-06_sge_varconv_freeze_fix.md`) — the fix was never applied here.
- **Substring gene-name collisions, empirically confirmed against the real 164-gene cardiac panel.** `main.sh:378` (`grep -q "$gene_name" "$GENES_BED_FILE"`, no `-w`) and `gene_coords.sh` (`grep -q "^$GENE_NAME" ...`, anchored but not word-bounded) can both false-positive. Checked directly against `resources/cardio_genes_loc.bed`:
  ```
  EMD is substring of LEMD2
  TNNI3 is substring of TNNI3K
  ```
  Both pairs are real, distinct genes in this exact panel. Processing `EMD` when only `LEMD2`'s BED line exists yet would make the existence check return a false "already present," silently skipping BED-coordinate generation for `EMD`.
- **State-aware checkpointing never validates content, only file existence + non-zero size** (`is_valid_file()`, `main.sh:335-337`). Combined with the two points above, a partially-failed run is indistinguishable from a fully successful one to the orchestrator, and simply re-running the same command will not retry the broken step.
- **`gene_name` is derived from the input filename by convention** (`cut -d'_' -f1`, uppercased, `main.sh:352`) — fragile for any filename not following `GENE_*.ext`, or for a cohort-level multi-gene VCF, which the pipeline otherwise claims to accept generally.
- **Transcript resolution silently defaults to the literal string `"UNKNOWN"`** (`main.sh:372-375`) when a gene isn't in `gene_transcript_mapping.txt`, and this string is passed as a real positional argument into `vcf2parsed.sh`.
- **`query_new_transcripts.py --auto-append` runs unconditionally (`|| true`) at the start of every single invocation** (`main.sh:340`) and mutates the shared, unlocked `resources/gene_transcript_mapping.txt` in place — concurrent runs race on this file, and its Ensembl API fallback accepts the first-listed transcript when no `is_canonical` flag is present, potentially writing a non-MANE transcript permanently into a resource all future runs depend on.

## 🟠 HIGH — Genome-build / liftover correctness (`src/python/universal_variant_converter.py`)

- **hg19→GRCh38 liftover silently drops all genotype/FORMAT data, asymmetrically with GRCh38 input.** The direct-VCF liftover branch re-emits only the 8 mandatory VCF columns (CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO); `parts[8:]` (FORMAT + all sample GT/GQ/DP columns) is never read or written, and the output header is hardcoded to 8 columns. The GRCh38-input passthrough path (no liftover needed) instead copies the file wholesale, genotypes intact. **The exact same cohort submitted as hg19 vs GRCh38 VCF gets different feature availability** — hg19 users silently lose genotypes, and therefore compound-het/inheritance/genotype-QC, with zero warning.
- **Liftover never handles strand flips or re-normalizes alleles.** Both liftover call sites use only `lift_res[0][0]` (chrom) and `lift_res[0][1]` (pos) from `pyliftover`; the returned strand flag (`lift_res[0][2]`) is never inspected, so REF/ALT are never reverse-complemented across a strand-flip region, and there is no REF-vs-actual-reference-base validation anywhere to catch the resulting mismatch.

## 🟡 MEDIUM — Structural variant handling is unfinished, undisclosed scaffolding

`main.sh` threads a `large_sv_vcf` path into `merge_vep_spip.sh`, which has real logic to `bcftools concat -a` it back into the final annotated VCF (`merge_vep_spip.sh:129-139`) — but **grep-confirmed zero producers of this file anywhere in the codebase** (`universal_variant_converter.py`, `variant_converter.sh`, and every `annotate_*.sh`). Either large/structural variants are never diverted at all and instead flow through VEP/SPiP/Pangolin/SpliceAI (tools built and validated for point mutations and small indels, with undefined behavior on large deletions/duplications/insertions), or this is scaffolding for a feature that was never finished. The README makes no mention of SV handling either way — there is no user-facing disclosure that this is absent.

## 🟡 MEDIUM — Downstream filtering semantics are correct-ish but undocumented, and duplicated across files

- **NaN-passthrough on every score filter is undocumented.** `filter_and_summarize.py`'s REVEL/AlphaMissense/SPiP/CADD filters (`mask = (val>=threshold) | val.isna()`, lines 113/127/141/155) always **pass** a variant with a missing score for that predictor. This is defensible (a `--min-revel` filter shouldn't wrongly drop a splice variant that REVEL never scores) but is invisible to a user, who is likely to assume `--min-revel 0.75` restricts the output to "high-confidence missense variants" — it does not; it only removes missense variants that were scored and scored low. A variant with **every** predictor score NaN (zero evidence at all) passes every filter and is indistinguishable in the output from one that passed on strong evidence.
- **Two disagreeing AF-column resolution orders.** `filter_and_summarize.py:91` checks generic `AF` **before** the pipeline's actual primary population-frequency field `gnomADv4_AF_grpmax_joint` (which isn't even in its candidate list), while `generate_clinical_prioritization_report.py:74-99` correctly prioritizes `gnomADv4_AF_grpmax_joint` first and visibly warns on fallback. Two different implementations of "what is the population frequency" in one pipeline (currently only affects the orphaned file from C4, but would become live the moment that's fixed).
- **`audit_run_results.py` checks row-count completeness for every intermediate predictor VCF but not the final deliverable parquet** — the one place in the whole audit chain a silent-attrition check is missing is the output that actually reaches a clinician.
- **`generate_clinical_prioritization_report.py` re-derives its own hardcoded evidence-count thresholds** (splice≥0.20/0.50, AlphaMissense≥0.564, REVEL≥0.75) independently of the upstream `PRIORITY_TIER` computation in `scoring.py` — if the latter's thresholds are ever tuned, the report's tier badge and its own "N variants pathogenic by consensus" counter can silently disagree.
- **Absence-from-gnomAD is conflated with AF=0.0** (`filter_and_summarize.py`'s `.fillna(0.0)` before the max-af comparison) — "not observed in a reference database" (which can reflect a coverage gap) is not the same as "confirmed absent," even though the practical filtering direction here happens to be the permissive/safe one.

## 🟢 LOW — Documentation, provenance, and code hygiene

- **CADD version mismatch**: `annotate_vep_vars.sh` invokes CADD **v1.7** in its VEP plugin path; README and the prior (2026-09-04) audit both state "CADD v1.6" — a direct Rule-of-Provenance violation, and an easy fix.
- **Dead/orphaned duplicate implementations** create maintenance traps: `vcf_parser_pd.py` (zero references anywhere, coexists with the actually-used `vcf_parser_pysam.py`); `annotate_alphamissense.py` (orphaned, unreachable from the live `filter_variants.py` import chain); top-level `parsing_dbnsfp.py` duplicates `modules/predictors.py::parse_dbnsfp` under a near-identical name.
- **Misleading log message**: `merge_vep_spip.sh:130-132` always prints "Concatenating 1 large structural variants..." regardless of true count (`head -n 1 | wc -l` can only ever yield 0 or 1) — moot today given the SV section above, but a good example of the general logging-accuracy gap.
- **Swapped `include_ad`/`include_dp` argument order** at one call site in `vcf_parser_pysam.py`'s VEP-annotated code path — dormant today, live landmine if per-sample AD/DP extraction is ever enabled there.
- **Inconsistent binary-path fallback robustness** across `annotate_*.sh` wrappers — some have 2–4-tier local/HPC fallback chains, `annotate_spip_vars.sh` has none.

---

## What actually happens for each of the three claimed input modes

The pipeline is framed around three input modes — variant-only, variant+existing-annotation, and full VCF (variant+genotype+annotation). Findings above were scattered across stages; collected here against that framing:

| Input mode | What the code actually does | What breaks |
| --- | --- | --- |
| **Variant-only** (bare CHROM/POS/REF/ALT or cDNA notation, no genotypes, no prior annotation) | `universal_variant_converter.py` builds a minimal VCF; goes through the full VEP/SPiP/Pangolin/SpliceAI/Branchpointer stack. This is the mode the pipeline is best-tested for. | Everything upstream of C1 still has to work (it currently doesn't). No genotype-dependent features apply here by definition, so C7 doesn't cost anything extra in this mode. |
| **Variant + existing annotation** (e.g. an `.xlsx`/`.pq` that already carries columns like a prior `REVEL`, `Consequence`, or `AF`) | `write_vcf` (`universal_variant_converter.py:266-279`) generically preserves every extra input column as a `String`-typed INFO field (`Description="Input metadata column {col}"`) — so pre-existing annotation values *are* carried through the VCF conversion step, not silently dropped there. But the pipeline then re-annotates from scratch (VEP/predictors), and `vcf_parser_pysam.py` + `fix_dup_columns` (`io_qc.py:45-70`) are what reconcile a name collision between the user's carried-over column and the pipeline's freshly computed one of the same name. `fix_dup_columns` only knows how to merge a plain column with its pandas-style `.1`/`.2` duplicate via `fillna` — **it has no explicit "always trust the fresh annotation" or "always trust the user's pre-existing annotation" policy; the winner is whichever name pandas/the parser happened to assign the non-suffixed slot to first.** For a specific, already-documented instance of this: `filter_and_summarize.py:91` checks a plain `AF` column *before* the pipeline's own `gnomADv4_AF_grpmax_joint` — exactly the column name a pre-annotated input is most likely to already carry, and exactly the collision this mode creates. | A user who supplies pre-annotated data expecting the pipeline to either respect or clearly override their existing scores gets neither a documented policy nor a warning when a name collides — the resolution is an implementation accident, not a decision. |
| **Full VCF** (variant + genotype + annotation, e.g. a trio) | Converter passes the VCF through largely as-is (for GRCh38 input) or drops FORMAT/genotype columns entirely during hg19 liftover (finding under "Genome-build / liftover correctness" above). Even when genotypes survive conversion, `vcf2parsed.sh` never requests `--add_gt/--add_gq/--add_dp/--add_ad` from the parser (C7) — so genotype data is discarded before annotation regardless of whether the VCF still has it. | This is the mode with the most silent capability loss: compound-het, per-sample inheritance, and genotype-quality filtering are all dead code in the actual orchestrated flow, despite dedicated implementations and passing unit tests (which test the functions directly, not the real parse step that would need to feed them). |

---

## Prioritized Remediation Roadmap

0. **Commit the `modules/` refactor in `shared/utils` right now, broken or not**, before anything else touches that working tree (`git add -A && git commit`, or at minimum a WIP tag/branch). There is currently zero version-controlled checkpoint of this architecture (C3b) — fixing C1 without first committing the current state risks losing the ability to bisect what changed, and leaves the shared library one accidental `git clean`/`checkout` away from reverting a week of unversioned work.
1. **Unbreak the pipeline.** Fix `disease_hpo.py` to actually export `get_full_gene_curation_dataframe`/`filter_by_custom_gene_list` (or update `filter_variants.py`'s import to the current real function names), and either implement or remove the `engine.py` import. Verify by running `python3 -c "import filter_variants"` **without** first importing `config.py`, and add that exact invocation as a CI/test-suite check so tests can no longer pass while production is broken (this addresses C1–C3 together).
2. **Wire the filtering module into the actual deliverables**, or remove the "Customizable Filtering Module" claim from the README until it is. Either point the R/report/dashboard scripts at `filtered/filtered_variants.pq`, or fold `filter_and_summarize.py`'s logic into the report scripts directly (C4).
3. **Decide on one QC policy and make it opt-out, not opt-in-and-unreachable.** Reconcile `build_qc_status`'s thresholds with `filter_variant_quality`/`filter_genotype_quality`'s thresholds into a single documented policy; expose `--pass-qc-only`/`--min-qual` through `main.sh`; stop unconditionally passing `--skip_quality_filter` (C5–C6).
4. **Fix the `--consequences` shell-quoting bug** (drop the embedded `\"..\"`; pass the value quoted properly through to `qsub -b y`) (C8).
5. **Decide whether genotype-aware analysis is in scope for the main flow.** If yes, add `--add_gt --add_gq --add_dp --add_ad` to the `vcf_parser_pysam.py` call in `vcf2parsed.sh` and re-verify compound-het/pedigree tests against the real pipeline invocation, not a hand-built DataFrame (C7).
6. **Add `set -e`/`exit $CMD_EXIT_CODE` to every `annotate_*.sh` and `merge_vep_spip.sh`**, matching the pattern already used correctly in `annotate_spliceai_vars.sh` and `variant_converter.sh`, so SGE job failure actually blocks downstream jobs.
7. **Add `-w` to the gene-name `grep` checks** in `main.sh` and `gene_coords.sh` (confirmed-real collisions: EMD/LEMD2, TNNI3/TNNI3K in the current panel).
8. **Resolve the ACMG/scoring policy conflicts**: pick one BA1/BS1 frequency threshold, cap PP3 at Supporting/Moderate per ClinGen SVI guidance, and either validate `VARIANT_PRIORITY_SCORE` against a labeled truth set or relabel it clearly as an unvalidated heuristic in every report that shows it.
9. **Either implement or remove the `large_sv_vcf` mechanism** and disclose the pipeline's actual SV-handling status (none) in the README until it's real.
10. **Fix the hg19 liftover genotype-drop and add strand-flip allele re-normalization**, or explicitly document that hg19 VCF input with genotypes is not fully supported.

---

## Notes on method

All findings above were verified against the actual code (file:line references given inline) or by direct execution/reproduction (`python3 -c "import filter_variants"`, the bash quoting reproduction in C8, `grep` against the live BED panel for the substring-collision claims). Where a fork's hypothesis was checked and found **not** to be a bug (e.g., the consequence-filter substring matching itself, the `filter_freq` direction, `gene_priority` casing consistency, a suspected pandas boolean-alignment crash), it is omitted here rather than reported as a finding — this document only contains defects that survived verification.

---

## Resolution Log

### 2026-09-10 — C1, C2, C3, C3b fixed and verified end-to-end

**Safety commit first.** Before any edit, the uncommitted `modules/` refactor in `/home/mruizp/data_lab_PGP/shared/utils/` (C3b) was committed as a checkpoint (commit `419c1ba`, "wip: checkpoint in-progress modules/ refactor of filter_variants.py"), excluding the untracked 1.2 GB `pangolin_db/` reference file (now added to that repo's `.gitignore` instead). This gives the in-progress refactor a rollback point for the first time since 2026-09-01.

**Root cause, confirmed precisely.** `get_full_gene_curation_dataframe`, `filter_by_custom_gene_list`, `convert_to_parquet`, and `DuckDBVariantEngine` are imported into `shared/utils/src/filter_variants.py` but **never called anywhere in that file** (confirmed by grep — dead imports). Cross-checked against `/home/mruizp/data_lab_PGP/pipelines/clinical_variant_prioritization/` (a separate, independent repo, forked out of this shared library on 2026-08-31 per `TODO.md`'s "ACMG calculation externalized to downstream app pipeline" note): that repo has its own complete, working, independently-evolved copies of these exact names (`src/modules/disease_hpo.py:107,138` genuinely define `get_full_gene_curation_dataframe`/`filter_by_custom_gene_list`; `src/engine/` genuinely implements the Parquet/DuckDB engine). Confirmed that repo does **not** import Python source from `shared/utils/src/` (only reads reference *data* from `shared/utils/data/`), so it is provably unaffected by this fix. The four names were leftover references to functionality that migrated to `clinical_variant_prioritization` and was never fully cleaned out of `shared/utils/src/filter_variants.py`'s import list during the ongoing `modules/` refactor.

**Fix applied**:
- `shared/utils/src/filter_variants.py`: removed `get_full_gene_curation_dataframe, filter_by_custom_gene_list` from both the `.modules.disease_hpo`/`modules.disease_hpo` import blocks, and removed the two dead `from .engine import ...` / `from engine import ...` lines entirely (both the `try` and `except ImportError` branches).
- `annotation_pipeline_new/src/python/config.py`: removed the non-functional monkey-patch shim (the "in-memory fallback" TODO.md had marked done but which was never on the production import path — C2) that tried to patch around this exact ImportError.

**Verified end-to-end, not just at import time**:
```
$ /home/mruizp/apps/miniforge3/envs/datasci/bin/python3 -c "import filter_variants"   # exact prod interpreter, no config.py shim
OK: filter_variants imports cleanly
```
Then ran the actual two-stage `vcf2parsed.sh` production path manually against a real annotated test VCF (`test_data/test_run/annotation/MYBPC3.annotated.vcf.gz`, 21 variants):
- Stage 1 (`vcf_parser_pysam.py`, `vcf_parser` conda env): `Processed 21 variants. Variants skipped = 0` → 239-line intermediate TSV.
- Stage 2 (`filter_variants.py`, `datasci` conda env, `--skip_quality_filter` as production always passes): exit 0 → 238-row × 573-column Parquet output with `QC_STATUS` (`PASS`), `PRIORITY_TIER` (`Tier 3 (VUS / Moderate Potential)`), `NEW_IMPACT` (`MODERATE`) all correctly populated.
- `pytest tests/python/ -q` (this repo): still `34 passed` — the fix didn't break anything the existing suite covers, though per C3 that suite still doesn't exercise the real `vcf2parsed.sh` invocation directly (a gap that remains open — see Remediation item 1's CI suggestion).

**Not yet done** (left for a follow-up pass, since the user asked to fix the import breakage first): committing these two fixes to git (pending user confirmation on commit scope/message), and everything from C4 onward in this document — the filtering-module disconnection, inert QC gating, ACMG/scoring issues, predictor-parsing issues, and infrastructure findings are all still open.

*(Update: both fixes above were committed — `shared/utils` commit `4315e8d`, `annotation_pipeline_new` commit `59cce44` — after user confirmation.)*

### 2026-09-10 (continued) — C8 and C9 fixed; user decision recorded on default QC policy

**User decision**: default `main.sh` runs (no flags) must not apply any QC filter — this was already the de facto behavior (`--skip_quality_filter` always passed at the `vcf2parsed.sh` stage), so it's now confirmed as intentional rather than left as an unresolved question. QC filtering becomes available, opt-in, through two new `main.sh` flags.

**C9 (run-level filename collision, discovered while implementing the C4 fix — logged above) — fixed**: `filter_and_summarize.py` gained an `--output-prefix` argument (default `None`, fully backward compatible for any other caller); `main.sh` now passes `--output-prefix "$gene_name"`, so outputs are `filtered/<GENE>_filtered_variants.<ext>` and `filtered/<GENE>_filtering_summary_metrics.tsv` — parallel to the existing `results/<GENE>.parsed.clean.pq` convention, no more run-level collision across genes.

**C8 (shell-quoting bug) — fixed**: `main.sh`'s `filter_flags` changed from a plain string (which embedded literal `\"..\"` characters that survived into argv) to a bash array, expanded as `"${filter_flags[@]}"` at the `qsub` call site. Verified correct with the same reproduction technique used to originally catch the bug (a `py()`/nested-function stand-in for the qsub call, confirming clean argv splitting through an extra layer of indirection matching `qsub -b y bash run_python_hpc.sh ...`).

**New opt-in QC flags added to `main.sh`**: `--min-qual <FLOAT>` and `--pass-qc-only`, threaded through to `filter_and_summarize.py`'s existing (previously unreachable) `--min-qual`/`--pass-qc-only` arguments. Neither is set by default, matching the user's decision.

**Deliberately not done in this pass** (flagged to the user rather than decided unilaterally, per the advisor's review): whether/how the three report generators (`plot_annotation_results.R`, `generate_interactive_report.py`, `generate_clinical_prioritization_report.py`) should consume the now-correctly-named filtered output — the original C4 finding (filtering is invisible to every deliverable) is **not yet fixed**. Silently swapping every report's input source to the filtered file the moment any filter flag is passed would trade one silent behavior (flags do nothing) for another (reports quietly narrow with no on-report indication of what was excluded). This needs an explicit choice between: (a) reports always stay on the full unfiltered set, and the filtered file remains a correctly-produced side artifact only, or (b) reports switch to the filtered set only when filter flags were actually passed, with a visible banner stating which filters were applied and how many variants were excluded. Also flagged for a decision: the `--min-qual`/`--pass-qc-only` filters both let a variant with a **missing** `QUAL`/`QC_STATUS` value pass through (NaN-passthrough / not-`"LOW_QUAL"`-passthrough) — defensible for the evidence-based score filters (REVEL/CADD/etc., where "not scored" shouldn't count as "failed"), but arguably backwards for an explicit QC gate a user just opted into, where "we don't know this variant's QC status" is not the same claim as "this variant passed QC."

**Verified**: `filter_and_summarize.py` re-run directly against real output from the earlier C1 verification (238-variant `MYBPC3` parquet) —
- No flags: `238 / 238` retained (confirms default behavior is unchanged).
- `--pass-qc-only --min-qual 30`: `238 → 191` retained (confirms the new opt-in QC flags do real, meaningful filtering when requested).
- `--max-af 0.01 --consequences missense_variant,stop_gained`: `238 → 3` retained, output correctly named `MYBPC3_filtered_variants.pq`/`MYBPC3_filtering_summary_metrics.tsv` (confirms C8 and C9 fixes together).
- `bash -n main.sh`, `pytest tests/python/ -q` (34 passed), `bash tests/bash/test_cli_args_main.sh` (all passed) — no regressions.

**User decisions on the two open questions**:
1. **C4 report-wiring**: reports (R plots, interactive report, clinical prioritization report) will **continue to always show the full unfiltered set**. `filtered/<GENE>_filtered_variants.<ext>` remains a separate, correctly-produced (post-C9) side artifact for anyone who wants a filtered view — it is not wired into any report. This is now the pipeline's intended design, not an open defect; C4 is considered resolved by this decision (documented, not code-changed further).
2. **QC missing-data handling**: "these variants should pass the filter but be flagged as such"; if the whole dataset lacks QC values, that must be reported, not silently skipped; missing values must never be coerced to a numeric 0.

**Implemented for decision 2** — `filter_and_summarize.py`'s `filter_dataframe()`:
- `--pass-qc-only`: rows with a missing `QC_STATUS` still pass (not excluded), get a new `QC_STATUS_evaluated=False` column, and a `logger.warning` naming the count. Rows with a real, present `QC_STATUS == "LOW_QUAL"` are still correctly excluded — only the missing-data case changed.
- `--min-qual`: identical treatment via a new `QUAL_evaluated` column; `qual_vals.isna()` is used directly (never `.fillna(0)`) so a missing QUAL is never made to look like a confirmed low value. Rows with a real QUAL below threshold are still correctly excluded.
- If `QC_STATUS`/`QUAL` is **entirely absent** from the input (not just missing per-row), a `logger.warning` fires and `stats["QC_STATUS_column_present"] = False` / `stats["QUAL_column_present"] = False` is recorded in `filtering_summary_metrics.tsv` — this case is no longer silent.

**Verified** with a direct unit-level test (not just a smoke run): a 4-row synthetic frame with one row of all-missing QC data confirms real sub-threshold rows are still excluded, the missing-data row is retained and flagged `False` in both new columns, and a second frame with no `QC_STATUS`/`QUAL` columns at all correctly logs both warnings and records `False` in the stats without crashing or filtering anything. `pytest tests/python/ -q` still 34/34 and `bash -n main.sh` clean after this change.

**Status as of this update**: C1, C2, C3, C3b, C8, C9 fixed and verified; C4 resolved by explicit design decision (no code change needed beyond C9); the QC-missing-data handling requested alongside C5/C6 is implemented. **Not yet committed** — user asked to hold off on committing this round so the `main.sh`/`filter_and_summarize.py` diff can be reviewed first. Everything else in this document (C5/C6's remaining wiring already covered by the new opt-in flags; C7 genotype-extraction gap; the ACMG/scoring findings; predictor-parsing findings; orchestration/liftover findings) is still open.
