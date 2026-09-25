# Pangolin tissue-of-origin patch (real WGS run, S223, in progress)

**2026-09-25**

## Context

Mid-way through the real WGS run (`S223`), stopped to answer a scientific-validity question before
touching anything else: is Pangolin actually being used in a way that's meaningful for a
cardiovascular-disease platform?

## FACT — what Pangolin's model actually computes

The vendored model files (`src/external/Pangolin-main/pangolin/models/final.<j>.<i>.3.v2`) are
exactly 4 GTEx tissues × 3 replicate models = 12 files. Tissue order confirmed against this same
vendored repo's `scripts/custom_usage.py:6-13` header comment: **Heart, Liver, Brain, Testis** — no
left-ventricle/atrial-appendage split, no other tissues. `pangolin.py`'s `compute_score()` computed
all 4 tissues' scores, then folded them via `np.argmin`/`np.argmax` **across all 4 tissues** into a
single loss/gain value per position (`pangolin.py:64-65`, pre-patch) — the reported score was "the
most extreme effect across any of the 4 tissues," with the winning tissue's identity computed
internally and then discarded before ever reaching the output.

`modules/local/pangolin.nf`'s invocation (`-m` default, no `-s`/`--score_exons`) only ever exercises
this default-fold branch — confirmed by reading the process script, no flags override it.

## FACT — a separate, pre-existing issue found while checking this (not fixed here)

`shared/utils/src/modules/splicing.py:34-38`'s `NEW_PANGOLIN_COLUMNS` already has
`Pangolin_heart_lv_score`/`Pangolin_heart_aa_score` — names implying real GTEx heart-subtissue
scores. `parse_pangolin` (`splicing.py:296-297`) writes the *same* generic any-tissue max/min value
into all three columns (`max_score`, `heart_lv`, `heart_aa`) — they're identical numbers under
different names, not derived from any tissue-specific computation (the underlying model has no
LV/AA split to derive them from in the first place). Not a data-loss bug — flagged, not fixed,
per Mario 2026-09-25: `carve-platform` `PLT-138` / `PRED-14` (renumbered from an initial `PRED-13` —
that ID turned out to already be claimed in this repo's own `BUG_TRACKER.md` for an unrelated
UTRAnnotator finding from concurrent work found while checking `git status`; not touched, just
avoided).

## What was built

Patched the vendored `src/external/Pangolin-main/pangolin/pangolin.py` (this repo's own copy, not
a shared/utils file) to track and report which of the 4 *real* tissues produced the winning
gain/loss score, additively:

- `compute_score()` now also returns `loss_tissue`/`gain_tissue` (per-position `argmin`/`argmax`
  index into the 4-tissue axis, computed before it gets folded away) alongside the existing
  `loss`/`gain` arrays.
- `process_variant()`'s default branch (the only one this pipeline exercises) builds a new
  `tissue_list` in parallel to the existing `scores_list`, one `gene|gain_tissue|loss_tissue` entry
  per gene.
- **The existing `Pangolin=` INFO field's format is byte-for-byte unchanged** — `scores_list`'s
  construction was not touched. The new information goes into a **separate new INFO field,
  `PangolinTissue`**, specifically so `shared/utils`'s `parse_pangolin` (which requires exactly 2
  `:`-delimited parts per entry) keeps working unmodified. A 3rd `:`-delimited component was
  considered and rejected — it would have silently broken that parser (`len(sub_parts) == 2` check
  would fail, `try/except ValueError` would swallow the resulting float-parse failure, and every
  variant would report `not_covered`).
- `main()`'s VCF path writes the new header line and field; the CSV path (never used by this
  pipeline, confirmed via `pangolin.nf`'s always-`.vcf` invocation) was kept consistent with the
  changed `process_variant()` return contract (now a tuple, not a bare string) so it isn't left
  silently broken for any other caller.
- `modules/local/merge.nf`: `bcftools annotate -a ... -c Pangolin` was **whitelist-only** — the new
  tag would have been silently dropped at the merge step. Changed to
  `-c Pangolin,PangolinTissue`.
- No change needed in `vcf_parser_pysam.py`: `info_columns = [k for k in vcf.metadata.info_dict if
  k != "CSQ"]` (`vcf_parser_pysam.py:201`) is a generic INFO-field pass-through — any new tag that
  survives the merge automatically becomes a table column with zero additional code.

## Verification

No GPU/model weights available in this sandbox (same constraint as every other Pangolin/SpliceAI
check this whole migration) — real inference can't run here. What was verified:
- `python3 -m py_compile` on the patched file.
- A standalone synthetic numpy test (4×5 fake tissue-score matrix, no torch/model dependency)
  confirming the `argmin`/`argmax`-based tissue selection matches an independently-computed
  column-wise max/min, and that the winning tissue label always corresponds to the actual winning
  value.
- Direct read-through confirming `scores_list`'s construction (the existing `Pangolin=` field) is
  untouched by the diff — only new `tissue_list` entries were added alongside it.
- Traced the new field's full path to the final table (`merge.nf`'s `-c` list →
  `vcf_parser_pysam.py`'s generic pass-through) to confirm it won't be silently dropped anywhere
  in between.

**Not yet verified**: real end-to-end behavior against actual model inference — needs a real
cluster run. This is the same class of residual risk every predictor change in this migration has
had until a real run confirmed it.

## What's next

One more thing is needed before relaunching `S223` — see the follow-up section below for the
`chunk_0017` crash (now fixed) and the wall-clock kill issue (not yet fixed).

## `chunk_0017` root cause found and fixed — real IUPAC ambiguity code, real GRCh38 primary chromosome

Pinpointed the exact trigger by replaying Pangolin's own reference-fetch logic (pure pysam, no
torch/model dependency, no GPU needed) against `chunk_0017.vcf.gz`'s actual records starting from
line 12060 (the point right after the last "skipping variant" warning in the failed task's log):
**`chr16:88373021 C>T`** (in `ZNF469`'s gene body). Its ~30kb fetch window (`pos ± (5000+10000)`,
matching `process_variant`'s exact fetch call) contains an `R` at `chr16:88366825`, 6,196bp
upstream — confirmed independently via `samtools faidx chr16:88366820-88366830` → `CCATCRCCATA`,
not an artifact of the replay script.

This is a **real, primary-chromosome** ambiguity code, not a decoy/ALT/HLA-contig artifact as
initially guessed (see the earlier section above) — the human reference genome does legitimately
carry a small number of true IUPAC-ambiguous positions even on primary chromosomes, inherited from
the original assembly process. Pangolin's `one_hot_encode` (`pangolin.py:23-31`, pre-fix) only
translated `A/C/G/T/N`; any other character (here, `R`) reached `int()` untranslated and crashed.

**Fix**: any character outside `A/C/G/T` now maps to the same `'0'` code `N` already used — i.e.
the same all-zero "no information" one-hot row, not a crash. Verified with the actual crashing
11bp window (`CCATCRCCATA`): encodes correctly on both strands (including through the
reverse-complement transform), and produces output *identical* to what the same window would
produce if the `R` were literally an `N` — confirming this is a direct, non-regressing extension
of the code's own existing convention, not new behavior invented for this case.

**Not verified**: whether any *other* chunk of the real WGS input hits a different non-ACGTN
position — the fix is generic (handles any character, not just `R` at this one position), so it
should cover those too, but that's only confirmed on a real cluster run.
