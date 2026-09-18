# 📋 Master Annotation Pipeline & Clinical Dashboard Roadmap

---

## 🔴 Priority 0: Foundational Fixes (block further predictor/parameter changes until done)

- [ ] **Run Provenance & Reproducibility Manifest** — confirmed gap, no code change yet
  - Verified 2026-09-15: there is currently **no way to reconstruct what actually ran for a past run**. Checked `main.sh`, `config/env.sh`, every `RUNS/<RUN_NAME>/` folder on disk, and `src/python/audit_run_results.py`:
    - No file anywhere records the CLI flags/thresholds a run was invoked with (`--build`, `--max-af`, `--min-revel`, `--consequences`, gene list, etc.).
    - No file records tool/model versions actually used (VEP version + cache version, CADD plugin version, REVEL/AlphaMissense version, SpliceAI package version + `-D` value, Pangolin model version, SPiP version, gnomAD release, reference FASTA/GENCODE build).
    - No git commit hash of the pipeline codebase is captured per run.
    - `config/env.sh` auto-detects/falls back across multiple binary locations (e.g. `PYTHON_SPLICEAI`, `BCFTOOLS_BIN`) with silent fallthrough — the actual binary resolved can differ between two runs of the same command on the same machine with no record of which one was used.
    - `audit_run_results.py`'s "version_lineage" only tracks output-file mtimes/re-run staleness, not tool/parameter versions — not a substitute for this.
  - This directly violates this project's own provenance requirement (every metric must be traceably tied to its source script, tool version, and database release) and is why the SpliceAI `-D` question below couldn't be answered by "check the last run" — there was nothing to check.
  - **Scope for the fix (not yet designed in detail)**: write a `run_manifest.json`/`.md` into each `RUNS/<RUN_NAME>/` at start of `main.sh` capturing: resolved CLI args, resolved binary paths from `env.sh` (not just env var names), pipeline git commit hash + dirty-tree flag, and per-predictor tool/model/database versions (queried at run time, not hardcoded docs — see existing `DOC-1` bug in `BUG_TRACKER.md` where docs already drifted from the real CADD version).
  - **Blocks**: the SpliceAI window fix immediately below (and any other predictor/parameter change) should not land until a run can at least record what value it used — otherwise this exact ambiguity (what did we actually run last month?) recurs immediately.
  - **Live example, 2026-09-16**: needed the original `--raw-dir` for `run_more_genes_20260817_1143` (to resume it after the `ORCH-10` VEP fix below) and it wasn't recorded anywhere — the user's own recollection (`_RAW/`) turned out to be empty and unrelated to this run. Only recovered it by brute-force forensics: found `RUNS/predictors_050826/input/by_gene_more_genes/` had the same gene count (154) as the run's audit report, then confirmed by opening `DMD.pq` and matching its row count (533,174) exactly against the converted `_tmp/DMD.vcf.gz`. That worked this time only because an unrelated run folder happened to still hold the input — it is not a real recovery mechanism and won't always work.

- [ ] **SpliceAI `-D` distance parameter: switch from `-D 10000` to `-D 4999`** (HIGH priority, not yet implemented)
  - Literature-reviewed 2026-09-15 (see chat/session notes; sources below) — this pipeline's current `-D 10000` is outside the range Illumina's own SpliceAI CLI accepts (`argparse choices=range(0, 5000)`, i.e. 0–4999 only). Our `src/python/annotate_spliceai.py` only reaches `10000` because it calls the library's `get_delta_scores()` directly instead of going through the CLI, bypassing that limit.
  - `-D` is not just a delta-*position* reporting radius: in `spliceai/utils.py`, `cov = 2*dist_var+1` and the delta-*score* argmax is taken over that same `cov`-length window — so a larger `-D` can change the reported score itself, not only where the tool says the site is. **Existing `-D 10000` output cannot simply be reused/truncated for a `-D 4999`-calibrated tool (e.g. SAI-10k-calc above) — a fresh `-D 4999` pass is required.**
  - `-D 4999` is independently used by: Moles-Fernández et al. 2021 (Cancers 13:3341, deep-intronic BRCA/Lynch screening), Canson et al. 2023/SAI-10k-calc (btad179), and the ClinGen SVI Splicing Subgroup calibration (Walker et al. 2023, AJHG) — the body that calibrates splicing evidence strength for ACMG/AMP.
  - No literature or tool precedent was found for a gene-specific/variable distance (e.g. larger window for TTN's large introns vs. small genes like CSRP3) — would be novel/unvalidated, not established practice, if pursued instead of a single fixed value.
  - **Not yet implemented per user instruction** — this changes the pipeline's default SpliceAI behavior pipeline-wide (all future runs, not just the SAI-10k-calc feature) and should wait until the provenance manifest above exists, so the change itself is actually recorded.

- [ ] **Scale every predictor job's `h_vmem`/thread count to gene size instead of one fixed value for every gene** (CRITICAL priority, not yet implemented — see `BUG_TRACKER.md` `ORCH-10`)
  - Originally scoped to VEP only; user expanded scope 2026-09-17 to cover all annotation tools (SpliceAI, Pangolin, SPiP, Branchpointer), not just VEP — `main.sh` hardcodes a fixed `h_vmem`/`-pe smp` for every one of these per-gene qsub calls (lines ~401-497), the same class of bug as `ORCH-10`.
  - Caused two real production failures on the same gene so far: `main.sh:442` submits every gene's VEP job with the identical `-l h_vmem=15G -pe smp 4` regardless of variant count. DMD (533,174 variants, 2.3x the panel's next-largest gene) was SIGKILLed once (job `4975105`, 2026-08-17, exit 137 at `maxvmem=56.885G`) and, resubmitted unchanged, then ran 19 days pinned at that same memory ceiling with zero output (job `4979674`, still stuck as of 2026-09-16).
  - Confirmed `h_vmem` is per-slot on this cluster (user confirmed 2026-09-17) — so `-l h_vmem=15G -pe smp 4` is really a 60GB total budget, consistent with the `56.885G` ceiling observed.
  - **Real historical data being gathered, not yet analyzed** (2026-09-18): `src/hpc/logger.sh`'s own per-step timing/memory instrumentation has never actually run (see `ORCH-11` — `LOG_DIR` is never exported by `main.sh`), so there is no local dataset to mine directly. Instead: extracted all 1,996 real SGE job IDs for every (gene, step) in `run_more_genes_20260817_1143` from local `_log/*/*.out`/`*.err` files into `resources/run_more_genes_20260817_1143_job_ids.tsv`, and wrote `resources/pull_qacct_snapshot.sh` to look up each job's real `maxvmem`/`cpu`/`ru_wallclock`/`slots` via `qacct -j` (the `-o <user>` broad-query form of `qacct` returned nothing on this cluster, likely an implicit date-window default — per-job-ID lookup works fine). **User will run this script; output lands at `resources/run_more_genes_20260817_1143_qacct.tsv`.**
  - **Important sequencing caveat (user, 2026-09-17)**: this snapshot reflects the *current* SpliceAI `-D 10000` behavior. Once the separate SpliceAI window fix above (`-D 10000` → `4999`) lands, SpliceAI's own memory/wallclock profile will change, making any SpliceAI-specific tiers derived from this snapshot stale. VEP/Pangolin/SPiP/Branchpointer tiers are unaffected by that change and can be derived from this snapshot as-is; SpliceAI's tiers should be re-derived from a fresh snapshot after the `-D` change lands, not before.
  - Fix scope (not yet designed): size `h_vmem` (and/or thread/fork count) per tool from the gene's input variant count at submission time, using real tiers derived from the `qacct` snapshot above rather than a fitted formula from too few points — or maintain a hardcoded override list for known outlier genes (DMD, RYR2, CACNA1C, LAMA2, NEBL, NF1, PRDM16, PCCA — all 5–20x the panel median).
  - **Blocked on the run-provenance manifest above** — a resource-sizing change like this is exactly the kind of per-run parameter that must be recorded going forward, or this incident repeats invisibly.

---

## 🚀 Phase 1: High-Performance Compute & Orchestration (Inference & Checkpointing)
- [x] **State-Aware Pipeline Checkpointing & Granular Flags**
  - [x] Implement `-s "$file"` file integrity checks across all pipeline stages (Steps 1–5)
  - [x] Add dynamic SGE dependency hold chaining (`active_predictor_jobs`)
  - [x] Add granular CLI override flags: `--overwrite-all`, `--force-spliceai`, `--force-vep`, `--force-pangolin`, `--force-spip`, `--force-branchpoint`, `--force-merge`, `--force-vcf2parsed`, `--force-reports`, `--skip-genes`
- [x] **Parallel SpliceAI VCF Chunking Engine**
  - [x] Implement streaming header-preserving VCF chunker (`src/python/split_vcf_chunks.py`)
  - [x] Integrate multi-worker parallel inference pool matching `$NSLOTS` in `src/hpc/annotate_spliceai_vars.sh`
  - [x] Add automated volume thresholding ($\ge 15,000$ variants) and individual chunk checkpointing
  - [x] Lossless `bcftools concat -a` merge and `tabix` indexing
  - [x] Unit test suite in `tests/python/test_vcf_chunking.py` (100% passed)
- [x] **Built-in Run Quality & Completeness Auditor**
  - [x] Implement `src/python/audit_run_results.py` with predictor completeness %, deep learning throughput, and schema diffing
  - [x] Integrate `--audit-run` and `--audit-only` CLI flags into `main.sh`
  - [x] Add automated Step 6 post-run SGE audit job (`audit_${RUN_NAME}`) chained to terminal tasks

---

## 🧬 Phase 2: Biological & Clinical Variant Prioritization
- [x] **gnomAD v4 Joint Population Frequency Integration**
  - [x] Fix UCSC `chr10` $\leftrightarrow$ Ensembl `10` contig mapping via `--synonyms chr_synonyms.txt` in `annotate_vep_vars.sh`
  - [x] Add `gnomADv4 AF grpmax` across all 4 dashboard tabs with fallback detection and warning alerts
- [x] **5' UTR Annotator & Multi-Tab Dashboard**
  - [x] Integrate 5' UTR consequence, uORF, and start/stop disruption annotations into `vcf_parser_pysam.py`
  - [x] Implement 4-tab Clinical Prioritization Report (`generate_clinical_prioritization_report.py`)
- [x] **Centralized Configuration & Shared Library Fallbacks**
  - [x] Implement `config/env.sh` and `src/python/config.py` path resolver
  - [x] Removed dead `get_full_gene_curation_dataframe`/`filter_by_custom_gene_list`/`engine` imports from `filter_variants.py` (leftover references to functionality that migrated to `clinical_variant_prioritization`; the previously-checked "in-memory fallback" in `config.py` was never on the production import path and did not actually fix this — see [walkthrough/20260910_full_technical_methodological_audit.md](./walkthrough/20260910_full_technical_methodological_audit.md))
  - [x] ACMG calculation externalized to downstream app pipeline; lightweight annotations retained in-pipeline
- [x] **Cardiomyopathy Domain & Cardiac Isoform / PSI Filter**
  - [x] Add cardiac ventricle Percent Spliced In (PSI > 85%) flag for `TTN` truncating variants (TTNtv in A-band)
  - [x] Map critical DCM structural domains (`BAG3` BAG domain, `LMNA` rod domain, `FLNC` Ig-like folds)

---

## 📊 Phase 3: Reporting & Clinical Dashboard Enhancements
- [x] **Multi-Gene Cohort Master Dashboard (`generate_cohort_master_dashboard.py`)**
  - [x] Aggregate top-tier prioritized variants across all single-gene Parquet tables into a unified master cohort dashboard
  - [x] Cohort-level summary metrics (distribution of High/Moderate candidates, gene breakdown, splicing vs missense burden)
  - [x] Interactive Plotly visualizations and 1-click CSV candidate exporter
  - [x] Comprehensive pytest coverage in `tests/python/test_generate_cohort_master_dashboard.py`

---

## 🔒 Phase 4: Containerization & Infrastructure Stability
- [x] **Unified Apptainer / Singularity SIF Container (`annotation_pipeline.sif`)**
  - [x] Create definition file `resources/containers/annotation_pipeline.def` packaging Python datasci, R, genomics binaries, and Plotly
  - [x] Automated build script `src/hpc/build_container.sh` with environment-resolved output paths
  - [x] Integrated `--use-container` and `--sif <PATH>` CLI flags in `main.sh` and execution wrapper in `config/env.sh`

---

## 🧬 Phase 5: Future Predictor Integrations
- [ ] **SAI-10k-calc Integration (SpliceAI-10k calculator, Canson et al. 2023, Bioinformatics btad179)**
  - Deterministic post-processor over raw SpliceAI delta scores: predicts splicing aberration type (pseudoexonization, partial/whole intron retention, partial exon deletion, (multi)exon skipping), aberration size (bp), reading-frame effect, and predicted amino acid sequence.
  - Belongs in the annotation pipeline (per-variant, deterministic, computed once — same shape as SPiP/Pangolin/Branchpointer), not in `clinical_variant_prioritization`; the app would just read the resulting columns.
  - **Inputs already available**: `SpliceAI_pred_DS_AG/AL/DG/DL` + `DP_AG/AL/DG/DL` from `src/python/annotate_spliceai.py`; RefSeq NM_ transcript IDs already in `resources/gene_transcript_mapping.txt`.
  - **Open items before implementation**:
    - Window mismatch: confirmed (see **Priority 0** at the top of this file) — SAI-10k-calc's thresholds were calibrated on `-D 4999` output, and `-D` changes the reported delta *score*, not just its position, so this needs the pipeline's `-D 4999` fix landed first, then its own dedicated pass.
    - SNV-only tool; this pipeline also annotates indels — new columns need an explicit `not_applicable` status-contract value for indels (not blank/0, to avoid the PRED-6-style "not scored" vs "confidently benign" conflation).
    - Wrap the published R implementation (https://github.com/adavi4/SAI-10k-calc) as a subprocess rather than reimplementing — the real decision logic is in a supplementary flowchart (Supplementary File S1) not available locally; porting it from the paper's summary text would mean guessing at the algorithm.
  - Source paper saved at `/home/mruizp/Downloads/btad179.pdf` (not yet copied into the repo).

