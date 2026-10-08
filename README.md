# Universal Genomic Variant Annotation, Filtering, and Visualization Pipeline

The active annotation workflow is the Nextflow whole-VCF pipeline ([annotate_vcf.nf](./annotate_vcf.nf)), which annotates GRCh38 panel, WES and WGS VCFs. See [Nextflow Whole-VCF Pipeline](#-nextflow-whole-vcf-pipeline) below.

> **Legacy cleanup — 2026-10-08T11:53:30+02:00:** the old Bash pipeline, report/dashboard scripts, downstream R analyses, their tests and the old general-purpose container were moved to `trash/legacy_pipeline_20261008_115330/`. Nothing was deleted. [Move manifest](./reports/legacy_cleanup_20261008/manifest.json) records original and current paths. Active annotation uses Nextflow; SPiP remains an R dependency.

> [!WARNING]
> **Every Pangolin score produced before 2026-09-30 is invalid.** The old shared database `shared/utils/pangolin_db/pangolin_grch38.db` is a mouse annotation (GENCODE M23), not human GRCh38 (`BUG_TRACKER.md` `MISC-12`, `carve-platform` `PLT-139`). It controls both gene assignment and exon masking. Both pipelines now expect a GENCODE 45 rebuild, `shared/utils/pangolin_db/gencode.v45.ensembl_canonical.grch38.db`, and fail at the Pangolin step until it is deployed. Re-run Pangolin on earlier results before using them ([walkthrough](./walkthrough/20260930_pangolin_db_gencode_v45_rebuild.md)).

---

## 🚀 Key Features & Integrated Predictors

1. **Universal Input & Coordinate Standardization**:
   - **cDNA / Transcript notation**: E.g. `NM_000257.3:c.526C>T`, `ENST00000343260:c.100A>G`, `c.1504C>T`.
   - **Genomic coordinates**: Standard `CHROM`, `POS`, `REF`, `ALT` in **hg19/GRCh37** (with automatic liftover to GRCh38) or **GRCh38**.
   - **File formats**: Direct `.vcf`, `.vcf.gz`, `.xlsx`, `.csv`, `.tsv`, `.pq`, or `.parquet`.
   - **Signed Intron Offset & Splice Side**: Calculates exact signed distance to nearest exon boundary (`-12` = 12bp intronic, `+5` = 5bp exonic) and splice side (`donor`, `acceptor`, `exonic`).

2. **State-of-the-Art Splicing & Functional Annotation Matrix**:
   - **Ensembl VEP 111 (GRCh38)**: CADD (v1.6), REVEL, AlphaMissense, UTRAnnotator, MaxEntScan, and gnomAD v4.1 population frequencies (`AF_joint`).
   - **Population fields**: custom gnomAD v4 frequencies, allele counts/denominators and alternate-homozygote counts (`gnomADv4_nhomalt_*`) are retained. Cached gnomAD v2/v3, VEP `MAX_AF*` and homozygote age histograms are excluded from new runs. 1000 Genomes fields remain. Mastermind and LoFtool remain enabled; dbNSFP requests all 4.8a fields except duplicate `TSL`/`APPRIS`, which come from VEP ([field list](./resources/vep_dbnsfp_fields.txt)).
   - **Population AF interpretation**: custom `gnomADv4_AF_grpmax_joint` falls back per row to `gnomADv4_AF_joint`. An allele with no gnomAD v4 record (`absent_from_gnomad`) earns PM2_Supporting and counts as rare, while `POPULATION_AF` stays null; site coverage is not checked. Unknown or invalid AF (no gnomAD annotation, or a record with `AN_joint` 0) earns no frequency evidence or rarity bonus; reported zero stays distinct. Tables expose `POPULATION_AF`, `POPULATION_AF_SOURCE` and `POPULATION_AF_STATUS` (`available_frequency`, `absent_from_gnomad`, `missing_frequency`, `invalid_frequency`). Rarity filters retain unknowns for review; they do not establish adequate population coverage. Older or generic AF columns are used only through an explicit frequency-column override (`--freq_col` in the parser; the standalone filter is archived).
   - **SPiP v2.1**: Multi-threaded empirical and machine learning splicing impact predictor (up to 12 cores) with decoded mechanisms (`Exon_skipping`, `Donor_disruption`, etc.).
   - **Pangolin Splicing Predictor**: PyTorch deep learning ensemble predicting splice site gain/loss within 4999 bp of each variant (`-d 4999`), with GENCODE 45 Ensembl_canonical gene models. Output produced before 2026-09-30 is invalid: see the warning above (`MISC-12`). Pangolin has one generic heart model, not separate left-ventricle/atrial-appendage models, so the downstream `heart_lv`/`heart_aa` columns are not subtissue-specific (`PRED-1`, `PLT-138`).
   - **Local SpliceAI at `-D 4999`**: Deep residual neural network predicting max delta scores within 4999 bp of each variant (the maximum the Illumina CLI accepts).
   - **SpliceAI source selection**: VEP lookup and custom prediction fields are exported separately. A valid custom score takes precedence, including zero; valid VEP is the fallback. `SPLICE_MAX_UNIFIED` stays null when neither is available, with `SPLICEAI_METHOD = Not Covered` and availability in `SPLICEAI_STATUS`.
   - **Branch Point Predictor Pair**:
     - **LaBranchoR**: Bidirectional LSTM predicting top human catalytic branch points across GRCh38 (206,249 branch points).
     - **Branchpointer**: Signal-feature model estimating branch point probability, U2 snRNA duplex free energy (kcal/mol), and motif disruption flags (`YES`/`NO`).
   - **SpliceVault**: Direct empirical RNA-seq aberrant event lookup across GTEx/SRA datasets with human-readable decoded events and sample counts.
   - **SpliceVarDB**: Direct lookup of experimentally validated splicing variants from minigene, RT-PCR, and RNA-seq assays.
   - **ClinVar** (VEP `--custom`, `ClinVarVCF*` columns): germline classification, review status and disease names for SNVs and indels from a dated ClinVar GRCh38 VCF (`clinvar_vcf` in `nextflow.config`, stored under `/data_lab_PGP/resources/annotation/clinvar/GRCh38/`; rebuild with `resources/build_clinvar_vcf.sh`).
   - **ClinVar source priority**: interpretation prefers a nonmissing `ClinVarVCF_CLNSIG`, then `clinvar_clnsig` (dbNSFP), then `CLIN_SIG` (VEP cache), then legacy aliases. Selection is per row; VUS/conflicting classifications do not trigger fallback. `CLINVAR_SOURCE`, `CLINVAR_SELECTED_CLNSIG`, `CLINVAR_SELECTED_REVIEW` and `ClinVar_Stars` expose the selected evidence. Review metadata is never borrowed from another release.

3. **Unified Status Contracts & Data Schema**:
   - All predictors adhere to strict, unified status contracts: `scored` (valid numerical prediction), `not_covered` (outside gene model/receptive field), or `error`.
   - For complete column documentation, see [resources/column_dictionary.md](./resources/column_dictionary.md).

4. **Parser Filtering and Annotation Interpretation**:
   - The retained parser supports explicit frequency filtering (`--freq`, `--freq_col`). The former standalone multi-predictor filter and its CLI are archived; these are not Nextflow launch options.

5. **Legacy Visualization & Clinical Report Suite (outdated)**:
   - **Static Publication Plots (R ggplot2)**: Consequence breakdown barplots, REVEL vs AlphaMissense correlation scatter plots, and AF spectrum histograms saved as PDF and PNG.
   - **Interactive 4-Tab Clinical Prioritization Dashboards (`generate_clinical_prioritization_report.py`)**:
     - **Multi-Track SVG Transcript Visualizer & Variant Position Map**: Interactive Plotly multi-track diagram mapping genomic/cDNA positions and amino acid coordinates across exon-intron boundaries, displaying priority tiers, SpliceAI Δ, and AlphaMissense scores.
     - **gnomAD v4.1 Joint PopMax Frequencies**: Full extraction of `gnomADv4_AF_grpmax_joint` with population ancestry labels (`AFR`, `AMR`, `EAS`, `NFE`, `SAS`, etc.).


---

## 📁 Repository Architecture & Codebase Organization

### Historical Execution Directory (`RUNS/<RUN_NAME>/`)
All execution outputs, intermediate files, logs, and figures are organized within a dedicated, isolated run folder inside `RUNS/`:

```
RUNS/<RUN_NAME>/
├── _tmp/                                     # Intermediate VCF files from variant conversion
├── annotation/                               # VEP, SPiP, Pangolin, SpliceAI, Branchpoint annotated VCFs
├── results/                                  # Final parsed clean tables (.pq, .tsv, .xlsx)
├── filtered/                                 # Filtered variant tables & summary metrics
├── plots/                                    # Publication-grade static plots (PDF & PNG)
├── reports/                                  # Standalone interactive HTML dashboards
└── _log/                                     # Standard error and stdout logs per step and gene
```

### Active Source Modules

- `src/python/parsing/`: versioned VCF-to-table parser and annotation interpretation.
- `src/python/`: predictor runners, schema/provenance, genotype extraction and reference builders.
- `src/tools/`: retained audit, reference and analysis utilities.
- `src/external/`: required third-party Pangolin code and model weights.
- `modules/`, `workflows/`, `annotate_vcf.nf`, `nextflow.config`: Nextflow orchestration.

Legacy code is retained under the dated trash directory listed above. Its scripts expect their original paths; use the move manifest to restore them before attempting any historical reproduction. Existing `RUNS/` outputs and prior `archive/` content remain in place.

---

## 🧬 Nextflow Whole-VCF Pipeline

[annotate_vcf.nf](./annotate_vcf.nf) takes one GRCh38 VCF of any size, runs every predictor across it, merges the results, and writes one final table. It is intended for Module 1 (`sarek_pipeline`) output.

**Validation so far:**
- **7-gene panel (`panel7_test.vcf.gz`):** 72/72 records match `main.sh`. The only remaining difference is SPiP float formatting (`PLT-045`, [walkthrough](./walkthrough/20260924_nextflow_migration_item11_parity_verdict.md)).
- **First real WGS run (`S223`):** completed on 2026-09-29 with exit code 0. Every predictor's reassembled output matches its input record count (5,085,977 broad-tier, 936,412 gene-restricted).

This doesn't establish that the scores are correct; see [Known limitations](#known-limitations).

### Running it

Launch through the wrapper, which writes a provenance manifest before and after the run and passes every argument to `nextflow run annotate_vcf.nf` unchanged:

```bash
bash run_annotate_vcf.sh -profile standard \
  --input_vcf /data_lab_PGP/.../<sample>.vcf.gz \
  [--run_id <label>] [--spip_tier restricted|broad] [--output_format pq|tsv] [-resume]
```

- **Where to launch:** from a host that can submit SGE jobs. The `nextflow` head process stays alive for the entire run; `S223` (WGS) took about 4 days (2026-09-25 → 2026-09-29), so use a session that survives disconnection.
- **Work directory:** new `standard` launches use `/mnt/scratch/mruizp/annotation_pipeline_new/work`. Published results remain under `nf_work/annotation_out/`; `local_dev` retains its default work location. The scratch path must be writable and visible to the cluster execution nodes. Existing work directories are not migrated.
- **Resuming older runs:** use the original launch directory and add `-work-dir /data_lab_PGP/pipelines/annotation_pipeline_new/work -resume` for runs whose tasks are stored in the project’s old `work/`. Preserve both the corresponding work files and `.nextflow/cache`.
- **Finding Nextflow:** the wrapper uses `nextflow` from `PATH`, falling back to `/opt/nextflow/nextflow`.

| Argument | Default | Description |
| --- | --- | --- |
| `--input_vcf <FILE>` | required | bgzipped VCF with a sibling `.tbi` index. **Must be GRCh38**, called against the GATK-bundle `Homo_sapiens_assembly38.fasta` that the predictors use. A pre-flight check fails the run before any predictor starts if the `##contig` lengths or the first 1000 REF alleles disagree with that FASTA. Without `##contig` lines it continues with a warning only if ≥20 REF alleles match (`MISC-14`). Both `chr1` and `1` contig names are accepted. The same check also fails the run on multiallelic records, records that are not left-aligned or parsimonious, gVCF input and duplicate sample names (fix: `bcftools norm -m -any -f <fasta>`); duplicate sites and ALT equal to REF only warn |
| `--run_id <LABEL>` | input filename without `.vcf.gz` | Output folder name and prefix for every published file |
| `--spip_tier` | `restricted` | Which tier SPiP runs in (see below) |
| `--output_format` | `pq` | Final table format: `pq` or `tsv` only |
| `--genotypes` | `auto` | `auto` writes the long genotype table when the input has samples with GT; `off` skips it |
| `--genotype_wide` | `false` | Also write `<run_id>.parsed.clean.withGT.<fmt>`: the main table plus `GT_<sample>` (HET/HOMALT/HOMREF/MISSING) and `GQ_/DP_/AD_<sample>` columns. Parquet only; joined in batches; fails above `--genotype_wide_max_samples` (default 10) or if any table row has no genotype row. On a WGS table (S223: 22.6 M rows × 591 columns) the output is as large as the final table |
| `-profile` | `standard` | `standard` = SGE + Singularity (production). `local_dev` = local miniforge environments, for `-stub-run` testing only |
| `-work-dir <PATH>` | scratch under `standard` | Override the task work directory; use the original directory when resuming an older run |
| `-resume` | off | Reuse cached tasks from the selected work directory and launch-directory cache |

### Two predictor tiers

| Tier | Records | Tools |
| --- | --- | --- |
| Broad | every record in the input | VEP 111 (with its plugins), Branchpointer/LaBranchoR |
| Gene-restricted | only records inside [`resources/v5_genes_loc.bed`](./resources/v5_genes_loc.bed) (`Complete_gene_list_V5` regions) | SpliceAI (`-D 4999`), Pangolin (`-d 4999`), SPiP (unless `--spip_tier broad`) |

Records outside the gene BED have **no** SpliceAI/Pangolin/SPiP annotation. A missing score there means the tool was not run, not that the variant is predicted benign.

### Outputs

Everything is published to `nf_work/annotation_out/<run_id>/`, in sub-folders (paths below are relative to it):

| File | Content |
| --- | --- |
| `reports/<run_id>.assembly_check.tsv` | Pre-flight input check result: contig-length and REF-allele comparison against `params.fasta`, record normalisation counts, genotype mode (`no_samples`, `samples_without_GT`, `single_sample`, `multi_sample`), sample names, declared FORMAT fields, plus any warnings |
| `pipeline_info/RUN_MANIFEST.json` | Launch/completion provenance: command line, git commit and dirty-tree status, resolved binaries, input/resource hashes, database releases read from the databases themselves (`databases`), tool versions (VEP image and version, env fingerprints), exit code |
| `tables/<table>.pq.provenance.json`, `genotypes/<table>.pq.provenance.json` | Written next to each published table when the run exits 0: run id, code commit, parameters, database releases, tool versions, and the full manifest's path and sha256. Never overwritten |
| Parquet footer key `carve.provenance` | The launch-time part of the same provenance (`pipeline_info/TABLE_STAMP.json`) stored in the file metadata of every published table, so it travels with a copied or uploaded `.pq`. Read with `pq.ParquetFile(f).metadata.metadata[b"carve.provenance"]` |
| `pipeline_info/RUN_DIRTY.<timestamp>.patch` | Only for a launch from a dirty tree: the uncommitted changes. Rebuild the launch-time code with `git checkout <commit> && git apply <patch>` |
| `predictors/<run_id>.annVEP.vcf.gz`, `.annBranchpoint.vcf.gz`, `.annSpliceAI.vcf.gz`, `.annPangolin.vcf.gz`, `.annSPiP.vcf.gz` | Per-predictor annotated VCFs (+ `.tbi`) |
| `annotation/<run_id>.annotated.vcf.gz` | All predictors merged into one VCF (+ `.tbi`) |
| `tables/<run_id>.parsed.clean.pq` | **Final table.** One row per variant × VEP transcript. Every transcript is kept and tagged with `TRANSCRIPT_PRIORITY_TIER`: 1 = curated 217-panel transcript, 2 = MANE Select, 3 = other |
| `genotypes/<run_id>.genotypes.pq` | Long genotype table, one row per variant × sample (every pair, incl. homozygous-reference and no-call), joinable to the final table on `Locus`. Written only when the input has samples with a `GT` field. A pair that is absent means no evidence, not homozygous reference. Columns in `resources/column_dictionary.md` §4 |
| `tables/<run_id>.parsed.clean.withGT.pq` | Only with `--genotype_wide`: the final table with `GT_<sample>` columns added (never replaces the final table) |

| `pipeline_info/execution_trace_<ts>.txt`, `execution_timeline_<ts>.html`, `execution_report_<ts>.html`, `nextflow_<ts>.log` | Nextflow's own per-task trace (duration, requested memory, peak RSS, node), timeline, report and log, one timestamped set per launch. The log is written by `run_annotate_vcf.sh`; a direct `nextflow run` leaves it in the repo's `.nextflow.log` |
| `reports/<run_id>.run_qc.<timestamp>.html` (+ `.json`, `.processes.tsv`) | Run QC report written by `src/python/build_run_qc_report.py` (the wrapper runs it after Nextflow ends; it never overwrites an earlier report): run summary, record counts per stage (input, each predictor VCF, merged VCF, table rows, distinct `Locus`, genotype rows), predictor `*_status` and `TRANSCRIPT_PRIORITY_TIER` counts, time / memory / nodes per process, failed or retried attempts with the last error line, provenance. A figure that cannot be determined is shown as `unavailable`, never as 0 |

Runs made before 2026-10-05 have every file directly in `<run_id>/`; the manifest reader, the QC script and `write_run_manifest.py` accept both layouts. The QC script also runs by hand on any run folder, with `--nextflow-log` when there is no trace: `python src/python/build_run_qc_report.py --run-dir nf_work/annotation_out/<run_id> [--nextflow-log .nextflow.log]`. Do not `-resume` a run of the flat layout with the new code: cached outputs are re-published into the sub-folders, so every large file would exist twice.

For scale: on `S223` the final table was 5.0 GB and the merged VCF 2.7 GB. Intermediate task files live in Nextflow's `work/` directory. On a WGS run this directory holds days of compute; deleting it forces a full recompute.

The manifest (and any dirty-tree patch) stays inside the run folder (`pipeline_info/`), so it moves with the run if the folder is relocated. The wrapper does not add run records to this repository.

### Known limitations

- **Pangolin output before 2026-09-30 is invalid:** the old reference database was a mouse annotation (`MISC-12` / `PLT-139`). The GENCODE 45 rebuild must be deployed to `shared/utils/pangolin_db/` before the Pangolin step can run, and it is not yet validated on the cluster.
- **Table built before 2026-10-01 may carry another gene's or transcript's predictor scores (`PRED-14`):** SpliceAI, Pangolin and SPiP values were attached to rows (variant x transcript) by position/symbol, not by the row's own gene/transcript. The matched parsers (`src/python/parsing/`, commit `ae5a38c`) need the HGNC complete set deployed at `/data_lab_PGP/resources/annotation/hgnc/hgnc_complete_set.2026-09-29.txt` (`params.hgnc_table`); the Nextflow `VCF_TO_TABLE` step fails without it by design. New per-row `*_match` status columns say why a score is empty.
- **SPiP only scores transcripts in its own RefSeq database (early-2020 snapshot):** newer accessions and some MANE transcripts get no strict SPiP result (about 79% of S223 rows are empty; labelled gene-level `SPiP_samegene_*` columns exist). Rebuilding the database is deferred (walkthrough `20261001_predictor_gene_transcript_attribution_audit.md`, Addendum 9).
- **Derived identity resources must be rebuilt when their inputs change** (inputs and hashes are in each `*.manifest.json`): `resources/spliceai_symbol_to_ensg.*.tsv` (`src/python/build_gene_identity_tables.py`: SpliceAI annotation + GENCODE gene models), `resources/enst_to_spip_nm.*.tsv` (`src/python/build_enst_spip_map.py`: GENCODE GTF + SPiP RefSeq database + HGNC; rebuild if SPiP's database or the HGNC file changes), `resources/dbnsfp_transcript_aligned_columns.txt` (`src/tools/derive_dbnsfp_dependency.py`: dbNSFP readme + raw lines; rebuild when dbNSFP is upgraded).
- **Genes missing from the gene BED:** 144 curated genes, including panel gene `TAZ`, get no gene-restricted predictors (`MISC-4`).
- **No cross-predictor record-count checks:** counts aren't compared between predictors or against the merge output (`MISC-6`).
- **Dirty-launch provenance has limits:** a launch from a dirty tree saves `RUN_DIRTY.<timestamp>.patch` next to the manifest (`MISC-13`). Untracked files larger than 5 MB, or outside the code paths (`*.nf`, `*.config`, `*.sh`, `modules/`, `workflows/`, `src/`, `config/`, `resources/`), are listed but not captured. Launching from a clean, committed tree is still preferred.
- **Assembly check samples, it doesn't exhaustively verify:** the pre-flight check compares header contig lengths and only the first 1000 REF alleles (`MISC-14`). VEP runs without `--check_ref`, so per-record REF disagreements beyond that sample are not detected.
- **Unprofiled resource requests:** CPU/memory/time are not yet sized from profiling data. Per-chunk wall-clock varied from minutes to more than 46 h on `S223` (`TODO.md`, Priority 0).

`main.nf` is the older per-gene Nextflow adapter. It is kept working but is no longer developed.

---

## 📊 Pipeline Outputs & Data Schema

All output tables contain standardized, clean columns. For the full column dictionary, renamed field mappings, and status contracts (`scored`, `not_covered`, `error`), see [resources/column_dictionary.md](./resources/column_dictionary.md).

Current runs publish beneath `nf_work/annotation_out/<run_id>/` into `annotation/`, `predictors/`, `tables/`, `genotypes/`, `reports/` and `pipeline_info/`. See the Nextflow Outputs section above. Historical `RUNS/` outputs remain in place.

---

## 📜 HPC Environment & Container Options

### Active VEP Container

Nextflow uses `params.vep_sif` from `nextflow.config` (the separate VEP image under `/data_lab_PGP/resources/sif_images/vep.sif`). The former general-purpose image and definition are archived in the dated trash directory; they are not the active VEP container.

### Option 2: Conda Environments (`/data_lab_PGP/shared/utils/conda_envs/`)
- **Data Science & Post-Processing**: `datasci/`
- **Genomics & VCF Tools**: `genomics/`
- **Coordinate Conversion & Liftover**: `liftover/`
- **SPiP v2.1 Predictor**: `spip_env/`
- **Pangolin Predictor**: `pangolin_env/` (`resources/conda/pangolin_env.yml`)
- **SpliceAI Predictor**: `spliceai_env/` (`resources/conda/spliceai_env.yml`)

---

## 🧪 Master Unit Testing Suite

Run retained Python test modules independently with the local-compatible runner:

```bash
bash tests/run_all_tests.sh
```

The runner isolates each module's Python imports to avoid the existing `modules` package collision. It no longer invokes retired R plotting or Bash launcher/container tests. A single test module can be run with `python3 -m pytest -q tests/python/<test_file>.py`. This checks retained Python logic; real annotation still needs the small cluster fixture run.

---

## 📖 Walkthrough Logs & System Troubleshooting

- [2026-08-05: Pipeline Initialization](./walkthrough/2026-08-05_pipeline_initialization.md)
- [2026-08-06: NFS Storage I/O Error Investigation & Fix Guide](./walkthrough/2026-08-06_nfs_io_error_investigation_and_fix.md)
- [2026-08-06: SGE Variant Converter Freeze Investigation & Fix](./walkthrough/2026-08-06_sge_varconv_freeze_fix.md)
- [2026-08-07: Predictors 050826 Run Audit & Verification](./walkthrough/20260807_run_audit_predictors_050826.md)
- [2026-08-07: Splicing Predictors Implementation Roadmap & Priority Guide](./walkthrough/20260807_splicing_tools_implementation_guide.md)
- [2026-08-07: Objective 1 Completion - SPiP v2.1 Enhancements](./walkthrough/20260807_objective1_spip_enhancements.md)
- [2026-08-07: Objective 2 Completion - SpliceVault Empirical RNA Integration](./walkthrough/20260807_objective2_splicevault_integration.md)
- [2026-08-07: Objective 3 Completion - Signed intron_offset & Splice Side Refactor](./walkthrough/20260807_objective3_signed_intron_offset.md)
- [2026-08-07: Objective 4 Completion - Pangolin Splicing Predictor Integration](./walkthrough/20260807_objective4_pangolin_integration.md)
- [2026-08-07: Pipeline Unit Testing Suite Implementation & Execution](./walkthrough/20260807_unit_testing_suite.md)
- [2026-08-12: Objective 5 Completion - Local SpliceAI at -D 10000 Integration](./walkthrough/20260812_objective5_spliceai_d10000.md)
- [2026-08-12: Objective 6 Completion - Branch Point Predictor Pair (Branchpointer + LaBranchoR)](./walkthrough/20260812_objective6_branchpoint_predictors.md)
- [2026-08-12: Pipeline Pre-Run Readiness & Verification Audit](./walkthrough/20260812_pipeline_pre_run_readiness_audit.md)
- [2026-08-12: Apptainer / Singularity Container Architecture & Build Guide](./walkthrough/20260812_apptainer_container_architecture_guide.md)
- [2026-08-13: Pipeline Run Execution Audit & Diagnostics (run_20260812_1750)](./walkthrough/20260813_last_pipeline_execution_audit.md)
- [2026-08-17: HPC Cluster Run Progress & Throughput Diagnosis (run_20260813_1028)](./walkthrough/20260817_pipeline_run_progress_audit.md)
- [2026-08-18: Multi-Evidence Variant Prioritization & Tiering Guide](./walkthrough/20260818_variant_prioritization_and_tiering_guide.md)
- [2026-08-18: High-Throughput SpliceAI Parallel VCF Chunking](./walkthrough/20260818_spliceai_parallel_vcf_chunking.md)
- [2026-08-19: Pipeline Session Handover Summary & Architecture Updates](./walkthrough/20260819_pipeline_session_summary.md)
- [2026-08-19: Comprehensive Dual Pipeline Execution Progress Audit](./walkthrough/20260819_dual_pipeline_execution_progress_audit.md)
- [2026-08-19: Source Directory Reorganization & Obsolete Script Archival](./walkthrough/20260819_src_directory_reorganization_and_cleanup.md)
- [2026-08-20: Interactive Multi-Track Transcript & Exon Lollipop Visualization](./walkthrough/20260820_transcript_visualization_implementation.md)
- [2026-09-04: Master Pipeline Architecture Audit & Improvement Roadmap](./walkthrough/20260904_pipeline_analysis_and_improvements.md)
- [2026-09-10: Full Technical & Methodological Audit (root-cause of broken filter_variants.py import, disconnected filtering module, inert QC, ACMG/scoring issues)](./walkthrough/20260910_full_technical_methodological_audit.md)
- [2026-09-10: Predictor & Annotation Inventory + Audit (full list of every predictor actually running, incl. SpliceVarDB evidence silently inert, duplicate REVEL/SIFT/PolyPhen sources, dead MaxEntScan)](./walkthrough/20260910_predictor_inventory_and_audit.md)
- [2026-09-10: ACMG Criteria & Priority Tier — Deprioritization Note (real classification now happens in the `clinical_variant_prioritization` web app; this pipeline's ACMG/tiering logic is not yet deprecated but is no longer the primary decision surface)](./walkthrough/20260910_acmg_tiering_deprioritization_note.md)
- [2026-09-18: DMD VEP OOM Incident, Run-Provenance Gap, and Resource-Sizing Handoff (live incident root-caused to a fixed `h_vmem` not scaled to gene size — `ORCH-10`; dormant timing/memory logger found — `ORCH-11`; SpliceAI `-D` window literature review; session handoff with concrete next steps)](./walkthrough/20260918_dmd_vep_oom_incident_and_resource_provenance_handoff.md)
- [2026-09-30: Pangolin annotation db rebuilt from GENCODE 45 (`MISC-12`)](./walkthrough/20260930_pangolin_db_gencode_v45_rebuild.md)
- [2026-10-01: Predictor score attribution audit (wrong-gene / wrong-transcript scores on overlapping loci; alias-vs-neighbour resolution)](./walkthrough/20261001_predictor_gene_transcript_attribution_audit.md)
- [2026-10-02: Slow and hung Nextflow tasks in the S223 and panel7 runs (`MISC-16`, `MISC-17`)](./walkthrough/20261002_nextflow_slow_and_hung_tasks_s223_panel7.md)

---

## 🔍 Quality Control & Run Auditor

Run the built-in auditor script anytime to verify output `.pq` files, column schemas, annotation completeness, and log health across all genes:

```bash
python3 src/python/audit_run_results.py --run-dir RUNS/<RUN_NAME> --raw-dir RUNS/predictors_050826/input/by_gene
```

See [BUG_TRACKER.md](./BUG_TRACKER.md) for the current status of every known defect (fixed, in progress, or open), and [walkthrough/20260910_full_technical_methodological_audit.md](./walkthrough/20260910_full_technical_methodological_audit.md) for full technical/methodological detail on each.

---

## ✉️ Contact
@Mario Ruiz - Lab PGP

The standalone legacy filter and its CLI tests are archived. The active parser retains `--freq_col` and the population availability contract described above.
