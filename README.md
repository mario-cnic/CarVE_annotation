# Universal Genomic Variant Annotation, Filtering, and Visualization Pipeline

This repository contains a production-ready, modular HPC bash and Python pipeline ([main.sh](./main.sh)) for parsing, annotating, filtering, and visualizing human genomic variants. It accepts diverse input formats (cDNA ENST/NM notation, hg19 or GRCh38 genomic coordinates, direct VCFs, Excel, TSV, CSV, or Parquet), standardizes all variants to **GRCh38**, runs state-of-the-art functional and splicing predictor ensembles (**Ensembl VEP 111, CADD v1.6, REVEL, AlphaMissense, SPiPv2.1, Pangolin, Local SpliceAI at `-D 10000`, Branchpointer, LaBranchoR, and SpliceVault empirical RNA evidence**), applies customizable variant filtering strategies, and generates both publication-grade static plots (PDF/PNG) and interactive HTML dashboards.

A second entry point, the Nextflow whole-VCF pipeline ([annotate_vcf.nf](./annotate_vcf.nf)), annotates one arbitrary-size GRCh38 VCF (panel, WES or WGS) in a single run. See [Nextflow Whole-VCF Pipeline](#-nextflow-whole-vcf-pipeline) below. It does not replace `main.sh`, which is unchanged.

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
   - **SPiP v2.1**: Multi-threaded empirical and machine learning splicing impact predictor (up to 12 cores) with decoded mechanisms (`Exon_skipping`, `Donor_disruption`, etc.).
   - **Pangolin Splicing Predictor**: PyTorch deep learning ensemble predicting splice site gain/loss in a 20kb window (`-d 10000`), with GENCODE 45 Ensembl_canonical gene models. Output produced before 2026-09-30 is invalid: see the warning above (`MISC-12`). Pangolin has one generic heart model, not separate left-ventricle/atrial-appendage models, so the downstream `heart_lv`/`heart_aa` columns are not subtissue-specific (`PRED-1`, `PLT-138`).
   - **Local SpliceAI at `-D 10000`**: Deep residual neural network predicting max delta scores across a 20kb intronic window.
   - **Branch Point Predictor Pair**:
     - **LaBranchoR**: Bidirectional LSTM predicting top human catalytic branch points across GRCh38 (206,249 branch points).
     - **Branchpointer**: Signal-feature model estimating branch point probability, U2 snRNA duplex free energy (kcal/mol), and motif disruption flags (`YES`/`NO`).
   - **SpliceVault**: Direct empirical RNA-seq aberrant event lookup across GTEx/SRA datasets with human-readable decoded events and sample counts.
   - **SpliceVarDB**: Direct lookup of experimentally validated splicing variants from minigene, RT-PCR, and RNA-seq assays.

3. **Unified Status Contracts & Data Schema**:
   - All predictors adhere to strict, unified status contracts: `scored` (valid numerical prediction), `not_covered` (outside gene model/receptive field), or `error`.
   - For complete column documentation, see [resources/column_dictionary.md](./resources/column_dictionary.md).

4. **Customizable Filtering Module**:
   - Customizable command-line thresholds for gnomAD AF (`--max-af`), REVEL (`--min-revel`), AlphaMissense (`--min-alphamissense`), SPiP (`--min-spip`), CADD (`--min-cadd`), and VEP consequences (`--consequences`).

5. **Automated Visualization & Clinical Report Suite**:
   - **Static Publication Plots (R ggplot2)**: Consequence breakdown barplots, REVEL vs AlphaMissense correlation scatter plots, and AF spectrum histograms saved as PDF and PNG.
   - **Interactive 4-Tab Clinical Prioritization Dashboards (`generate_clinical_prioritization_report.py`)**:
     - **Multi-Track SVG Transcript Visualizer & Variant Position Map**: Interactive Plotly multi-track diagram mapping genomic/cDNA positions and amino acid coordinates across exon-intron boundaries, displaying priority tiers, SpliceAI Δ, and AlphaMissense scores.
     - **gnomAD v4.1 Joint PopMax Frequencies**: Full extraction of `gnomADv4_AF_grpmax_joint` with population ancestry labels (`AFR`, `AMR`, `EAS`, `NFE`, `SAS`, etc.).

6. **Standalone Clinical Reporting Web Application ([run_web_app.sh](./run_web_app.sh))**:
   - Drag-and-drop web UI built with Streamlit for clinicians and researchers.
   - Upload any variant file (`.vcf`, `.vcf.gz`, `.parquet`, `.pq`, `.xlsx`, `.csv`, `.tsv`).
   - Automated ACMG/ClinGen priority tiering, live interactive filtering, Plotly charts, embedded 4-tab clinical report preview, and export center (`.html`, `.xlsx`, `.pq`).

---

## 📁 Repository Architecture & Codebase Organization

### Execution Directory (`RUNS/<RUN_NAME>/`)
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

### Source Modules (`src/`)
The codebase is organized into clean, functional modules:

```
src/
├── python/                                   # Core active Python pipeline modules (Step 1-6)
├── hpc/                                      # Core SGE cluster drivers & execution wrappers
├── R/                                        # Core plotting & statistical libraries
├── downstream/                               # Cohort burden & association analysis scripts
├── tools/                                    # Reference builders, backfills, & exploratory utilities
└── external/                                 # Third-party source dependencies (Pangolin-main)

config/
└── env.sh                                    # Centralized shell environment & binary auto-detection

archive/
└── deprecated_scripts/                       # Safely archived legacy & superseded scripts
```

---

## 💻 Usage & CLI Parameters

### Basic Run
Process all variant tables inside your target directory. Outputs are stored in `RUNS/<RUN_NAME>/`:

```bash
bash main.sh --raw-dir /data_lab_PGP/raw_data/my_cohort/
```

### Custom Timestamped Run (Preserving Past Runs)
Specify a custom, isolated run folder inside `RUNS/` (e.g. `RUNS/run_20260813_1100/`):

```bash
bash main.sh \
  --raw-dir RUNS/predictors_050826/input/by_gene \
  --run-name "run_$(date +%Y%m%d_%H%M)" \
  --build GRCh38 \
  --output-format pq
```

### Advanced Usage with Filtering & Build Options

```bash
bash main.sh \
  --raw-dir RUNS/predictors_050826/input/by_gene \
  --run-name cohort_valencia_hg19 \
  --build hg19 \
  --gene SCN5A \
  --output-format pq \
  --max-af 0.001 \
  --min-revel 0.75 \
  --min-alphamissense 0.80 \
  --min-spip 0.50 \
  --consequences "missense_variant,stop_gained,frameshift_variant,splice_acceptor_variant,splice_donor_variant"
```

### CLI Parameters:

| Flag | Type | Description |
| --- | --- | --- |
| `--raw-dir <DIR>` | Directory | **(Required for custom runs)** Path to input variant directory containing `.xlsx`, `.csv`, `.tsv`, `.pq`, or `.vcf` files |
| `--run-name <NAME>`| String | Custom run folder name in `RUNS/` (default: folder name of `--raw-dir`) |
| `--test` | Flag | Run end-to-end test suite on 50-variant test dataset into `test_data/test_run/` with `--overwrite-all` |
| `--audit-run <DIR>`| String | Run master auditor to inspect quality, completeness %, and schema consistency of a run folder |
| `--gene <GENE>` | String | Filter execution to specific gene (case-insensitive) |
| `--build <BUILD>` | String | Input assembly build: `hg19` or `GRCh38` (default: `GRCh38`) |
| `--output-format <FMT>` | String | Output format: `pq` (Parquet), `tsv`, or `xlsx` |
| `--max-af <FLOAT>` | Float | Maximum gnomAD allele frequency threshold |
| `--min-revel <FLOAT>` | Float | Minimum REVEL pathogenicity score |
| `--min-alphamissense <FLOAT>` | Float | Minimum AlphaMissense score |
| `--min-spip <FLOAT>` | Float | Minimum SPiP splicing score |
| `--min-cadd <FLOAT>` | Float | Minimum CADD phred score |
| `--consequences <LIST>` | String | Comma-separated list of target VEP consequences |
| `--overwrite-all` | Flag | Overwrite intermediate files and re-run all steps |

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
- **Finding Nextflow:** the wrapper uses `nextflow` from `PATH`, falling back to `/opt/nextflow/nextflow`.

| Argument | Default | Description |
| --- | --- | --- |
| `--input_vcf <FILE>` | required | bgzipped VCF with a sibling `.tbi` index. **Must be GRCh38**, called against the GATK-bundle `Homo_sapiens_assembly38.fasta` that the predictors use. A pre-flight check fails the run before any predictor starts if the `##contig` lengths or the first 1000 REF alleles disagree with that FASTA. Without `##contig` lines it continues with a warning only if ≥20 REF alleles match (`MISC-14`). Both `chr1` and `1` contig names are accepted |
| `--run_id <LABEL>` | input filename without `.vcf.gz` | Output folder name and prefix for every published file |
| `--spip_tier` | `restricted` | Which tier SPiP runs in (see below) |
| `--output_format` | `pq` | Final table format: `pq` or `tsv` only |
| `-profile` | `standard` | `standard` = SGE + Singularity (production). `local_dev` = local miniforge environments, for `-stub-run` testing only |
| `-resume` | off | Reuse cached tasks from Nextflow's `work/` directory |

### Two predictor tiers

| Tier | Records | Tools |
| --- | --- | --- |
| Broad | every record in the input | VEP 111 (with its plugins), Branchpointer/LaBranchoR |
| Gene-restricted | only records inside [`resources/v5_genes_loc.bed`](./resources/v5_genes_loc.bed) (`Complete_gene_list_V5` regions) | SpliceAI (`-D 10000`), Pangolin (`-d 10000`), SPiP (unless `--spip_tier broad`) |

Records outside the gene BED have **no** SpliceAI/Pangolin/SPiP annotation. A missing score there means the tool was not run, not that the variant is predicted benign.

### Outputs

Everything is published to `nf_work/annotation_out/<run_id>/`:

| File | Content |
| --- | --- |
| `<run_id>.assembly_check.tsv` | Pre-flight assembly check result: contig-length and REF-allele comparison against `params.fasta`, plus any warnings |
| `RUN_MANIFEST.json` | Launch/completion provenance: command line, git commit and dirty-tree status, resolved binaries, input/resource hashes, tool versions, exit code |
| `<run_id>.annVEP.vcf.gz`, `.annBranchpoint.vcf.gz`, `.annSpliceAI.vcf.gz`, `.annPangolin.vcf.gz`, `.annSPiP.vcf.gz` | Per-predictor annotated VCFs (+ `.tbi`) |
| `<run_id>.annotated.vcf.gz` | All predictors merged into one VCF (+ `.tbi`) |
| `<run_id>.parsed.clean.pq` | **Final table.** One row per variant × VEP transcript. Every transcript is kept and tagged with `TRANSCRIPT_PRIORITY_TIER`: 1 = curated 217-panel transcript, 2 = MANE Select, 3 = other |

For scale: on `S223` the final table was 5.0 GB and the merged VCF 2.7 GB. Intermediate task files live in Nextflow's `work/` directory. On a WGS run this directory holds days of compute; deleting it forces a full recompute.

The wrapper stages `RUN_MANIFEST.json` with `git add -f` but does not commit it. The manifest stays inside the run folder, so it moves with the run if the folder is relocated.

### Known limitations

- **Pangolin output before 2026-09-30 is invalid:** the old reference database was a mouse annotation (`MISC-12` / `PLT-139`). The GENCODE 45 rebuild must be deployed to `shared/utils/pangolin_db/` before the Pangolin step can run, and it is not yet validated on the cluster.
- **Genes missing from the gene BED:** 144 curated genes, including panel gene `TAZ`, get no gene-restricted predictors (`MISC-4`).
- **No cross-predictor record-count checks:** counts aren't compared between predictors or against the merge output (`MISC-6`).
- **Incomplete provenance for dirty launches:** a launch from a dirty working tree records only the names of changed files, so the exact code can't be reconstructed (`MISC-13`). Launch from a clean, committed tree.
- **Assembly check samples, it doesn't exhaustively verify:** the pre-flight check compares header contig lengths and only the first 1000 REF alleles (`MISC-14`). VEP runs without `--check_ref`, so per-record REF disagreements beyond that sample are not detected.
- **Unprofiled resource requests:** CPU/memory/time are not yet sized from profiling data. Per-chunk wall-clock varied from minutes to more than 46 h on `S223` (`TODO.md`, Priority 0).

`main.nf` is the older per-gene Nextflow adapter. It is kept working but is no longer developed.

---

## 📊 Pipeline Outputs & Data Schema

All output tables contain standardized, clean columns. For the full column dictionary, renamed field mappings, and status contracts (`scored`, `not_covered`, `error`), see [resources/column_dictionary.md](./resources/column_dictionary.md).

Output files are organized into structured subdirectories under `results/`:

- **Parsed Clean Annotations**: `results/<GENE>.parsed.clean.pq` (or `.tsv`, `.xlsx`)
- **Filtered Variant Tables**: `filtered/filtered_variants.pq`
- **Filtering Metrics Summary**: `filtered/filtering_summary_metrics.tsv`
- **Static Publication Plots**: `plots/consequence_distribution.pdf`, `pathogenicity_correlation.pdf`, `allele_frequency_spectrum.pdf`
- **Interactive HTML Dashboard**: `reports/<GENE>_interactive_dashboard.html`

---

## 📜 HPC Environment & Container Options

### Option 1: Apptainer / Singularity SIF Container (Recommended for HPC)
Eliminates NFS I/O latency and provides a single immutable execution image:

- **Container Definition**: [`resources/containers/annotation_pipeline.def`](./resources/containers/annotation_pipeline.def)
- **Build Container Image**:
  ```bash
  qsub src/hpc/build_container.sh
  ```
- **Architecture Guide**: [walkthrough/20260812_apptainer_container_architecture_guide.md](./walkthrough/20260812_apptainer_container_architecture_guide.md)

### Option 2: Conda Environments (`/data_lab_PGP/shared/utils/conda_envs/`)
- **Data Science & Post-Processing**: `datasci/`
- **Genomics & VCF Tools**: `genomics/`
- **Coordinate Conversion & Liftover**: `liftover/`
- **SPiP v2.1 Predictor**: `spip_env/`
- **Pangolin Predictor**: `pangolin_env/` (`resources/conda/pangolin_env.yml`)
- **SpliceAI Predictor**: `spliceai_env/` (`resources/conda/spliceai_env.yml`)

---

## 🧪 Master Unit Testing Suite

The repository includes a comprehensive unit testing suite covering Python, R, and Bash CLI components:

- **Run all unit tests**:
  ```bash
  bash tests/run_all_tests.sh
  ```

- **Run Python unit tests directly (`pytest`)**:
  ```bash
  python3 -m pytest -v tests/python/
  ```

- **Run R unit tests directly (`Rscript`)**:
  ```bash
  Rscript tests/R/test_statsJPO.R
  Rscript tests/R/test_plot_annotation_results.R
  ```

- **Run Bash unit tests directly (`sh`)**:
  ```bash
  bash tests/bash/test_logger.sh
  bash tests/bash/test_cli_args_main.sh
  ```

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
