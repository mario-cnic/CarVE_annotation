# Universal Genomic Variant Annotation, Filtering, and Visualization Pipeline

This repository contains a production-ready, modular HPC bash and Python pipeline ([main.sh](./main.sh)) for parsing, annotating, filtering, and visualizing human genomic variants. It accepts diverse input formats (cDNA ENST/NM notation, hg19 or GRCh38 genomic coordinates, direct VCFs, Excel, TSV, CSV, or Parquet), standardizes all variants to **GRCh38**, runs state-of-the-art functional and splicing predictor ensembles (**Ensembl VEP 111, CADD v1.6, REVEL, AlphaMissense, SPiPv2.1, Pangolin, Local SpliceAI at `-D 10000`, Branchpointer, LaBranchoR, and SpliceVault empirical RNA evidence**), applies customizable variant filtering strategies, and generates both publication-grade static plots (PDF/PNG) and interactive HTML dashboards.

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
   - **Pangolin Splicing Predictor**: PyTorch deep learning ensemble predicting splice site gain/loss in a 20kb window (`-d 10000`) with cardiac left ventricle (`heart_lv`) and atrial appendage (`heart_aa`) score deltas.
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

5. **Automated Visualization Suite**:
   - **Static Publication Plots (R ggplot2)**: Consequence breakdown barplots, REVEL vs AlphaMissense correlation scatter plots, and AF spectrum histograms saved as PDF and PNG.
   - **Interactive HTML Dashboards (Plotly)**: Standalone HTML report containing interactive charts with hover tooltips and filtering metric summaries.

---

## 📁 Repository Architecture & RUNS Folder

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
| `--raw-dir <DIR>` | Directory | **(Required)** Path to input variant directory containing `.xlsx`, `.csv`, `.tsv`, `.pq`, or `.vcf` files |
| `--run-name <NAME>`| String | Custom run folder name in `RUNS/` (default: folder name of `--raw-dir`) |
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

---

## 🔍 Quality Control & Run Auditor

Run the built-in auditor script anytime to verify output `.pq` files, column schemas, annotation completeness, and log health across all genes:

```bash
python3 src/python/audit_run_results.py --run-dir RUNS/<RUN_NAME> --raw-dir RUNS/predictors_050826/input/by_gene
```

---

## ✉️ Contact
@Mario Ruiz - Lab PGP
