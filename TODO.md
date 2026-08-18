# 📋 Master Annotation Pipeline & Clinical Dashboard Roadmap

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
- [ ] **Automated ACMG / AMP In-Silico Evidence Code Engine**
  - [ ] Rule engine for `PVS1` (null LoF in intolerant gene with `pLI >= 0.90` / `LOEUF < 0.35`)
  - [ ] Rule engine for `PM2_Supporting` (gnomAD v4 joint AF $< 0.0001$ or absent)
  - [ ] Rule engine for `PP3` (concordant in-silico: `SpliceAI > 0.50` + `Pangolin > 0.50` OR `AlphaMissense > 0.56` + `REVEL > 0.70`)
  - [ ] Rule engine for `BP4` (concordant benign predictions)
  - [ ] Add `ACMG_Evidence_Codes` and `ACMG_Suggested_Tier` columns to Parquet outputs and reports
- [ ] **Cardiomyopathy Domain & Cardiac Isoform / PSI Filter**
  - [ ] Add cardiac ventricle Percent Spliced In (PSI > 85%) flag for `TTN` truncating variants (TTNtv in A-band)
  - [ ] Map critical DCM structural domains (`BAG3` BAG domain, `LMNA` rod domain, `FLNC` Ig-like folds)
- [ ] **ClinGen Dosage Sensitivity & Gene Constraint Metrics**
  - [ ] Integrate ClinGen Haploinsufficiency (`HI`) and Triplosensitivity (`TS`) scores into KPI header banners
  - [ ] Display gnomAD v4 gene constraint metrics (`LOEUF`, `pLI`, missense $Z$-score)

---

## 📊 Phase 3: Reporting & Clinical Dashboard Enhancements
- [ ] **Multi-Gene Cohort Master Dashboard (`cohort_master_dashboard.html`)**
  - [ ] Aggregate top-tier prioritized variants across all 73 genes into a unified cross-gene dashboard
  - [ ] Cohort-level summary metrics (distribution of High/Moderate candidates, splicing vs missense burden)
  - [ ] Interactive gene selector and global variant search across the entire cohort
- [ ] **1-Click Exportable Filtered Candidates (.xlsx / .tsv)**
  - [ ] Add an in-browser "Download Filtered Candidates (.xlsx)" button to HTML reports
  - [ ] Auto-generate an executive Excel summary sheet alongside HTML reports
- [ ] **Visual Splice Junction Delta-Score Track**
  - [ ] Interactive lollipop / delta-score plot showing donor/acceptor gains and losses relative to exon-intron boundaries in the Splicing Tab

---

## 🔒 Phase 4: Containerization & Infrastructure Stability
- [ ] **Unified Apptainer / Singularity SIF Container (`annotation_suite.sif`)**
  - [ ] Package PyTorch, TensorFlow, Pangolin, SpliceAI, and VEP utilities into a single immutable `.sif` image
  - [ ] Eliminate all conda path, shebang, and Python ABI mismatches across heterogeneous HPC compute nodes
