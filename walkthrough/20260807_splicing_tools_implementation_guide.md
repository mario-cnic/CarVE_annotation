# Step-by-Step Implementation Guide for Splicing Predictors & Regulatory Annotation Tools

**Date**: 2026-08-07  
**Author**: Bioinformatics Pipeline Architecture Team  
**Based on Document**: [resources/splicing_predictor_evidence_review.html](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/splicing_predictor_evidence_review.html)  
**Target Reference Build**: GRCh38 (Ensembl VEP 111 / gnomAD v4.1)

---

## 1. Executive Summary & Design Architecture

This implementation guide translates the comprehensive evidence review ([resources/splicing_predictor_evidence_review.html](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/splicing_predictor_evidence_review.html)) into an actionable, prioritized engineering roadmap for the human genomic variant annotation pipeline ([main.sh](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/main.sh)).

### Key Architectural Decisions

1. **Prioritized Phased Rollout**: 7 core implementation priorities ordered by scientific leverage, license cleanliness, compute efficiency, and structural prerequisite dependencies.
2. **Three-Layer Hybrid Region Model**:
   - **Layer 1 (Analysis Substrate)**: Continuous geometric & sequence features (`intron_offset_signed`, `distance_to_donor`, `distance_to_acceptor`, `intron_length`, `splice_side`).
   - **Layer 2 (Functional Class)**: Biological mechanism (Canonical $\pm 1,2$, Donor Motif, Acceptor Motif, Branch Point / PPT, Deep Intronic $>50\text{ bp}$).
   - **Layer 3 (Reporting Skin)**: Derived 5-class side-aware grid used exclusively for presentation and ACMG/AMP evidence mapping (preventing metric pooling artifacts).
3. **Strict Status Contract (Eliminating "NULL-as-Zero")**:
   Every score column must carry explicit status semantics:
   - `scored`: Valid numerical prediction calculated.
   - `not_applicable`: Variant outside the tool's biological mechanism or strand (e.g. Branchpointer evaluated on a donor-side variant).
   - `not_covered`: Variant inside target class but outside tool's genomic window/receptive field.
   - `error`: Computational or data parsing exception.

---

## 2. Implementation Roadmap & Priority Matrix

| Priority | Tool / Feature | Type | Primary Mechanism / Advantage | License & Provenance | Compute & Hardware | Integration Path | Status / Blockers |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **P1** | **SPiP v2.1** | Predictor | Cryptic sites, branch point loss, exonic regulatory elements ($\Delta t$ ESRseq) | MIT (Clean) / GRCh38 Native | CPU (12 cores), Rscript | `main.sh` (Step 04) | 1 Blocker: SourceForge dataset download |
| **P2** | **SpliceVault** | Empirical RNA | Real transcriptome aberrant splicing events (exon skipping, cryptic junctions) | Open / Ensembl FTP | Tabix lookup (<1 sec) | VEP Plugin or Tabix Join | No blockers |
| **P3** | **Signed `intron_offset`** | Core Refactor | Signed distance (`+` donor, `-` acceptor) & side classification | Internal Pipeline Logic | Pure Python / Bash | `parse_variant_inputs.py` | Prerequisite GATE for P4–P6 |
| **P4** | **Pangolin** | Predictor | Splice site usage & native cardiac tissue channel (`Heart_Left_Ventricle`) | GPL-3.0 / GRCh38 Native | CPU (PyTorch), ~400 MB DB | Python CLI Step | Build DB locally from GENCODE GTF |
| **P5** | **SpliceAI `-D 10000`** | Predictor | Re-run SpliceAI locally to 10 kb receptive field (replaces 50 bp precomputed cap) | CC BY-NC 4.0 / Local | CPU (TensorFlow), GRCh38 FASTA | `annotation_vep.sh` / Python | No blockers (shared FASTA with Pangolin) |
| **P6** | **Branch Point Pair** | Predictor | Branchpointer (18–44 nt) + LaBranchoR (70 nt precomputed hg19$\rightarrow$GRCh38) | Open / Bioconductor & Precomputed ISM | R / `liftOver` | Acceptor-gated columns | Do NOT live-install LaBranchoR; use precomputed ISM |
| **P7** | **SpliceVarDB** | Validation | Experimental status lookup (minigene / RT-PCR literature database) | Open DB (`splicevardb.org`) | API / Tabix query | Validation Column | 1 Blocker: Domain allowlist approval |
| *Defer* | **MetaSplice** | Predictor | Multi-tier ensemble | GPU-dependent / hg19-native | High (Docker, 100GB+ DB) | External only | DEFERRED: GPU & Docker missing, hg38 divergence |
| *Defer* | **AbSplice Precomputed** | Predictor | Tissue-specific aberrant splicing | Restricted / Zenodo | Heavy Snakemake | Capture-side mask | DEFERRED: 100 bp component-zero flaw |
| *Watch* | **OpenSpliceAI** | Predictor | Open-source SpliceAI alternative | GPL-3.0 / GRCh38 Native | CPU (PyTorch) | Drop-in replacement | Ready as fallback if CC BY-NC 4.0 blocks usage |
| *Watch* | **Evo 2 / Borzoi** | Foundation | Genome foundation models | GPU-gated | High (24–80 GB VRAM) | Subagent skill | Gated on GPU server availability |

---

## 3. Detailed Tool Specifications & Step-by-Step Implementation

---

### Priority 1: SPiP v2.1 (Splicing Prediction in PIPelines)

#### 1. Description & Provenance
SPiP v2.1 is an R-based deep learning predictor that evaluates complex splicing alterations across entire genes. It is natively compiled for GRCh38 and carries an **MIT License**, making it the most publication-clean tool in the suite.

#### 2. Biological Mechanism & Gap Addressed
Addresses three critical gap mechanisms in a single execution pass:
- Cryptic / *de novo* splice site creation (`posCryptMut`, `classProbaCryptMut`).
- Branch point region disruption (`mutInPBarea` via BPP across -18 to -44 nt).
- Exonic Splicing Regulatory (ESR) sequence alterations (`deltaESRscore` via $\Delta t$ ESRseq).

#### 3. Environment & Resource Requirements
- **Conda Environment**: `/home/mruizp/data_lab_PGP/shared/utils/conda_envs/spip_env/`
- **Packages**: `Rscript`, `foreach`, `doParallel`, `randomForest`.
- **Compute**: 12 CPU cores (`nproc --all`, set `OMP_NUM_THREADS=8`). Execution time: ~369 variants in <3 minutes.

#### 4. Blocker Resolution
- **Blocker**: `transcriptome_hg38.RData` is hosted on SourceForge (not currently on local allowlist).
- **Resolution**: Request one-time domain allowlist approval for `sourceforge.net` or manually mirror `transcriptome_hg38.RData` to `/home/mruizp/data_lab_PGP/shared/utils/spip_resources/`.

#### 5. Integration Step in `main.sh`
```bash
# Step 04: SPiP Annotation
Rscript /home/mruizp/data_lab_PGP/shared/utils/SPiP/SPiP_v2.1.R \
  -i RUNS/${RUN_NAME}/_tmp/input_variants.vcf \
  -o RUNS/${RUN_NAME}/annotation/spip_output.tsv \
  -g GRCh38 \
  -t 12
```

#### 6. Column Contract & Status Semantics
- `SPiP_score`: Float [0.0 - 1.0] (`classProbaComplex`).
- `SPiP_recommendation`: String (`High Risk`, `Medium Risk`, `Low Risk`).
- `SPiP_mechanism`: String (`cryptic_donor`, `cryptic_acceptor`, `bp_disruption`, `esr_alteration`).
- **Status Semantics**: `'scored'`, `'not_covered'` (if variant is non-gene/intergenic), `'error'`.

#### 7. ACMG/AMP Threshold Calibration
- `SPiP_score >= 0.50`: PP3 / PM4 support for splicing disruption.
- `SPiP_score >= 0.85`: Strong PVS1_Moderate candidate when combined with SpliceAI/Pangolin.

---

### Priority 2: SpliceVault (Empirical RNA Aberrant Splicing Evidence)

#### 1. Description & Provenance
SpliceVault provides **empirical RNA-seq evidence** derived from >300,000 human RNA-seq samples (GTEx, SRA). It predicts exact aberrant event types (exon skipping, cryptic donor/acceptor usage) and their empirical frequencies.

#### 2. Biological Mechanism & Gap Addressed
Supplies direct empirical transcriptomic validation. Unlike sequence-only neural networks, it quantifies what actual human cells do when a nearby splice junction is disrupted. Fully supports indels.

#### 3. Environment & Resource Requirements
- **Dataset Size**: 843 MB (`SpliceVault_data_GRCh38.tsv.gz`, index 730 KB).
- **Location**: Ensembl FTP (`ftp.ensembl.org`, fully allowlisted).
- **Compute**: Zero model inference; instant tabix indexing lookup (<0.1s per variant).

#### 4. Integration Step in Pipeline
Integrate via official Ensembl VEP 111 plugin or custom python tabix join in `src/python/annotate_splicevault.py`:
```bash
python3 src/python/annotate_splicevault.py \
  --input-parquet RUNS/${RUN_NAME}/results/${GENE}.parsed.clean.pq \
  --db /home/mruizp/data_lab_PGP/shared/utils/splicevault/SpliceVault_data_GRCh38.tsv.gz \
  --output-parquet RUNS/${RUN_NAME}/results/${GENE}.parsed.clean.pq
```

#### 5. Column Contract & Status Semantics
- `SpliceVault_event_type`: String (`exon_skipping`, `cryptic_acceptor`, `cryptic_donor`, `intron_retention`).
- `SpliceVault_sample_count`: Integer (Count of RNA-seq samples witnessing event).
- `SpliceVault_fraction`: Float [0.0 - 1.0] (Relative proportion of aberrant read count).
- **Status Semantics**: `'no_events_found'` (covered in dataset, no aberrant splicing witnessed), `'event_detected'`, `'not_covered'`.

---

### Priority 3: Signed `intron_offset` & Splice Side Refactor

#### 1. Description & Architectural Significance
**Prerequisite structural refactor** for the core variant parsing module (`src/python/parse_variant_inputs.py`). Currently, `intron_offset` is recorded as an absolute positive integer, making `c.123+5` and `c.123-5` indistinguishable (`5`), conflating donor-side and acceptor-side variants.

#### 2. Biological Mechanism & Gap Addressed
Features like Branch Point (BP) and Polypyrimidine Tract (PPT) are strictly **acceptor-side-only** biological structures (-10 to -50 nt upstream of exon start). Without a signed offset, acceptor-specific tools process donor-side variants, creating false positives and breaking `not_applicable` state enforcement.

#### 3. Engineering Implementation in Python
Modify `src/python/parse_variant_inputs.py`:
```python
def compute_signed_intron_offset(hgvs_c_str: str) -> tuple[int, str]:
    """
    Parses HGVS cDNA string (e.g. c.123+5C>T, c.456-12A>G)
    Returns: (signed_offset, splice_side)
    """
    import re
    match = re.search(r'c\.\d+([+-])(\d+)', hgvs_c_str)
    if not match:
        return (0, 'exonic')
    sign, val = match.groups()
    offset = int(val) if sign == '+' else -int(val)
    side = 'donor' if sign == '+' else 'acceptor'
    return (offset, side)
```

#### 4. Column Contract
- `intron_offset_signed`: Integer (Positive for donor/downstream intron, negative for acceptor/upstream intron, 0 for exonic).
- `splice_side`: Enum (`donor`, `acceptor`, `exonic`).

---

### Priority 4: Pangolin Splicing Predictor

#### 1. Description & Provenance
Pangolin is a deep learning model (GPL-3.0 licensed) trained on tissue-specific splice site **usage** rather than just site recognition. It provides native tissue-specific loss/gain scores for cardiac tissues (`Heart_Left_Ventricle`, `Heart_Atrial_Appendage`).

#### 2. Biological Mechanism & Gap Addressed
Evaluates splice site strength changes in cardiac-relevant tissue context. Offers an unbounded distance parameter (`-d`), rendering deep intronic regions accessible.

#### 3. Local Annotation DB Construction (Bypassing Dropbox Blocker)
- **Do NOT download prebuilt databases from Dropbox** (violates network allowlist rules).
- Build the database locally using the GENCODE v44 GTF:
```bash
python3 /path/to/pangolin/create_db.py \
  --gtf /home/mruizp/data_references/GRCh38/gencode.v44.annotation.gtf \
  --fasta /home/mruizp/data_references/GRCh38/GRCh38.primary_assembly.genome.fa \
  --output /home/mruizp/data_lab_PGP/shared/utils/pangolin_db/pangolin_grch38.db
```

#### 4. Execution Step
```bash
pangolin \
  RUNS/${RUN_NAME}/_tmp/input_variants.vcf \
  /home/mruizp/data_references/GRCh38/GRCh38.primary_assembly.genome.fa \
  /home/mruizp/data_lab_PGP/shared/utils/pangolin_db/pangolin_grch38.db \
  RUNS/${RUN_NAME}/annotation/pangolin_output.vcf \
  -d 10000 \
  --distance-cutoff 10000
```

#### 5. Column Contract & Status Semantics
- `Pangolin_max_score`: Float [0.0 - 1.0].
- `Pangolin_heart_lv_score`: Float [0.0 - 1.0] (Heart Left Ventricle specific delta).
- `Pangolin_heart_aa_score`: Float [0.0 - 1.0] (Heart Atrial Appendage specific delta).
- **Status Semantics**: `'scored'`, `'not_covered'` (distance > 10,000 nt), `'error'`.

---

### Priority 5: Local SpliceAI Re-run at `-D 10000`

#### 1. Description & Observability Gap Resolution
Public precomputed SpliceAI VCFs (e.g. Illumina 2019 precomputed tables) truncate predictions at **50 bp** from exon boundaries. However, SpliceAI's underlying neural network receptive field is **10,000 nt**. Currently, all deep intronic variants (>50 bp) read from precomputed tables are unscored-by-construction (`NULL`), not negative.

#### 2. Execution Protocol
Re-run SpliceAI locally with parameter `-D 10000` using the local PyTorch/TensorFlow environment:
```bash
spliceai \
  -I RUNS/${RUN_NAME}/_tmp/input_variants.vcf \
  -O RUNS/${RUN_NAME}/annotation/spliceai_10k.vcf \
  -R /home/mruizp/data_references/GRCh38/GRCh38.primary_assembly.genome.fa \
  -A grch38 \
  -D 10000
```

#### 3. Column Contract & Provenance Enhancement
Extend existing pipeline SpliceAI schema with provenance fields:
- `SpliceAI_DS_max`: Float [0.0 - 1.0] (Max delta score among AG, AL, DG, DL).
- `SpliceAI_DP_max`: Integer [-10000 to +10000] (Delta position).
- `SpliceAI_distance_provenance`: String (`computed_10k`, `precomputed_50bp`).
- **Status Semantics**: `'scored'`, `'not_covered'` (distance > 10,000 nt).

---

### Priority 6: Branch Point Predictor Pair (Branchpointer + LaBranchoR Precomputed)

#### 1. Biological Justification & Pair Rationale
The Branch Point (BP) region is not a fixed 15 bp box. Lariat RNA-seq studies demonstrate that while 90% of branch points lie in the 19–37 nt window (Mercer 2015), ~17% extend from 45 to 70 nt upstream (Taggart 2017, Pineda & Bradley 2018).
- **Branchpointer**: Focuses on the core 18–44 nt acceptor window.
- **LaBranchoR**: Extends reach up to 70 nt upstream.

#### 2. Execution & Precomputed Data Integration
- **Branchpointer**: Installed via Bioconductor (`biocLite("branchpointer")`) and executed via R script on acceptor-side variants (`splice_side == 'acceptor'`).
- **LaBranchoR**: **DO NOT attempt live installation** (GitHub repo unmaintained since Feb 2019 on obsolete Keras stack). Ingest the **precomputed ISM tables** (70 nt upstream of GENCODE exons) and perform an hg19$\rightarrow$GRCh38 liftover via `CrossMap` or `liftOver`.

#### 3. Column Contract & Acceptor Gating
- `Branchpointer_prob`: Float [0.0 - 1.0].
- `LaBranchoR_score`: Float [0.0 - 1.0].
- **Status Semantics**:
  - `splice_side == 'acceptor'`: `'scored'` or `'not_covered'`.
  - `splice_side == 'donor'` or `'exonic'`: `'not_applicable'`.

---

### Priority 7: SpliceVarDB Experimental Validation Lookup

#### 1. Description & Provenance
SpliceVarDB is a curated repository of experimentally tested genetic variants with documented splicing outcomes from published minigene assays and patient RT-PCR studies.

#### 2. Blocker Resolution & Integration
- **Blocker**: Requires domain allowlist approval for `splicevardb.org` / `compbio.ccia.org.au`.
- **Query Script**: Perform automated tabix or REST API join against SpliceVarDB release 2026.

#### 3. Column Contract & Benchmarking Warning
- `SpliceVarDB_status`: String (`tested_splice_altering`, `tested_normal`, `not_in_database`).
- `SpliceVarDB_assay`: String (`minigene`, `rt_pcr`, `rna_seq`).
- `SpliceVarDB_pmid`: String (PubMed Identifier).
- **CRITICAL CIRCULARITY WARNING**: SpliceVarDB is frequently used as a benchmark gold-standard set. **Do NOT use SpliceVarDB entries to benchmark or train local predictors.**

---

## 4. Derived 5-Class Side-Aware Region Grid & ACMG/AMP Mapping

To report results without introducing metric pooling artifacts, the pipeline implements a derived reporting layer (Layer 3) mapped directly from Layer 1 continuous geometry:

```mermaid
flowchart TD
    V[Variant Coordinate] --> G{intron_offset_signed}
    G -->|offset == 0| R1[1. Exonic Region]
    G -->|+1, +2 or -1, -2| R2["2. Canonical Dinucleotides (±1,2)"]
    G -->|+3 to +6 or -3 exonic| R3[3. Donor Motif Region]
    G -->|-3 to -14| R4[4. Acceptor Proximal / PPT Region]
    G -->|-15 to -44| R5[5. Branch Point Region]
    G -->|> |50| bp| R6[6. Deep Intronic Region]
```

### ACMG/AMP Splicing Evidence Calibration Rules

| Region Class | Distance / Bounds | Primary Predictors | Moderate / Strong Cutoff | ACMG/AMP Code |
| --- | --- | --- | --- | --- |
| **Canonical Dinucleotide** | $\pm 1, 2$ | VEP Consequence / All | Null allele in LOF gene | **PVS1** (Very Strong) |
| **Donor / Acceptor Motif** | $+3..+6$ / $-3..-14$ | SpliceAI $\ge 0.50$, SPiP $\ge 0.50$, Pangolin $\ge 0.50$ | 2+ Predictors Concordant | **PP3** / **PVS1_Moderate** |
| **Branch Point** | $-15..-44$ nt (Acceptor) | Branchpointer $\ge 0.50$, SPiP `mutInPB` | 2+ Predictors Concordant | **PP3** / **PM4** |
| **Deep Intronic** | $> 50$ bp | SpliceAI (10k) $\ge 0.20$, SPiP $\ge 0.50$, Pangolin $\ge 0.20$ | Concordant Cryptic Call | **PP3** / **PVS1_Optional** |

---

## 5. Verification & Testing Plan

### Automated Pipeline Audit & Schema Verification
Execute the audit engine to verify column presence, non-null percentages, and status semantics across all priority tools:
```bash
python3 src/python/audit_run_results.py \
  --run-dir RUNS/predictors_050826 \
  --raw-dir RUNS/predictors_050826/input/by_gene
```

### Expected Output Criteria
1. `intron_offset_signed` and `splice_side` present with 100% completeness.
2. Acceptor-only tools (`Branchpointer`, `LaBranchoR`) display status `'not_applicable'` for all `donor` and `exonic` variants.
3. SpliceAI distance provenance correctly distinguishes between `precomputed_50bp` and `computed_10k`.
4. Zero unhandled `NULL` values in numerical score columns.

---

## 6. Document History & Provenance

- **2026-08-07**: Implementation Guide initialized from comprehensive evidence review ([resources/splicing_predictor_evidence_review.html](file:///home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/resources/splicing_predictor_evidence_review.html)).
- **Maintainer**: Bioinformatics Pipeline Development Team.
