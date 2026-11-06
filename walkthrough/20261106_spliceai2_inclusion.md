# SpliceAI2: inclusion study for Module 2

2026-11-06. Kanvas PR-02. Sources read on this date: the GitHub README (`Illumina/SpliceAI2`), the Hugging Face dataset page (`illumina-ai/SpliceAI2-data`) and the bioRxiv abstract (doi 10.64898/2026.10.08.757673). The gated files were not downloaded and the full preprint text was not read.

## What it is (FACT, from the sources)

- Illumina's successor to SpliceAI (Jaganathan et al., preprint posted 2026-10-11, not peer reviewed). It predicts splice-site usage, junctions and full transcripts from sequence and gives variant effects quantitatively.
- Benchmarks (cryptic splicing, sQTL, MPRA) and the rare-disease result (Genomics England; 15% of the excess burden explained, about half deep intronic) were all run by the developers. No independent evaluation was found.
- **Live scoring:** `pip install spliceai2` (PyTorch). The README says a CUDA GPU is required. The input is a TSV (`chrom pos ref alt strand`), not a VCF. The output is a TSV with `spliceai2_summary_score` plus 260 columns; 8 of them map to SpliceAI's `DS_*`/`DP_*`. The README shows no window/distance option.
- **Precomputed scores:** GRCh38, GENCODE v48, protein-coding genes. They cover every possible SNV (~4 billion) and population-observed indels (~150 million). The dataset is ~350 GB and gated on Hugging Face (accept terms); it is also served through DRAGEN Annotation.
- **Licence:** "SpliceAI2 Model Terms of Use". It covers academic and non-commercial research only, forbids redistribution and third-party access without Illumina's written consent, and requires attribution. SpliceAI 1.3.1, which we run now, is GPLv3.
- **Thresholds suggested by the authors** (SpliceAI2 / SpliceAI equivalent): 0.1 / 0.2 high recall, 0.25 / 0.5 balanced, 0.5 / 0.8 high precision. They are not ClinGen-calibrated; the ClinGen SVI splicing calibration (PP3/BP4) is for SpliceAI.

## What it means here (INFERENCE)

- **It does not replace SpliceAI for ACMG evidence.** PP3/BP4 thresholds exist only for SpliceAI, so SpliceAI 1.3.1 stays the scored predictor; SpliceAI2 would be informational, like ReMM and AVI.
- **The thesis fit is good on paper.** The model was built for deep-intronic and quantitative effects, which is the non-coding splice-hotspot question (MYBPC3 is protein-coding, so it is covered). Developer-only benchmarks mean it should be one more predictor in the comparison, not the reference.
- **Live scoring needs GPU nodes.** Our SpliceAI runs on CPU (TensorFlow). Whether the SGE cluster has CUDA nodes and a queue for them is unknown here.
- **Transcript set mismatch.** The precomputed scores use GENCODE v48, while VEP, Pangolin and the SpliceAI symbol map use GENCODE 45 / Ensembl 111. Gene attribution by symbol or ENSG works for almost all genes, but transcript-level fields would not line up row for row.
- **Licence questions for you.** (1) Is clinical reporting to families covered by "research"? (2) Redistribution: Module 3 tables shared outside the lab (collaborators, the thesis repository) would contain SpliceAI2 values. Both need an answer before the scores enter tables that leave the lab.

## Options (PROPOSAL)

| | What | Cost | Covers |
|---|---|---|---|
| A | Precomputed lookup, like ReMM/AVI: download only `precomputed_scores_v2.0/` GRCh38, convert to one bgzipped, tabix-indexed file, annotate in `VEP_ANNOTATE` (`--custom`) or in a small lookup step; add `SpliceAI2_score`, `SpliceAI2_DS_*`/`DP_*`, `SpliceAI2_status` (`scored` / `not_covered`) | download and disk (~hundreds of GB, size of the scores folder unknown until the gate is accepted); no GPU | all SNVs and population indels in protein-coding genes; private indels and non-coding genes stay `not_covered` |
| B | A + live scoring of the misses on GPU nodes, cached like SpliceAI/Pangolin (`annotation_cache.py`) | GPU queue, PyTorch env, TSV wrapper with strand from the gene | everything A misses |
| C | Do not integrate; score selected variants ad hoc for the thesis comparison | none in the pipeline | only what is asked |

**Recommendation:** A, informational, with `SpliceAI2_status` keeping "not covered" separate from 0 (as SpliceAI/ReMM/AVI already do). B only if the cluster has GPUs and private indels matter. Before any code:

1. You accept the Hugging Face terms (personal account) and answer the two licence questions above.
2. Check the size of `precomputed_scores_v2.0/` and its format (the dataset page does not say), then choose where it lives under `/data_lab_PGP/resources/annotation/`.
3. Confirm the build of the files (the README says GRCh38) by checking reference alleles against our FASTA on a sample, as for ReMM.

## Sources

- GitHub README: https://github.com/Illumina/SpliceAI2
- Dataset: https://huggingface.co/datasets/illumina-ai/SpliceAI2-data
- Preprint: https://www.biorxiv.org/content/10.64898/2026.10.08.757673v1
- Illumina press release: https://www.illumina.com/company/news-center/press-releases/press-release-details.html?newsid=29236
