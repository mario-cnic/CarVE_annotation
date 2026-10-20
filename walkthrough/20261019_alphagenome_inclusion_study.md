# AlphaGenome / AlphaGenome Atlas: inclusion study for Module 2

2026-10-19T14:54:56+02:00. Study only: nothing downloaded, no pipeline file changed. One local VEP 111 run on synthetic data (below).

## Decision proposed

**PROPOSAL.** Add the precomputed **AVI SNV score** file to `VEP_ANNOTATE` with Ensembl's `AVI.pm` plugin. This follows the AlphaMissense/ReMM pattern: a local tabix lookup, informational only. Do **not** call the AlphaGenome API or the Atlas API from the pipeline.

Mario decides on two points:

1. Download the 88.5 GB file. The site requires a Google sign-in, so the download has to be done from a browser.
2. Whether to also take the two non-commercial files, AlphaGenome merged splicing scores and the AVI SHAP feature importances. Either could serve the thesis hotspot map; neither belongs in the default run.

## What changed since the 2026-09 plan in Module 3

`../clinical_variant_prioritization/docs/ALPHAGENOME_INTEGRATION_PLAN.md` (§2, §3, §4.4) says the Atlas has "no separate bulk/offline export". **That is no longer true.** The rest of that plan, the on-demand API panel in the app, still stands as a Module 3 matter. The plan doc is not edited from this repo.

## Facts (checked 2026-10-19)

**Downloads page** (`https://deepmind.google.com/science/alphagenome/downloads`, read in a signed-in browser):

| File | Licence section on the page | Size |
|---|---|---|
| `avi_scores_snvs_tabix.zip`: "AVI scores and Phred-scaled scores in Tabix format" | Permissive use (commercial and non-commercial) | 88.5 GB |
| `combined_splicing_snvs_tabix.zip`: "AlphaGenome merged splicing scores in Tabix format" | Non-commercial only | 20.6 GB |
| `avi_feature_importances_snvs_tabix.zip`: "AVI SHAP feature importance scores in Tabix format" | Non-commercial only | 283.9 GB |

- Motif datasets: "will be made available soon".
- An unauthenticated `curl -I` on the download links returns HTTP 500. INFERENCE: a Google sign-in is required.
- Not yet seen, because nothing was downloaded: the column layout, the chromosome naming, the version and the zip contents.

**Terms** (AlphaGenome Services Additional Terms, last modified 2026-09-08, PDF `storage.googleapis.com/alphagenome/terms/20260908/...`):

- Use is limited to individuals and non-commercial organisations. Universities and research institutes qualify, so CNIC does.
- The body of the terms exempts AVI (as a "Permissive Use Downloadable Artifact") from the bans on commercial use and ML training, and from the Output-Terms notice.
- The "Key things to know" summary box also lists "clinical purposes" among the bans AVI is exempt from. **The binding body does not.**
  - "AlphaGenome Assets", which explicitly include AVI Scores, "are for theoretical modelling only. They are not intended, validated, or approved for clinical use."
  - Prohibition 6e, with no AVI carve-out, forbids use for "making decisions in domains that affect material or individual rights or well-being (e.g., healthcare)".
  - Read: AVI is permissive for commercial use and training, and still non-clinical. Module 3 must keep a research-only disclaimer on any AVI display.
- The "AVI Score Feature Breakdown" is not part of the AVI Score. So the SHAP file stays under all restrictions.
- The services must not be used "to transmit, store, or process health information subject to United States HIPAA regulations".

**What AVI is** (DeepMind blog, 2026-09-08):

- A single score per SNV that combines AlphaGenome and **AlphaMissense** predictions, over about 9 billion SNVs, coding and non-coding.
- No indels.
- The blog gives no formula and no calibration.

**`AVI.pm`** (`github.com/Ensembl/VEP_plugins`, last commit 2026-09-08, Apache-2.0):

- Reads a bgzipped, tabix-indexed TSV with columns `chr pos ref alt raw phred`.
- Writes `AVI_PHRED` (0 to about 80; 20 = top 1% genome-wide) and, optionally, `AVI_RAW` (about -1.3 to +4.0) and `AVI_LINK` (an Atlas URL).
- Returns nothing for non-SNVs.
- It is not in `/data_lab_PGP/resources/annotation/VEP/VEP_plugins/`. AlphaMissense.pm and PromoterAI.pm are.

**VEP 111 compatibility, tested locally:**

- `vep.sif` reports `ensembl-vep : 111.0`.
- The plugin needs `get_user_params` and `params_to_hash`. Both exist in the container, in `BaseVepTabixPlugin.pm:163` and `BaseVepPlugin.pm:115`.
- Test run: synthetic 3-row AVI TSV, 4 MYBPC3 sites, local cache 111. Command and files are in the session scratchpad, not kept. Results:
  - The ALT is matched per allele: `chr11:47350100 G>A` got the G>A row (25.1/1.5), not the G>C row.
  - A SNV missing from the file and a 1-bp deletion both come out with **empty** fields.
  - `chr`-prefixed input with a `chr`-prefixed file works.

**Earlier lab work:** `/data_lab_PGP/resources/annotation/alphagenome/` holds July 2025 live-API results for MYBPC3, LMNA, FLNC, DSP and DSC2: per-variant CSVs and plots. MYBPC3 has 52 files. They come from the 2025 model and API, not from the Atlas, so they are not comparable with AVI. The folder also contains a `.env` file that was not opened. It should not be copied anywhere.

## Fit with this pipeline

**How it would be wired (PROPOSAL):**

- `--plugin AVI,file=${params.vep_avi},raw=1` in `modules/local/vep.nf`, next to `AlphaMissense` (`vep.nf:45`).
- `params.vep_avi` in `nextflow.config`.
- `AVI.pm` copied into the VEP_plugins folder.
- File under `/data_lab_PGP/resources/annotation/AlphaGenome/AVI/<release>/`, like ReMM. `/references/` is read-only.
- Do not use `links=1`: it adds a long URL to every transcript row.

**Status contract** (same as `parse_remm`):

- `AVI_status`: `scored` / `not_covered` / `not_applicable` (non-SNV, decided by `VARIANT_CLASS`) / `not_annotated` / `error`.
- The plugin leaves the fields empty both for an indel and for a missing SNV, so the parser must decide the status from `VARIANT_CLASS`, never from the empty field.
- Record the file version in the manifest. DeepMind changed the recommended scorer quantiles on 2026-06-18 (forum announcement), so PHRED values can move between releases.

**Evidence policy:**

- Informational only. It does not touch `NEW_IMPACT`, `PRIORITY_TIER` or the ACMG criteria.
- PHRED 20 means "top 1% of all possible SNVs genome-wide". That is not a Mendelian or cardiac threshold.
- For missense variants, AVI includes AlphaMissense, which is already a column here. Treat the two as correlated, never as two independent PP3 lines.

**Cost (INFERENCE):**

- The file is about the size of CADD's `whole_genome_SNVs.tsv.gz` (87.5 GB, same 6-column layout), and CADD already runs in the same `VEP_ANNOTATE` call. Expect a time increase comparable to one more CADD lookup, not a new process.
- Disk: about 88.5 GB zipped, plus the extracted bgzip if the zip does not already hold `.gz` + `.tbi`. `/data_lab_PGP` has 1.4 TB free (`df`, 2026-10-19).

**Data egress:**

- The local-file route sends nothing outside.
- The API routes (live inference, Atlas `query`/`annotate`) send cohort variant coordinates to Google. The global rules forbid that without explicit authorisation, and the Module 3 plan does not cover it.
- With the local file, whether compute nodes have network access (the Module 3 plan's main blocker) no longer matters.

## For the thesis hotspot map (MYBPC3 first)

**INFERENCE.** Taken from a local tabix query over the MYBPC3 locus, the AVI file gives every possible SNV a score with no API call and no compute. This is a predictor-agnostic saturation map: the same shape as `tabix` on CADD or ReMM over the same window, so the predictors can be compared side by side. It is the strongest reason to take the file, more than the per-sample annotation.

The two non-commercial files could add mechanism:

- the splicing scores, as a fourth splicing predictor next to SpliceAI, Pangolin and SpliceVault;
- the SHAP breakdown, to show which of the 18 feature categories drives a hotspot.

That is research use (allowed for CNIC), but not clinical, so it belongs in an analysis script rather than in the default Module 2 table. Any published figure or table derived from the splicing or SHAP files needs a conspicuous notice that it is subject to the AlphaGenome Output Terms of Use, plus a note of any modifications made (terms, prohibition 5). Figures derived from AVI alone are exempt. The cost of the SHAP file (283.9 GB) should be weighed only once it is clear the map needs it.

## Downloaded files (update 2026-10-19, evening)

Mario downloaded all three zips by browser into `/data_lab_PGP/resources/annotation/AlphaGenome/`.

Each zip was processed with `resources/prepare_alphagenome_downloads.sh` in two steps:
- `check`: lists the contents and writes the sha256 next to the zip.
- `extract`: unzips into `<zip_stem>/`, verifying every file's CRC.

The zips are kept. DeepMind publishes no checksums. The AVI and splicing tables are dated 2026-08-27, the SHAP table 2026-09-01. The only model version string in any download is in the SHAP file name: `combined_ag_cond_linear_ensemble_20260417` (INFERENCE: the AVI model behind the AVI scores too).

| Data | Extracted file | Zip sha256 |
|---|---|---|
| AVI | `AVI/avi_scores_snvs_tabix/alphagenome_variant_impact_score_snvs.tsv.gz` (+ `.tbi`), 88,470,182,310 bytes | `a237c198bc1c033da127fe0257129fdf246e7aa311f8fb584bb053ffa0c7966b` |
| Splicing | `combined_splicing_snvs_tabix/combined_alphagenome_splicing_snvs.tsv.gz` (+ `.tbi`), 20,635,924,975 bytes | `21011fe376455c7692f7f8f24ae9898d8f758060a91a24b15c35ac1c6cffa7a3` |
| SHAP | `avi_feature_importances_snvs_tabix/combined_ag_cond_linear_ensemble_20260417_feature_importance_indels_with_am_snvs.tsv.gz` (+ `.tbi`), 283,872,737,358 bytes | `fc65ceae024bf18eb06a4de85f7e2f981e6fef066819df96f6edb8932e18ab33` |

Run times on the workstation over the mount:
- AVI: `check` 21 min, `extract` 63 min.
- Splicing: `check` 6 min, `extract` 6 min.
- SHAP: `check` 50 min, `extract` 111 min.

### AVI (FACT, pysam tabix queries)

- **Format**
  - Header: `#CHROM POS REF ALT raw_score PHRED`, the column order `AVI.pm` parses, so the plugin reads the file as delivered.
  - Contigs: `chr1`–`chr22`, `chrX`, `chrY`. There is no `chrM`, so mitochondrial variants will be `not_covered`.
  - Coordinates are 1-based, starting at chr1:10001.
- **Coverage**
  - Genome-wide, intergenic positions included: chr5:90,000,001 is scored.
  - MYBPC3 gene span in GENCODE v46 (chr11:47,331,405–47,352,702): all 21,298 positions × 3 ALTs = 63,894 SNVs scored.
  - 8,427 of those SNVs (13.2%) have PHRED ≥ 20.
- **Comparison with the Atlas website** (signed-in browser, public example variants):

  | Variant | File | Website | Agreement |
  |---|---|---|---|
  | chr11:5225450 A>T (HBB, the `AVI.pm` example) | 21.49414 | 21.5 | matches after rounding |
  | chr11:47342906 C>A (the highest PHRED in MYBPC3) | 59.63756 | 59.7 | off by about 0.06 |

  The second difference is not a rounding effect. Possible causes are a different quantile table on the website or extreme-tail precision. It is too small to matter for ranking.
- **Feature breakdowns shown on the website**
  - HBB variant: the top feature is Cactus 241-way phyloP (32%). Splicing, RNA-seq, polyadenylation and PhastCons are also listed.
  - MYBPC3 C>A: the top feature is "Protein Termination" (45%), a binary VEP stop-gained feature.
  - **Consequence for the evidence policy:** AVI also folds in conservation and the VEP consequence class. For many variants it is therefore correlated with phyloP, CADD and our own consequence-based `NEW_IMPACT`, not only with AlphaMissense. It stays informational.

### Splicing (FACT)

- **Format**
  - Header: `#CHROM POS REF ALT alphagenome_splicing`, one value per ALT.
  - The values look like 0–1.
- **Coordinates**
  - Contigs: `chr1`–`chr22`, `chrX`, `chrY`.
  - Positions are 1-based. The REF at chr1:65409–65411 matches the GATK GRCh38 FASTA.
- **Coverage: gene regions only**
  - The MYBPC3 region is dense, chr11:47,331,329–47,352,704.
  - Positions in gene deserts are absent (chr5:90,000,001, chr13:55,000,001, chrX:100,000,001).
  - For the pipeline, an intergenic SNV is `not_covered`, never 0.
- **Integration**
  - VEP `--custom` cannot read this TSV.
  - Two options: a small plugin derived from `AVI.pm` (5 columns instead of 6), or a pysam tabix lookup after VEP.

### SHAP (FACT, pysam tabix queries)

- **Format**
  - Header: `#CHROM POS REF ALT` + 18 feature contributions: `MERGED_SPLICING`, 9 `MAX_ABS_*` AlphaGenome modalities (ATAC, contact maps, DNase, ChIP TF, ChIP histone, CAGE, PRO-cap, RNA-seq, polyadenylation), `ALPHAMISSENSE`, `CACTUS_241_WAY`, `PHASTCONS_470_WAY`, `PROTEIN_TERMINATION`, `START_LOST`, `STOP_LOST`, `IS_INSERTION`, `IS_DELETION`.
  - Signed values (some negative).
  - Despite `indels` in the file name, the indel flags are 0 and no row has a multi-base REF or ALT in MYBPC3: SNVs only.
- **Coverage**: same as AVI. Contigs `chr1`–`chr22`, `chrX`, `chrY`; starts at chr1:10001; intergenic sites scored; MYBPC3 gene span 21,298 positions / 63,894 SNVs, the same counts as AVI.
- **Example**: chr11:47342906 C>A (MYBPC3 top AVI) has `PROTEIN_TERMINATION` 1.549, `MERGED_SPLICING` 1.062, `CACTUS_241_WAY` 0.570, matching the website's "Protein Termination" as top feature.
- **Relation to the AVI raw score**: on 129 SNVs (MYBPC3, chr5 intergenic), `raw_score − sum(SHAP)` ranged from −0.063 to −0.031. INFERENCE: roughly additive, but not an exact constant offset, so the SHAP sum must not stand in for the AVI score.

## Implementation (2026-10-20)

- **Plugin**: `AVI.pm` taken from Ensembl VEP_plugins at commit `86bc0f033647be0c8db879e0cad317c158bd83de` (2026-09-08, "adding link param"). Installed as an untracked file in the shared clone `/data_lab_PGP/resources/annotation/VEP/VEP_plugins/` (branch `release/116`), sha256 `6dd3f3d964c3259832473c50fe53ed0e7cc37b9558cd2940e89d53705df91f76`.
- **Pipeline**:
  - `modules/local/vep.nf`: `--plugin AVI,file=${params.vep_avi},raw=1` after AlphaMissense.
  - `nextflow.config`: `vep_avi` points at the extracted AVI table.
  - `run_annotate_vcf.sh`: `vep_avi` recorded in the run manifest as an external resource.
  - `resources/all_but_old_gnomad_vep_cols.txt`: `AVI_PHRED`, `AVI_RAW` kept.
  - `src/python/parsing/modules/predictors.py`: `parse_avi` turns both columns into floats and adds `AVI_status`: `scored` / `not_covered` (SNV absent from the file, never 0) / `not_applicable` (non-SNV by `VARIANT_CLASS`) / `not_annotated` (columns absent) / `error` (unparsable value; scores blanked).
  - `filter_variants.py`: `parse_avi` runs after `parse_remm`; the three columns sit next to the ReMM ones.
  - `scoring.py`: availability audit entry "AlphaGenome AVI (informational)". No effect on `NEW_IMPACT`, `PRIORITY_TIER` or ACMG.
- **Checks (FACT)**:
  - `tests/python/test_avi_status.py`: 4 passed. The full local runner has 2 failing modules (`test_provenance_stamp`, `test_tag_transcript_priority`, parquet metadata comparisons); both fail the same way on committed HEAD, so they are not caused by this change.
  - `nextflow config` resolves `vep_avi`; `bash -n run_annotate_vcf.sh` passes.
  - Local VEP 111 (`vep.sif`, offline cache 111) with the shared plugin dir and the real file: chr11:47342906 C>A 59.63756 / raw 3.378, C>G 34.16779 / 1.92, chr11:5225450 A>T 21.49414 / 0.6643, chr5:90000001 A>G 10.35304 / 0.1648, all equal to the file; the same value on every transcript row and on the intergenic row. A 1-bp deletion (chr11:47342950 TG>T) and chrM:3000 A>G got empty fields. A rerun adding `--regulatory --canonical --mane --biotype --numbers` gave the HBB SNV the same value on its RegulatoryFeature row as on its transcript rows.
- **Provenance**: `write_run_manifest.py` has no release reader for AVI (`database_releases`), and the path holds no version. It records a `.sha256` sidecar for files over 1 GB, so one was written next to the table: `alphagenome_variant_impact_score_snvs.tsv.gz.sha256`. The download date (2026-10-19), the zip sha256 and the model id (`combined_ag_cond_linear_ensemble_20260417`, from the SHAP file name, INFERENCE for AVI) are recorded only here.
- **Not checked**: a cluster fixture run, and the added VEP time per partition.

## Remaining work

1. ~~Download, sha256, inspect, compare with the website~~: done for all three (sections above).
2. ~~Keep or delete the zips~~: Mario deleted them on 2026-10-20. The `.contents.txt` and `.sha256` sidecars stay next to the extracted folders.
3. ~~Implement~~: done 2026-10-20 (section below). Still to do: the fixture cluster run (Mario), and the VEP time cost on a WGS.
4. Tell Module 3 that its plan's bulk-access statement is outdated. Its AVI columns would be two new ones, informational only.

## Sources

- AlphaGenome downloads, terms, Atlas: `https://deepmind.google.com/science/alphagenome/{downloads,terms}` (read 2026-10-19).
- Ensembl blog, 2026-09-08: "AlphaGenome Variant Impact scores integrated into Ensembl VEP".
- `https://github.com/Ensembl/VEP_plugins/blob/main/AVI.pm`.
- AlphaGenome community forum, thread 964 "Offline availability of precomputed scores" (DeepMind reply, 2026-09-10); thread 929 "Updating variant score quantiles".
- DeepMind blog, 2026-09-08: "AlphaGenome Atlas: a predictive map of every possible DNA letter change in the human genome".
- `github.com/google-deepmind/science-skills`, `skills/alphagenome_variant_impact_score/SKILL.md` (Atlas API `region` saturation scan, ~10 variants/s for `annotate`).
