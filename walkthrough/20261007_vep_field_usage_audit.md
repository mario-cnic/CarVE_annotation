# VEP field usage audit: redundant population annotations and oversized CSQ

Written 2026-10-07T14:59:23+02:00. Scope: audit requested by Mario; annotation configuration and existing results are unchanged. FACT = read or measured; INFERENCE = interpretation; PROPOSAL = future change.

## Verdict

**INFERENCE:** Cached gnomAD is redundant for the current v4-based workflow and should stop being generated for new runs, after frequency fallbacks are made explicit. The larger field-count problem is `--plugin dbNSFP,...,ALL`; the larger population-field storage block is custom gnomAD itself. Removing only cached gnomAD will not make the table compact.

“Never be used” cannot be established from source inspection. This audit separates identifiable consumers, compatibility fallbacks, mechanical parsing, and no explicit code reference. Generic column selection/export and future scientific questions can still use an unreferenced field. Proposed removals are relative to the current Mendelian prioritization and cohort-burden workflows.

## Evidence and reproducibility

- Baseline: `08c7bbeef751276eabec60368ee4849941fe3819`, initially clean tree. Module 3 inspected read-only; commit/status in the manifest.
- Producer: `modules/local/vep.nf:81` (`--everything`), `:90` (`dbNSFP,ALL`), `:96` (custom gnomAD). Deprecated bash equivalent: `src/hpc/annotate_vep_vars.sh:47-62`.
- `nextflow.config:101` selects `resources/all_but_old_gnomad_vep_cols.txt`. Despite its name, the file includes **all 20 cached gnomAD frequency fields**.
- Observed CSQ header: existing `panel8_trio_layout` annotated VCF. Storage: footer of existing family table, 23,686,438 rows × 614 columns, 23 row groups. These are different historical artifacts; field presence in each is recorded separately. This is not a claim that the family run emitted every field in the newer VCF header.
- Reference metadata: read-only local mapping of VEP 111 and dbNSFP 4.8a. No full reference/data scans or large hashes.
- Script: [audit_vep_field_usage.py](../src/tools/audit_vep_field_usage.py). Outputs: [field inventory and references](../reports/vep_field_audit_20261007_final/fields.tsv), [unreferenced candidates](../reports/vep_field_audit_20261007_final/unused_candidates.tsv), [summary](../reports/vep_field_audit_20261007_final/summary.json), [manifest](../reports/vep_field_audit_20261007_final/manifest.json). Earlier directories from this session are preliminary snapshots; `_final` is authoritative.

The script reads VCF/dbNSFP headers and Parquet metadata, then searches tracked source in Module 2 and Module 3 (Python, R, shell, Nextflow, frontend). References include comments, producer declarations and transformation lists: they are **not** proof of interpretation use. Generic `af`/`freq` matching and arbitrary-column display/export can access unnamed fields. Parsing, scoring, reports, candidate-set and burden code were additionally read to interpret dependencies. No patient rows were read/exported.

## Inventory (FACT)

**664 CSQ slots, 659 unique names; all 456 dbNSFP header fields emitted.** The parser selector has 563 entries but only 559 unique names. Of observed CSQ fields, 506 have no explicit source-code token reference; 399 of those are dbNSFP fields. These screening counts are not automatic deletion counts.

Storage is the sum of existing compressed Parquet column chunks, in decimal MB. It is not a measured rewrite saving or runtime improvement.

| Group | Unique CSQ names | Family-table columns | MB | Assessment |
| --- | ---: | ---: | ---: | --- |
| Cache gnomAD exomes/genomes | 20 | 20 | 297.98 | Superseded for current frequency workflow; report fallbacks remain |
| Cache `MAX_AF`, `MAX_AF_POPS` | 2 | 2 | 30.16 | Older-source aggregate, not custom v4; fallback references remain |
| Cache 1000 Genomes AF | 6 | 6 | 107.26 | Generic `AF` used by legacy filters/plots |
| dbNSFP population fields | 228 | 106 | 16.29 | No explicit semantic consumer identified; 134 unselected by parser |
| dbNSFP rank scores | 57 | 49 | 6.71 | Five have downstream references; 52 pruning candidates |
| Custom gnomAD v4 | 78 | 78 | 2471.64 | Twelve names explicitly referenced, including source identifier |

Cache gnomAD plus its old maximum = **328.14 MB**, 5.46% of 6.013 GB compressed column chunks. With 1000 Genomes, the older cache frequency block = **435.40 MB**, 7.24%. Custom v4 = **41.10%**; its 66 names without explicit references occupy **1927.38 MB**. Proportions are specific to this historical table.

## Cached gnomAD: retire in new annotation (PROPOSAL)

The cache's `111_GRCh38/info.txt` records **gnomADe r2.1.1**, **gnomADg v3.1.2**, **1000 Genomes phase3**, **ClinVar 202306**. Custom gnomAD is configured as `v4.1/gnomAD.v4.1.vcf.gz`; the previous manifest investigation could not determine its release from header content. Distinguish configured release identity from independently verified file provenance.

Retire these cache names:

- `gnomADe_AF`, `gnomADe_AFR_AF`, `gnomADe_AMR_AF`, `gnomADe_ASJ_AF`, `gnomADe_EAS_AF`, `gnomADe_FIN_AF`, `gnomADe_NFE_AF`, `gnomADe_OTH_AF`, `gnomADe_SAS_AF`.
- `gnomADg_AF`, `gnomADg_AFR_AF`, `gnomADg_AMI_AF`, `gnomADg_AMR_AF`, `gnomADg_ASJ_AF`, `gnomADg_EAS_AF`, `gnomADg_FIN_AF`, `gnomADg_MID_AF`, `gnomADg_NFE_AF`, `gnomADg_OTH_AF`, `gnomADg_SAS_AF`.
- `MAX_AF`, `MAX_AF_POPS`: custom gnomAD v4 **does not feed VEP's `MAX_AF`**.

VEP 111's [Config.pm](https://raw.githubusercontent.com/Ensembl/ensembl-vep/release/111/modules/Bio/EnsEMBL/VEP/Config.pm) confirms `--everything` enables `af`, `af_1kg`, `af_gnomade`, `af_gnomadg`, `max_af`. Removing explicit gnomAD flags does nothing here: they are not explicitly present. Replace `--everything` with individual options, omitting the three gnomAD/maximum flags. Keep `--check_existing` for identifiers and existing-variation information.

Conservative first-stage equivalent, retaining 1000 Genomes and other `--everything` behavior:

```text
--sift b --polyphen b --ccds --hgvs --symbol --numbers --domains
--regulatory --canonical --protein --biotype --af --af_1kg --pubmed
--uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna
```

Remove the retired names from the selector too. Selector-only pruning shrinks the table but leaves VEP generating the annotations in the VCF. `--fields` restricts output; avoiding options/plugins addresses annotation generation.

**Dependencies:**

- Module 2 `parsing/modules/acmg.py:130`: `MAX_AF` fallback when custom grpmax column is absent. `parsing/modules/scoring.py:31-38`: first name containing `af`/`freq` if expected field absent.
- Module 2 report `generate_clinical_prioritization_report.py:77-90`: old-cache fallback if v4 grpmax absent **or entirely empty**. Module 3 `src/report_generator.py:80-98`: cached genome/exome AF for older tables without v4 columns.
- Legacy Module 2 `filter_and_summarize.py:129`, `generate_interactive_report.py:106`, `src/R/plot_annotation_results.R:105`: prefer generic `AF`. In CSQ this is 1000 Genomes, not gnomAD (related to existing MISC-1).
- Current Module 3 uses v4 grpmax with joint AF when grpmax undefined (`src/modules/features.py:121-127` and candidate-set code). Keep both fields.

New v4 output should preserve missing values rather than silently use older releases or measured zero. Historical-table compatibility can remain separately labelled.

## dbNSFP `ALL`: largest excess field count

**FACT:** Population block = 122 `gnomAD_exomes_*`/`gnomAD_genomes_*`, 48 `ExAC_*`, 36 `ALFA_*`, 12 `1000Gp3_*`, 4 `ESP6500_*`, two each of `ALSPAC_*`, `TWINSUK_*`, `UK10K_*`. None has an explicit semantic consumer. dbNSFP's readme specifies gnomAD **v4.0.0**, not cache v2/v3: these are two distinct redundancy problems. The dbNSFP readme also records June 3, 2024, GENCODE 29/Ensembl 94 variant determination and ClinVar 20240215.

**PROPOSAL:** Replace `ALL` with a positive field allowlist, excluding these 228 fields first. Of them, 134 are already unselected for the table but still generated in CSQ. The [Ensembl plugin documentation](https://static.ensembl.org/info/docs/tools/vep/script/vep_plugins.html) supports requesting named fields.

Additional candidates: unconsumed conservation/fitCons, eQTL/sQTL, archaic-human alleles, hg18/hg19 coordinates, alternate HGVS copies and unused predictor/rank outputs. Exact names and references are in the inventory. Do not infer redundancy solely from absence in scoring: downstream research has consumers too.

**Retain until coordinated replacement:**

- `Ensembl_transcriptid` and retained transcript-aligned score dependencies; `parse_dbnsfp_by_row_transcript` prevents attribution to the wrong transcript. Update the aligned-column list with the allowlist.
- `REVEL_score`, `MetaRNN_score`, `MetaLR_score` and attribution/status outputs where used.
- Downstream predictor inputs: `AlphaMissense_score`, `BayesDel_addAF_score`, `BayesDel_noAF_score`, `MPC_score`, `VEST4_score`, `ClinPred_score`, `PHACTboost_score`, and parsing/display inputs.
- Five rank scores in `src/downstream/01_build_master.R:125-134`, `02_cohort_enrich.R:219-223`: `REVEL_rankscore`, `AlphaMissense_rankscore`, `BayesDel_addAF_rankscore`, `MetaRNN_rankscore`, `PHACTboost_rankscore`. Some uses are coercion/order lists, still an existing output contract.

Standalone CADD/REVEL/AlphaMissense plugins are not automatically equivalent to dbNSFP versions: coverage, releases, indel support, naming and attribution differ. Address case collisions (`CADD_PHRED`/`CADD_phred`, `CADD_RAW`/`CADD_raw`) explicitly.

## Other candidates

| Annotation | Evidence | Proposal |
| --- | --- | --- |
| `Mastermind_MMID3`, `Mastermind_counts` | No consumer identified; 30.90 MB; repeated selector entries | Remove plugin if literature-index output is outside product scope |
| `LoFtool` | Producer references only; 14.23 MB | Remove plugin; keep used `pLI_gene_value` |
| `CADD_RAW` | No consumer identified; 62.38 MB | Suppress output, retain `CADD_PHRED`; may not save computation |
| `PolyPhen`, `DOMAINS`, `UNIPARC`, `GENE_PHENO`, `miRNA`, `TRANSCRIPTION_FACTORS` | No explicit consumer | Optional metadata pruning after scope review |
| Cache/dbNSFP ClinVar copies | Older releases; some fallbacks used | Standardize new interpretation on dated `ClinVarVCF*`; retain labelled old-table compatibility |

Do not remove UTRAnnotator, SpliceRegion/SpliceAI, SpliceVault, SpliceVarDB, MaxEntScan, transcript identity/MANE, core consequence/HGVS solely because individual details are not scored: attribution, coverage and evidence remain useful. Dropping `--regulatory` changes consequence/feature rows. Avoid `--pick` or other transcript-row pruning in this cleanup.

## Custom gnomAD: compact core and deliberate research detail

**FACT:** Eleven data fields have explicit consumers (each table name prefixed `gnomADv4_`): `AF_joint`, `AF_grpmax_joint`, `grpmax_joint`, `faf95_joint`, `fafmax_faf95_max_joint`, `AF_joint_nfe`, `faf95_joint_nfe`, `AC_joint`, `AN_joint`, `AC_joint_nfe`, `AN_joint_nfe`. AC/AN fields are not useless: `src/downstream/02_cohort_enrich.R:89-96`, `03_burden.R`, and spatial analysis use reference counts/denominators.

Three unconsumed homozygote age fields occupy **250.53 MB**: `age_hist_hom_bin_freq_joint` (215.01 MB), `age_hist_hom_n_smaller_joint` (20.27 MB), `age_hist_hom_n_larger_joint` (15.25 MB). Strong first pruning candidates.

**PROPOSAL:** Keep `exomes_filters`, `genomes_filters`, `nhomalt_joint` and suitable ancestry/count/coverage fields if population QC/recessive interpretation is intended, despite limited code usage. Detailed ancestry AF/FAF and separate exome/genome fields should be intentional optional research detail. An 11-field core preserves explicit consumers but may be too narrow for broader QC.

The custom request repeats `nhomalt_joint`, `exomes_filters`, `genomes_filters`: 80 requests, 77 unique data fields. Remove duplicates. CSQ repeats `TSL`, `APPRIS` from VEP/dbNSFP; the selector repeats these and both Mastermind fields. Later `fix_dup_columns` coalesces suffixed columns, which can conceal source differences.

## Next steps and checks

1. Remove repeated requests/selector entries and resolve source collisions.
2. Replace `--everything` minus cached gnomAD/maximum; update selector and frequency consumers. Retire 1000 Genomes only after fixing generic `AF` consumers.
3. Replace dbNSFP `ALL` with a reviewed allowlist; preserve attribution dependencies and research inputs.
4. Prune Mastermind/LoFtool/age histograms if excluded by scope; define compact and optional custom gnomAD detail.
5. Validate on a small cluster fixture in a **new run directory**: identical variant keys and transcript/consequence rows, retained field values and score/status outputs, except intentional changes. Include absent gnomAD, grpmax missing/joint present, SNVs/indels and transcript matching. Inspect CSQ collisions; test old-table compatibility separately. Notify Module 3 before changing the contract.

**Executed locally:** audit script, header/reference/footer inspection, source review, syntax and aggregate consistency checks. **Not executed:** VEP/predictor workload, cluster jobs, table rewrite. Runtime/memory improvement and actual rewritten size remain unmeasured. Configuration is unchanged.

## Addendum 2026-10-07T15:15:12+02:00: approved cleanup implemented

Mario's decision supersedes the proposals above: remove cached gnomAD and VEP maximum AF, remove homozygote age histograms but keep alternate-homozygote counts, remove duplicates; **keep Mastermind and LoFtool** because they will be used. Other dbNSFP pruning is pending Mario's review.

Changes:

- `modules/local/vep.nf` and deprecated bash producer: explicit VEP 111 options replace `--everything`, excluding `af_gnomade`, `af_gnomadg`, `max_af`. All other bundled options, including 1000 Genomes and regulatory annotations, remain.
- Custom gnomAD: the three `age_hist_hom_*` requests are removed; joint, ancestry, grpmax, exome and genome `nhomalt_*` counts remain. `nhomalt_joint`, `exomes_filters`, `genomes_filters` are each requested once. 74 unique data fields plus the custom identifier remain.
- Parser selector: the 22 cache/maximum fields and three age fields are removed; repeated `TSL`, `APPRIS`, and Mastermind entries are deduplicated. 534 unique entries remain.
- `resources/vep_dbnsfp_fields.txt`: explicit list of all 454 dbNSFP 4.8a header fields except duplicate `TSL`/`APPRIS`. No score or population field is pruned. Both producers read this list; `nextflow.config` exposes it and `run_annotate_vcf.sh` records it as a manifest resource. The rendered values participate in the VEP task command/cache identity.
- `TSL`/`APPRIS` now come only from VEP, so they are removed from `dbnsfp_transcript_aligned_columns.txt`. This prevents dbNSFP transcript mismatches from blanking VEP metadata. No shared-utils or Module 3 code was edited.
- README and column dictionary describe the new schema. Existing result files are untouched. Legacy consumers' fallbacks remain for historical tables; current output still has their preferred custom v4 fields, including missing values.

Local validation:

- The explicit dbNSFP list was compared directly with the real database header: precisely `TSL`/`APPRIS` excluded, no other difference.
- Focused schema/row-order suite: **12 passed**. Attribution suite, run separately to avoid the known shared `modules` import collision: **26 passed**, including a new synthetic transcript-mismatch test preserving VEP `TSL`/`APPRIS`.
- Nextflow 24.04.2 command render/execution with **fake VEP and tabix only**, synthetic header-only input and isolated `/tmp/vep_cleanup_render_0ui4bjy7`: one local mock task succeeded. Argument checks confirm 454 dbNSFP fields, 74 unique custom data fields, retained homozygote counts and excluded cache-frequency flags; generated command passes `bash -n`. The installed capsule launcher could not parse Java `11.0.32.1`; the cached Nextflow dependency classpath was used directly. No dependency installation or predictor inference occurred.
- Bash producers/wrapper syntax and Git whitespace checks passed.

The removed cache/maximum and age columns accounted for **578.67 MB** of the historical family table's compressed column chunks. This is attributable storage, not measured size after a rewrite. Runtime, memory, actual new CSQ header and retained annotation-value parity remain cluster-unverified. Use a fresh fixture run/output directory; do not rewrite old results for this cleanup.
