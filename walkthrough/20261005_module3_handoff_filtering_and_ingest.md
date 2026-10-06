# Module 3 handoff: filtering strategy and ingest of WGS-size tables

Written 2026-10-05T11:05:00+02:00 in `annotation_pipeline_new` (Module 2). Module 3 and `carve-platform` were not edited from here. Labels: FACT (measured or read now), INFERENCE, PROPOSAL.

## Question

WGS parsed tables are about 5-6 GB. Can Module 3 load and prioritise them, and where should filtering live (separate filter module, Module 3, all variants, split by analysis type)?

## Facts

- Real tables: `S223.parsed.clean.pq` 22,624,536 rows x 591 columns, 5.0 GB, 22 row groups, 5,068,263 distinct `Locus`. `F2005_union.parsed.clean.pq` 23.7 M rows x 614 columns, 6.0 GB (3 samples).
- Rows are variant x VEP transcript. Row counts of the S223 table, from `src/tools/row_funnel.py` (reads 4 columns, ~10 s):

| Step | Rows | Share |
| --- | --- | --- |
| all rows | 22,624,536 | 100% |
| `TRANSCRIPT_PRIORITY_TIER` <= 2 | 2,073,038 | 9.16% |
| + `gnomADv4_AF_grpmax_joint` <= 0.001 or missing | 12,995 | 0.057% |
| + `IMPACT` HIGH or MODERATE | 249 | 0.001% |

  Tier 3 rows are 91% of the table (tier values: 3 = 20,551,498; 2 = 2,036,246; 1 = 36,792). One sample, rows not variants; the AF and impact thresholds are illustrative, not validated for cardiomyopathy. Missing AF is kept (1.3% of rows).
- Module 3 ingest today (verified now): `src/api/routers/variants.py:72` `content = await file.read()`; `:83` `pd.read_parquet(temp_path)`; `src/engine/parquet_converter.py:100` `df.to_parquet(...)` of a full DataFrame. `DuckDBManager` (`src/engine/duckdb_manager.py:33`) calls `convert_to_parquet` when given a DataFrame and queries `read_parquet(...)` afterwards (`:220`, `:278`).
- `docs/DUCKDB_NATIVE_INGESTION_ROADMAP.md` (Module 3, line 82) proposes `bcftools view -i 'AF<threshold'` as a prefilter and says WGS scale is out of scope (line 134).
- Accepted decisions: `DEC-0003` (all VCF filtering and conversion is Module 2; interpretation is Module 3). D-5 verdict (c): both row shapes, tag, do not filter. Prep notes `carve-platform/reviews/2026-09-24_seguimiento/prep-notes.md` ("Architecture review for discussion - filtering and backend at 100 WGS") hold the provisional backend recommendation this note follows.
- The genotype table is records x samples (F2005: 5,263,518 records x 3 = 15,790,554 rows).
- In the Python used for this check, `import duckdb` failed once with a `libstdc++` `CXXABI_1.3.15` error (the interpreter and library path differed between two invocations; cause not investigated). Benchmark Module 3 in its own environment.

## Proposal (not decided; needs Mario's verdict and a `DEC-` record in `carve-platform`)

Keep every variant in a lossless store; Module 2 tags, Module 3 filters at query time. No separate filter module (new cross-cutting structure), and no destructive prefilter anywhere.

Roadmap, in order:

1. **Module 3 - ingest by path.** Open a Module 2 run's `.pq` in place with DuckDB (column projection, row-group pruning); no `file.read()`, no pandas load, no rewrite. Needs an allowlisted run root (same concern as `PLT-002`). Scoring must not run at upload.
2. **Module 3 - default view.** Tier 1-2 rows by default (`PLT-136`), other tiers reachable. Remove the `AF<threshold` prefilter from `DUCKDB_NATIVE_INGESTION_ROADMAP.md`.
3. **Module 2 - flag columns, tag only.** QC status, frequency bin with missing kept distinct from zero, consequence class, panel membership, tier reason (`PLT-135`). Optionally a materialised candidate Parquet next to the full one. Depends on `PLT-130` (sites / variant_gene / variant_transcript split) and `PLT-040` (chromosome partitioning). Column names not designed yet.
4. **Module 3 - case-specific filters.** De novo, compound heterozygous, recessive and X-linked on the long genotype table (`gt_status = alt_plus_other_allele` meaning in `20261002_genotype_tracking.md`); fix the headerless `.ped` read (`pd.read_csv` needs `header=None`, three `pedigree.py` copies).
5. **Named analysis views** over the one store: diagnostic panel, secondary findings, splicing / non-coding, cohort carrier search. CNV/SV stays a separate track (`DEC-0005`).
6. **Both - benchmark 1, 10, up to 100 WGS** on the intended hardware: ingest time, peak RAM, scratch, query latency (patient/family review, cohort carrier lookup), candidate-view recall on reviewed cases. Up to 100 WGS is a target from the prep notes; D-6 is still open.
7. **Both - genotype table at cohort scale** (INFERENCE: records x samples reaches hundreds of millions of rows at 100 samples). A sparse design is needed, with no-call kept distinct from hom-ref.

## Open for Mario

- Accept the tag-in-Module-2 / query-in-Module-3 split and record it as a `DEC-` (owner: `carve-platform`).
- Which flag columns and thresholds are the default candidate view (needs cardiomyopathy-specific justification; none is validated here).

## Not done

Nothing was changed in Module 3 or `carve-platform`. No cluster run. The funnel covers one sample's rows.

## Addendum 2026-10-05T12:40:00+02:00: ClinVar for Module 3 (answers `PLT-075`, question 5)

Details and checks: `20261005_clinvar_custom_annotation.md`. Module 2 implemented this; Module 3 and `carve-platform` were not edited from here.

**FACT, why `CLIN_SIG` / `clinvar_clnsig` fail** (S223, 22.6 M rows):
- `clinvar_clnsig`, `clinvar_review`, `clinvar_id`, `clinvar_trait` come from dbNSFP 4.8a (ClinVar 20240215 per its readme). Filled on 11,446 rows, all SNVs (99.4% missense or splice). Never on indels. `panel8_trio`: 0 of 1,015 rows.
- `CLIN_SIG` comes from the VEP cache 111 (ClinVar 202306 per its `info.txt`). Any variant type, but no review status, and `&`-joined multi-values including conflicting mixes.
- A row can carry two different ClinVar releases. 1,368 S223 rows have a "pathogenic" term in `CLIN_SIG`; only 300 also have `clinvar_clnsig`; 163 of the 1,368 are indels.
- `_extract_clinvar_status` (`acmg.py:15`, same list in Module 2's vendored copy) has no `CLIN_SIG` and breaks on the first column that exists. `ClinVar_Stars` is read by `scoring.py:85` and nothing populates it.

**New in Module 2 output (next run onwards; existing tables do not have it):**
Source: ClinVar GRCh38 VCF release 2026-09-28 (`clinvar_20260928.vcf.gz`, md5 `707723f0a08d1d711eb1c1de660a9098`), VEP `--custom`, exact allele match, SNVs and indels.

| Column | Meaning |
|---|---|
| `ClinVarVCF` | ClinVar Variation ID |
| `ClinVarVCF_CLNSIG` | Aggregate germline classification |
| `ClinVarVCF_CLNSIGCONF` | Counts per classification for conflicting records, e.g. `Uncertain_significance_(1)&Benign_(1)` |
| `ClinVarVCF_CLNREVSTAT` | Review status |
| `ClinVarVCF_CLNDN` | Disease names |
| `ClinVarVCF_CLNVC` | Variant type |
| `ClinVarVCF_ALLELEID` | ClinVar Allele ID |

**Semantics (FACT, local VEP test):**
- No ClinVar record: all seven empty. Record without a germline classification (243,032 of 4,555,206 records have `CLNSIG=.`): `ClinVarVCF` filled, `ClinVarVCF_CLNSIG` empty. Do not treat that as benign or as "no data".
- Values are variant-level, repeated on every transcript row of the variant.
- Encoding: VEP turns `|` and `,` into `&`. `CLNREVSTAT` is `criteria_provided&_single_submitter`, not with commas; `CLNSIG` multi-values are `Benign&risk_factor`, `Pathogenic&_low_penetrance`; `CLNDN` names that contain commas break into several `&` pieces, so do not split `CLNDN` on `&`. Match review status against the known strings.
- `CLNSIG` values in this release (counts of records): Uncertain_significance 2,357,513; Likely_benign 1,132,231; Benign 193,977; Pathogenic 183,134; Conflicting_classifications_of_pathogenicity 178,063; Likely_pathogenic 128,251; Benign/Likely_benign 82,464; Pathogenic/Likely_pathogenic 43,769; plus small groups (`drug_response`, `risk_factor`, `Affects`, `Pathogenic|other`, `VUS-high`, ...). `CLNREVSTAT` values: criteria_provided,_single_submitter 3,262,207; ..._multiple_submitters,_no_conflicts 740,754; criteria_provided,_conflicting_classifications 177,766; no_assertion_criteria_provided 101,674; reviewed_by_expert_panel 22,228; no_classification_provided 6,690; practice_guideline 51; and two tiny "no classification" groups.

**PROPOSAL for Module 3 (their call):**
1. Read `ClinVarVCF_CLNSIG` explicitly, with a status per row (no record / record, no germline classification / classified), not "first existing column". Whether `CLIN_SIG` stays as a fallback is a separate choice; mixing it with `clinvar_clnsig` mixes ClinVar releases.
2. Stars from `ClinVarVCF_CLNREVSTAT` in one place only (Module 2 or Module 3). The usual ClinVar scheme is practice_guideline 4, reviewed_by_expert_panel 3, criteria_provided multiple submitters no conflicts 2, criteria_provided single submitter 1, criteria_provided conflicting classifications 1, every "no assertion / no classification" status 0; confirm it against NCBI's review-status page, not checked here.
3. Pathogenic matching on the new column: `Pathogenic`, `Likely_pathogenic` and `Pathogenic/Likely_pathogenic` count; a `Conflicting...` value must not; values such as `Pathogenic&risk_factor` need an explicit rule.
4. Panel tables and the S223/F2005 tables lack these columns until Module 2 re-runs VEP or a backfill is designed.

**Not verified:** nothing has run on the cluster; the parser stage after the column parse (`filter_variants.py`) was not run on these columns.

## Addendum 2026-10-05T13:15:00+02:00: provenance file next to each table (stale vs valid)

Details and checks: `20261005_manifest_databases_and_table_provenance.md`. Module 3 and `carve-platform` were not edited from here; this is a cross-module contract and needs a `DEC-` record.

**Problem (FACT):** the tables carry no version stamp (Parquet footers hold only `ARROW:schema` and `pandas`). S223's Pangolin columns hold mouse gene IDs (`ENSMUSG...`, `PLT-139`); nothing in the table says so.

**New in Module 2 (runs launched after this change, exit code 0):** next to every `tables/*.pq` and `genotypes/*.pq` a file `<table>.pq.provenance.json`, JSON, never overwritten:
- `run`: `run_id`, launch and completion times, `exit_code`, path and `sha256` of the full `pipeline_info/RUN_MANIFEST.json`.
- `code`: commit and dirty flag. `parameters`: SpliceAI and Pangolin distance, SPiP tier. `input`: path and sha256.
- `databases`: per database `release` (read from the database's own content), `details`, `name_hint` (version tokens from the file name only). `release: null` means not determinable, not unversioned. Releases readable today: Pangolin db (GENCODE 45 for the current db, sha256 `55831b7679132a7dd868f8aa65a5ae6cdfa298993a29460c846f276b3d5eeb8c`), ClinVar (`fileDate`), dbNSFP (version, date, ClinVar release), VEP cache (`source_*` list incl. ClinVar, dbSNP, GENCODE). gnomAD, SpliceVarDB, CADD, REVEL, AlphaMissense, SpliceAI scores, MasterMind, LoFtool, pLI, SpliceVault, UTRAnnotator: identity and name hint only.
- `resources`: path, size, mtime, sha256 (null for files over 1 GiB).
- `tools`: VEP image (SIF labels, VEP version), bcftools, Nextflow, and per env a fingerprint and key package versions (Pangolin 1.0.2, torch, ...).

**Linking (FACT):** the file is tied to its table by name, size, mtime and the full manifest's sha256; the table itself is not hashed. It only travels with the table when Module 3 reads tables by path. An uploaded `.pq` alone loses it; embedding the stamp in the Parquet file metadata is an open Module 2 follow-up.

**PROPOSAL for Module 3 (their call):**
1. Treat a table as valid for a predictor only when its provenance file is present and the predictor's database matches an accepted release (Pangolin: the sha256 above); no file = unknown, not valid. Existing runs (S223, F2005_union) have no provenance file.
2. Keep the `ENSMUSG` check on the raw Pangolin columns as a hard stale flag for tables without a provenance file.
3. Show the release (or "unknown") next to the score, instead of a bare value.

## Addendum 2026-10-05T14:15:00+02:00: the stamp is also in the Parquet footer

Module 2 now also writes the provenance into the file metadata of every published table, so a plain upload carries it. Details: addendum in `20261005_manifest_databases_and_table_provenance.md`. Needs a `DEC-` record in `carve-platform` (not written from here).

- Key `carve.provenance`, JSON text, in the Parquet footer of `tables/<run>.parsed.clean.pq`, `tables/<run>.parsed.clean.withGT.pq` and `genotypes/<run>.genotypes.pq`. Read: `pq.ParquetFile(path).metadata.metadata[b"carve.provenance"]` (footer only) or DuckDB `parquet_kv_metadata`.
- Fields: `run_id`, `launched_at`, `code{commit,dirty,dirty_patch_sha256}`, `parameters` (incl. `spliceai_distance`, `pangolin_distance`, `spip_tier`), `input{path,sha256}`, `databases{<name>{release,name_hint,details}}`, `resources{<name>{path,size_bytes,mtime_utc,sha256,hash_skipped}}` (Pangolin: `resources.pangolin_db.sha256`), `tools{vep_image,bcftools,nextflow_version,envs}`, `full_manifest`.
- The stamp holds launch-time facts only. The sidecar `<table>.pq.provenance.json` adds the table's size and mtime and the run's completion (exit code, manifest sha256).
- Absent key = unknown, not valid: runs made before this change (S223, F2005_union), runs started without `run_annotate_vcf.sh`, and the legacy per-gene path.
- The key is plain text; anyone can copy it onto another file. It protects against mistakes, not tampering.
- Order of trust for Module 3: footer stamp (works for uploads and by-path), else the sidecar next to the table (by-path only), else `ENSMUSG` heuristic on the raw Pangolin columns (stale flag only).

## Addendum 2026-10-05T16:05:00+02:00: scope is the 3-sample, 6 GB family table, not 100 WGS

Mario, 2026-10-05: the near-term target is tables like `F2005_union` (2 WES + 1 WGS, 6 GB), already waiting in Module 2 for filtering and prioritisation. Measurements and plan: `20261005_family_table_candidate_funnel.md`. Short version: DuckDB queries the 6 GB file in place in 0.1-5 s at about 1 GB RAM; a candidate selection in SQL gives 5,402 rows (1,020 loci) from 23.7 M rows (illustrative thresholds); only those rows need Module 3's pandas scoring. Replaces the 100-WGS framing above for planning; the 100-WGS items stay as later design targets.
