# ClinVar custom annotation (`ClinVarVCF*`)

Written 2026-10-05T12:05:00+02:00. Prompted by Module 3's `PLT-075` (ClinVar evidence effectively off for Module 2 tables). Labels: FACT (measured or read now), INFERENCE, PROPOSAL.

## Why

FACT, S223 table (22,624,536 rows, `nf_work/annotation_out/S223/S223.parsed.clean.pq`):

- `clinvar_clnsig` (dbNSFP 4.8a, ClinVar 20240215 per its readme) is filled on 11,446 rows, all SNVs, 99.4% with missense or splice consequences. dbNSFP covers nsSNVs and splice-site SNVs only, so no indel ever gets it. `panel8_trio`: 0 of 1,015 rows.
- `CLIN_SIG` (VEP cache 111, ClinVar 202306 per the cache's `info.txt`; cache dirs 110/111/112 hold 202301/202306/202310) has no review status and is multi-valued. 1,368 rows have a "pathogenic" term (conflicting mixes included), only 300 of them also have `clinvar_clnsig`, 163 are indels.
- `_extract_clinvar_status` (`src/python/parsing/modules/acmg.py:15`, same list in Module 3) has no `CLIN_SIG` in its column list and stops at the first column that exists.
- Neither release is in `RUN_MANIFEST.json` (it records VEP cache version 111 only).

## What changed (Module 2 only)

- `modules/local/vep.nf`: `--custom file=${params.clinvar_vcf},short_name=ClinVarVCF,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNSIGCONF%CLNREVSTAT%CLNDN%CLNVC%ALLELEID`. The short name is `ClinVarVCF`, not `clinvar`, because `clinvar_CLNSIG` would differ from dbNSFP's `clinvar_clnsig` only by case (the `PLT-021` hazard).
- `nextflow.config`: `clinvar_vcf` = `${references_dir}/genomes/Homo_sapiens/annotations/clinvar/GRCh38/clinvar_20260928.vcf.gz`.
- `resources/all_but_old_gnomad_vep_cols.txt`: seven new columns (the parser keeps only whitelisted CSQ columns).
- `run_annotate_vcf.sh`: `clinvar_vcf` added to the manifest resource list (size, mtime, sha256 if under 1 GiB, `.md5` sidecar content).
- `resources/build_clinvar_vcf.sh`: downloads NCBI's current VCF to a local directory, checks NCBI's md5, names it by the header `fileDate`, indexes it, writes a `.md5` sidecar, refuses to overwrite, prints the manual deploy commands.
- `resources/column_dictionary.md` section 5, README annotator list, `tests/python/test_vep_custom_columns.py`.

## Checks

Static (local, ran): `pytest tests/python/test_vep_custom_columns.py` 3 passed (ClinVar columns whitelisted, param declared and used, no case-insensitive clash); `bash -n` on the build script; `nextflow config -profile local_dev` shows `clinvar_vcf` (Nextflow 24.04.2); `nextflow run annotate_vcf.nf -profile local_dev -preview` compiled.

Build script (local, ran): produced `clinvar_20260928.vcf.gz`, md5 `707723f0a08d1d711eb1c1de660a9098` = NCBI's; the second run refused to overwrite.

VEP (local, ran, NOT the pipeline): the pipeline's own `vep.sif` with cache 111, the pipeline FASTA and only the new `--custom` plus `--check_existing`, on a 15-record chr-prefixed VCF made from ClinVar records (MYBPC3, one chrX, one chrM) in three small runs. The test VCFs and outputs are in the session scratchpad, not in the repo.

| Case | Result |
|---|---|
| Pathogenic SNV, likely pathogenic SNV, benign, VUS, conflicting | Right `CLNSIG`, `CLNREVSTAT`, `ALLELEID` |
| Deletion `CG>C`, duplication `A>AT`, indel `CC>AG` (all Pathogenic) | Matched (the old sources never match indels) |
| Same position and REF, other ALT (`A>T`, `A>G`, `A>C`) | Each gets its own ClinVar record: allele-specific |
| `chr11:47340005 T>A`, no ClinVar record in the window | All seven columns empty |
| Input `chr`-prefixed, ClinVar file unprefixed | Matched, including `chrM` to `MT` and `chrX` (VEP `--synonyms`) |
| chrX record with no germline classification | `ClinVarVCF` and `ClinVarVCF_ALLELEID` filled, `CLNSIG` empty |
| Vendored `vcf_parser_pysam.py` on the VEP output with the new whitelist (`vcf_parser` env) | The seven columns present, values identical to CSQ, no case-insensitive column clash |

## Not verified

- The cluster run and a full-size run. (The Nextflow-rendered `VEP_ANNOTATE` command was later checked with a real render: identical to the cluster's `.command.sh` apart from the new ClinVar line, FASTA path and file names; see `20261005_manifest_databases_and_table_provenance.md`.) The standalone VEP call above has the same `--custom` text but fewer flags (no `--everything`, no plugins, no `--fork`); a clash between the new `--custom` and the other plugins is untested.
- Time and memory cost of a ClinVar lookup on 5 M records (the file is 198 MB; expected small next to SpliceAI/Pangolin, not measured).
- `filter_variants.py` and the second parse stage on these columns (the table-level step after the parser).

## Location (updated 2026-10-05T12:30:00+02:00)

`/references` is read-only for the agent and for Mario, so the file lives in the project's resource tree, next to LaBranchoR and SpliceVault: `/data_lab_PGP/resources/annotation/clinvar/GRCh38/clinvar_20260928.vcf.gz` (+ `.tbi`, `.md5`). Built there directly by `resources/build_clinvar_vcf.sh` (download plus NCBI md5 check, no copy from local storage): md5 `707723f0a08d1d711eb1c1de660a9098`, 198,125,411 bytes; `bcftools view -r 11:47332113` on it returns the three ClinVar records. `nextflow.config` `clinvar_vcf` points there; the container already binds `/data_lab_PGP`. Not yet read by a cluster node. A newer release means a new file name and an edit of `clinvar_vcf`; old files stay. No manual deployment step remains.

## Remaining (PROPOSAL)

- Derive a `ClinVar_Stars` value from `ClinVarVCF_CLNREVSTAT` (0 to 4) in Module 2 or Module 3: one place, not both. Module 3's scoring reads `ClinVar_Stars` (`scoring.py:85` in the vendored copy).
- Module 3 reads `ClinVarVCF_CLNSIG` explicitly (status per row), not the first existing column. Their call.
- Tell Module 3: ClinVar fields are multi-value with `&`, commas became `&` in `CLNREVSTAT` and `CLNDN`, and `CLNDN` must not be split.
- Existing parsed tables (S223, F2005) do not have these columns; they need a re-run of VEP or a backfill (not designed).
