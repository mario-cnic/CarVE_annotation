# Gene constraint from gnomAD v4.1.2 (pLI, LOEUF, missense Z)

2026-11-04T16:35:07+01:00. Code base `e59316e` plus the uncommitted change described here. Mario approved replacing the VEP pLI plugin with gnomAD v4.1 constraint as one source (2026-11-04). FACT unless labelled.

## Why

- `pLI_gene_value` came from the VEP pLI plugin with `plI_gene.txt`, ExAC r0.3 (March 2016): its values match `fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt` in the same folder (AGRN 0.17335, NOC2L 1.3e-19). The plugin matches by gene symbol. Example: EMSY is `C11orf30` in ExAC (pLI 0.99999), so every EMSY row had no pLI.
- `LOEUF_gene_value` and `mis_z_gene_value` were never emitted (no commit adds them). Module 3 reads them (`features.py:77-85`), so its gene prior used pLI only and its mis_z ≥ 3 branch (`scoring.py:399-401`) never fired.

## Source

- `gs://gcp-public-data--gnomad/release/4.1.2/constraint/gnomad.v4.1.2.constraint_metrics.tsv.bgz` (https://storage.googleapis.com/…; 574,538,527 bytes, uploaded 2026-09-23). The GCS md5 matched; sha256 `a3fbca52…138a858f`. Downloaded 2026-11-04 into the session scratchpad only. It is public and can be downloaded again; keeping a copy on lab storage is Mario's call.
- The gnomAD changelog (2026-09-28) says v4.1.2 is a constraint-only release. It fixes the rank and bin fields of v4.1.1; "the constraint metric values themselves (including LOEUF) are unchanged". Field descriptions come from `gnomad.v4.1.1.gene_constraint_ht_README.md`; no v4.1.2 README exists.
- Transcripts are GENCODE v39, GRCh38. VEP 111 uses GENCODE 45; the join is by Ensembl gene ID, which is stable across the two releases, so transcript versions do not matter.
- The LOEUF field is `lof.oe_ci.upper` (high-confidence pLoF). It is the field the rank and bin columns use, and the one Ensembl's `LOEUF.pm` reads for v4. On unflagged MANE Select transcripts, 15.06% fall below the v4.1.1 percentile-15 cut (0.44896); `lof_hc_lc.oe_ci.upper` gives 15.00%.
- gnomAD guidance (v4.1.1 news post, 2026-03): for a pLoF-constrained cut use LOEUF < 0.45, which replaces v2's 0.35 and sits at about the 15th percentile. The help page prefers continuous use, and says that v2 LOEUF thresholds do not give an equivalent number of genes on v4.

## What changed

- `src/python/build_gene_constraint_table.py` builds `resources/gnomad_v4.1.2_gene_constraint.grch38.tsv` and a `.manifest.json` beside it. The table has one row per ENSG: its Ensembl MANE Select transcript, or the canonical one when there is no MANE Select.
  - Counts: 20,076 genes; 18,411 by MANE Select and 1,665 by canonical. RefSeq rows (Entrez `gene_id`) are dropped.
  - Columns: pLI (`lof.pLI`), LOEUF, LOEUF decile, lof Z, mis_z, `constraint_flags`, `gene_flags`.
  - Values of flagged transcripts are kept: 2,235 genes carry a flag. 1,820 genes have no pLI or LOEUF and 87 no mis_z. These stay empty, never 0.
- `annotate_gene_constraint` (`predictors.py`) joins on VEP `Gene` (ENSG only, no symbol fallback). It writes the names Module 3 reads:
  - `pLI_gene_value`, `LOEUF_gene_value`, `mis_z_gene_value` (float64)
  - `gene_constraint_flags`, `gene_constraint_transcript`
  - `gene_constraint_status`: `matched`, `not_in_table`, `not_applicable` or `not_annotated`
- Wiring:
  - The table is loaded once per process.
  - `--gene-constraint-table` is passed by `TSV_TO_TABLE`. The new param `gene_constraint_table` is in `run_annotate_vcf.sh`'s resource list, so the run manifest hashes it.
  - `--plugin pLI` is removed from `vep.nf`, along with `vep_pli` and `pLI_gene_value` in the VEP column list.
  - **VEP tasks re-run instead of resuming** on the next `-resume`.
- M2's own gates are unchanged (`pLI ≥ 0.9` in `acmg.py` PVS1 and `scoring.py`); they now read gnomAD v4.1.2 pLI. Moving them to LOEUF is interpretation, which is Module 3's (DEC-0008).

## Checks (local)

- `tests/python/test_gene_constraint.py`: 5 passed. Covers the join, NaN kept, statuses, the refusal of a repeated gene, and the build choosing MANE over canonical, Ensembl only, flags and the manifest. Full runner: only the two known local parquet-metadata failures (`test_provenance_stamp`, `test_tag_transcript_priority`).
- `concat_fixture_1` rebuilt locally: VCF → TSV with the new column list, then `filter_variants.py` whole and with `--chunk-rows 20`.
  - The whole and chunked outputs are equal, frame and schema.
  - Versus the published table, 20 genes: 19 `matched`, LINC00200 `not_in_table`.
  - EMSY: pLI none → 1.0 (LOEUF 0.27). Its 13 rows move M2 Tier 4 → 3 (score 15 → 25) through M2's +10 pLI bonus.
  - PRKCQ: 0.97 → 0.004 (LOEUF 0.58). ITGB1: 0.91 → 1.0.
  - No `ACMG_CRITERIA` change. The other 2 tier changes are the ITGB1 UTR rows of PRED-13.
- `F2005_union`, per gene, read-only (old = the table's ExAC value):
  - 61,106 ENSG on transcript rows; 19,394 are in the constraint table. The 41,712 not in it are mostly non-coding genes, 5.9 M rows.
  - pLI present: 16,830 genes → 17,909. 890 lose it: 819 have no gnomAD estimate, 71 are not in the table. 1,969 gain it.
  - pLI ≥ 0.9: 3,132 genes → 3,628 (740 drop out, 1,236 enter). LoF rows whose M2 PVS1 gate flips: 261 lost, 537 gained, of 8,580.
  - Curated "Prioritary" genes (`Complete_gene_list_V4.csv`): 2,926 in F2005. 69 have no v4.1.2 pLI (19 not in the table: H19, MEG3, KCNQ1OT1, snoRNAs; the others have no estimate: small ribosomal protein genes, KCNE5); ExAC lacked 214. pLI ≥ 0.9: 828 → 904. Examples: PTEN 0.98 → 0.34 (LOEUF 0.60), FHOD3 1.00 → 0.000 (LOEUF 0.62), BAG3 and CALM1-3 now ≥ 0.9.
- Not run: the cluster fixture.

## Consequence for Module 3 (INFERENCE from its code)

When a table built with this code reaches Module 3, its gene prior applies v2-era cuts to v4 values (LOEUF < 0.35 / < 0.60, mis_z ≥ 3; `scoring.py:399-401`). PVS1 (`acmg.py:311`) gates on the new pLI. No Module 3 code change is needed for this to happen. On v4.1.1 MANE Select, LOEUF < 0.36 is about the 10th percentile and < 0.60 about the 25th (gnomAD's table). Proposed record: `carve-platform/decisions/DEC-0012-gene-constraint-source.md`.

## Pending

- Commit, then **one** cluster fixture covering this change and PRED-13, with `concat_fixture_1`'s flags: `bash run_annotate_vcf.sh -profile standard --input_vcf test_data/raw_vcfs/clinvar_af_fixture.vcf.gz --annotation_cache_dir /mnt/scratch/mruizp/annotation_cache_test --chunk_size 5 --gene_restricted_chunk_size 2 --run_id <id>`.
- Pass, against `nf_work/annotation_out/concat_fixture_1/tables/concat_fixture_1.parsed.clean.pq`, joined on (`Locus`, `Feature`), 199 rows:
  - 6 new columns: `UTRAnnotator_status`, `LOEUF_gene_value`, `mis_z_gene_value`, `gene_constraint_flags`, `gene_constraint_transcript`, `gene_constraint_status`;
  - `pLI_gene_value` differs on 181 rows, `PRIORITY_TIER`/`VARIANT_PRIORITY_SCORE`/`TIER_RATIONALE` on 15 (2 ITGB1 UTR rows, 13 EMSY rows), `SCORE_EXPLANATION` on 18, `NEW_IMPACT` on 2 (ITGB1); nothing else (`TRANSCRIPT_PRIORITY_TIER` is present in both cluster tables);
  - the order of keys inside `5UTR_annotation` may differ, and a UTRAnnotator class can rarely flip (`PRED-15`);
  - gene rows are `matched` or `not_in_table` with the resource table's values; the VEP CSQ header has no `pLI_gene_value`; the run manifest lists `gene_constraint_table` with its sha256.
- The local expected table was built in the session scratchpad and will not survive; to rebuild it from the repo root:
  ```bash
  SU=/data_lab_PGP/shared/utils; OUT=<scratch dir>
  PYTHONPATH="$SU/conda_envs/vcf_parser/lib/python3.13/site-packages" $SU/conda_envs/vcf_parser/bin/python src/python/parsing/vcf_parser_pysam.py \
    --input nf_work/annotation_out/concat_fixture_1/annotation/concat_fixture_1.annotated.vcf.gz --output $OUT/scratch.tsv \
    --vep_columns resources/all_but_old_gnomad_vep_cols.txt --add_info --add_vep --overwrite --logging_level WARNING
  $SU/conda_envs/datasci/bin/python3 src/python/parsing/filter_variants.py --input $OUT/scratch.tsv --output $OUT/whole.pq \
    --hgnc-table /data_lab_PGP/resources/annotation/hgnc/hgnc_complete_set.2026-09-29.txt --disease-gene-list $SU/data/Complete_gene_list_V4.csv \
    --spliceai-symbol-map resources/spliceai_symbol_to_ensg.grch38.gencode45.tsv --enst-spip-map resources/enst_to_spip_nm.grch38.gencode45.tsv \
    --dbnsfp-aligned-columns resources/dbnsfp_transcript_aligned_columns.txt --gene-constraint-table resources/gnomad_v4.1.2_gene_constraint.grch38.tsv \
    --skip_quality_filter   # add --chunk-rows 20 for the chunked check
  ```
- Tables built before this change (e.g. `F2005_union` 2026-11-03, `F1830_union_avi`) keep ExAC pLI and have no LOEUF; DEC-0012's Module 3 review applies to tables built from this code.
- Module 3 threshold review before the next LAB promotion or M3 export of a table built with this code (Mario decides).
