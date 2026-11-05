# Review session 2026-11-04/05: what was fixed, how it was checked, what is left

Written 2026-11-05 in `CarVE_annotation` (Module 2), from a cloud session with no cluster access. Labels: FACT (code read, test or command run in this session), INFERENCE, PROPOSAL.

Per-defect detail (failure scenario, fix, test): `BUG_TRACKER.md`. Table and contract changes for Module 3: `walkthrough/20261105_module3_table_contract_changes.md`. This note does not repeat them.

## Branches

| Branch | Base | Content | State |
| --- | --- | --- | --- |
| `ccr-8536b091-fabr1g` | `main` | Every fix below except the LaBranchoR work, one commit per entry | Pushed; not merged; not run on the cluster |
| `pred18-labranchor-positional` | `ccr` | PRED-18 option (a): LaBranchoR as a positional, informational annotation; the Branchpointer columns are removed | Pushed, rebased on `ccr` |
| `labranchor-ism` | `pred18` | PRED-20 (ISM builder rewrite) and the allele-aware `LaBranchoR_ISM_*` columns | Pushed, rebased on `pred18`; needs the GRCh38 ISM table built on the cluster (`TODO.md` step 7d) |

Merge order: `ccr`, then `pred18`, then `labranchor-ism`. They stack, so merging only `ccr` is a valid stopping point.

## What changed, by area

- **SpliceAI and the score cache**
  - PRED-16: a SpliceAI exception fails the chunk instead of writing the record unscored.
    - `annotation_cache.py migrate` copies the scored sites of the old cache version into the new one.
  - MISC-35: a resumed run lists the cache segments of its own version.
  - MISC-41: `-resume` reruns a task whose code changed. A `# code ${params.code_md5.<key>}` line is in every process that runs repo code. Cache version ids are unchanged.
  - MISC-44: the same for large resources read by path (`# resources ${params.resource_id.<key>}`: path, size and mtime of every file), and `codeMd5` now follows a symlinked directory.
- **Gene-restricted tier input**
  - PRED-17 / SV-1: symbolic, breakend and `*` records are kept out of SpliceAI/Pangolin/SPiP.
  - MISC-22: an input with no site in the gene BED still gives one, header-only, chunk.
  - MISC-4: gene BED builder: case-insensitive and Ensembl-ID matching, plus a primary-assembly fallback. The BED itself is not rebuilt.
  - MISC-9: `GENE_SUBSET`'s real script is tested with bcftools.
- **Integrity checks**
  - MISC-6: CONCAT and MERGE compare CHROM/POS/REF/ALT, not only record counts.
  - MISC-2: the run QC report compares the table's distinct `Locus` with the merged VCF.
  - MISC-34: the input check refuses an input that already carries predictor INFO tags.
  - MISC-39: the input check refuses a symbolic ALT without an `INFO/END` header line.
- **Parser** (VCF_TO_TSV, TSV_TO_TABLE)
  - Silent losses that now fail or are kept:
    - MISC-28: a record without CSQ is kept.
    - MISC-29: a per-record error fails the task.
    - MISC-36: SpliceVault/MaxEntScan errors fail the task.
  - Wrong values:
    - MISC-33: `HGVSc` "nan".
    - PRED-7: SpliceVault cryptic-site position.
    - PRED-8: intron offset from UTR, non-coding and range HGVS.
    - MISC-37: categories of an all-null first chunk.
    - MISC-30/31/32: genotype columns not used by the pipeline.
  - Robustness and dead options:
    - MISC-21 / PRED-2: one-row dbNSFP crash; `parse_dbnsfp` removed.
    - PRED-19.
    - MISC-27.
    - MISC-26: `withGT` pandas metadata.
    - MISC-43: `parse_missense` shrank `MISSENSE_COLUMNS`.
- **Scoring and ACMG**
  - ACMG-4: one PVS1 veto list.
  - ACMG-7: the explanation and rationale follow the veto.
  - ACMG-8: one definition of the core panels.
  - ACMG-9: `AGL`, not `AGGL` (Mario's decision).
- **Launcher and provenance**
  - MISC-23: relative paths, and `--outdir` honoured.
  - MISC-42: the manifest and the cache take the params Nextflow resolves for the actual command line.
  - MISC-38: a republished table gets its own provenance sidecar.
  - ACMG-6: the disease-gene list's sha256 is in the manifest.
- **Configuration**
  - MISC-40: `local_dev` resource limits, so a `-stub-run` works on a laptop.
- **Tools outside the run**
  - LIFT-1/2: hg19 liftover keeps samples and checks strand, span and REF.
  - PRED-20: ISM builder, on `labranchor-ism`.
- **Tests and docs**
  - C3 and MISC-24: tests import what production imports.
  - DOC-1: CADD v1.7.
  - README: SV paragraph, liftover, `--outdir`, resume.
  - New end-to-end and fuzz tests: VCF → TSV → table (whole and chunked) → transcript tier.
- **Corrections to earlier tracker text**
  - MISC-5 was fixed by `7027950`. This was first misattributed on a shallow clone.
  - PRED-3/4/5 are confirmed fixed.
  - PRED-7 diagnosis corrected: the defect was the reference point, not the strand.
  - MISC-37 severity lowered to Low.

## How it was checked

- **Unit and integration tests** (FACT, `ccr` at `ac6cf35`):
  - pandas 2.3.3 venv: 393 passed, 1 skipped.
    - The skipped test needs `bcftools annotate -W`, which the local bcftools 1.19 lacks; the cluster versions have it.
  - Default `python3` with pandas 3.0.5: 9 failed (MISC-25, open). The cluster `datasci` env is pandas 2.2.
  - Where a tracker entry says a test "fails on the old code", that was run against the previous version of the file, not assumed.
- **End-to-end and fuzz** (FACT):
  - A synthetic VCF goes through `vcf_parser_pysam.py` → `filter_variants.py` (whole, and chunked at 1, 2, 3, 4 and 100 rows) → `tag_transcript_priority.py`.
  - Chunked output equals whole output (frame and Arrow schema).
  - Seeded random records give the same result.
  - The end-to-end tests found PRED-19 (MISC-27's `transform` test) and MISC-37 (a first record without CSQ, chunks of 2).
- **Local Nextflow** (FACT, Nextflow 25.10.2, `-stub-run`, `-profile local_dev`, synthetic chr14 input):
  - The full DAG completes with and without the score cache (5 hits / 16 misses).
  - The CONCAT identity checks pass.
  - `GENE_SUBSET` leaves out `<CNV>`.
  - After editing a script, `-resume` reruns only the tasks using it (MISC-41).
  - Through `run_annotate_vcf.sh`: the manifest records the command-line FASTA with its sha256, the QC report is written, and a resume is 39/39 cached.
  - Stub runs execute no predictor, VEP or real Python step. They check wiring, channels, caching and the launcher only.
- **Not checked here**: anything that needs the cluster.
  - Real SpliceAI/Pangolin/SPiP/VEP.
  - Real data sizes and run time (MISC-6 key comparison on a WGS).
  - The cache migration on the real cache.
  - The BED rebuild.
  - The ISM table build.
  - Module 3's reads (repo not available).

## What is left

- **Cluster checks**: `TODO.md`, item "Cluster checks for branch `ccr-8536b091-fabr1g`", steps 1-8 with 3b, 5b and 7b-7d.
  - Step 3, the SpliceAI cache migration, must run before the first run with this branch. Otherwise every gene-restricted site is recomputed.
  - Then flip the 🟡 entries in `BUG_TRACKER.md`.
- **Module 3 and `carve-platform`**: `walkthrough/20261105_module3_table_contract_changes.md` (draft DEC text inside).
- **Decisions for Mario** (left untouched, by instruction):
  - ACMG-1/2/3: which genes and criteria.
  - ACMG-5.
  - ACMG-6: version the curation list in git or not.
  - PRED-1, PRED-5, PRED-9 (two REVEL sources), PRED-10, PRED-15.
  - MISC-25: support pandas 3 or pin 2.x.
  - A threshold for `LaBranchoR_ISM_*`, if they should ever affect tiering.
  - Whether a resumed run should rewrite `TABLE_STAMP.json`.
- **Done after this note was first written**: MISC-44, resource identity in the `-resume` key (`TODO.md` ordered item 6).
  - Large references read by path (FASTA, VEP cache and plugin data, ...) now enter the task hash by path, size and mtime.
  - Checked with local stub runs. Its launch-time cost on the cluster is step 7a of the checklist.
  - Still open in item 6: conda envs and the provenance stamp on resume.
