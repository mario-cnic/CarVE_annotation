# Manifest: database releases, tool versions, and a provenance file next to each table

Written 2026-10-05T13:10:00+02:00. Prompted by: Module 3 cannot tell a stale Pangolin result (mouse db, S223) from a valid one, and the manifest only recorded paths and hashes. Labels: FACT (measured or read now), INFERENCE, PROPOSAL.

## State before (FACT)

- Launch manifest recorded each resource's path, size, mtime and sha256 (up to 1 GiB), env fingerprints, `vep_sif_sha256`, `pangolin_distance`. Not recorded: any release read from a database's own content, the VEP software version, and most databases VEP reads (their paths were literal in `modules/local/vep.nf`).
- The tables carry no stamp: the Parquet footers of `S223` and `F2005_union` hold only `ARROW:schema` and `pandas`. S223's manifest has no `external_provenance` at all.

## What changed

- `nextflow.config` + `modules/local/vep.nf`: every database file VEP reads is now one param (`vep_cadd`, `vep_utrannotator`, `vep_spliceai_snv`, `vep_spliceai_indel`, `vep_revel`, `vep_mastermind`, `vep_loftool`, `vep_pli`, `vep_dbnsfp`, `vep_alphamissense`, `vep_maxentscan_dir`, `splicevardb_vcf`, `splicevault_file`, `gnomad_vcf`; plus `clinvar_vcf`). `run_annotate_vcf.sh` passes them all to the manifest as resources.
- `src/python/write_run_manifest.py`:
  - `launch.databases`: for each database, `release` (read from the database's own content), `details`, and `name_hint` (version-looking tokens from the file name, kept apart because a file name is not evidence). Readers: Pangolin db (`directives` table + `<db>.BUILD_MANIFEST.txt`), VEP cache (`info.txt` sources: ClinVar, dbSNP, GENCODE, SIFT, ...), dbNSFP (readme: version, release date, ClinVar release), VCF headers (`fileDate`, `source`, `reference`) for ClinVar, gnomAD, SpliceVarDB. A reader that cannot determine a release returns `release: null` and a reason; nothing raises.
  - `launch.tool_versions.vep_image`: SIF labels (read from the file header) and VEP's own version from `vep --help` inside the container (best effort).
  - `completed`: when the exit code is 0, writes `<table>.provenance.json` next to every `tables/*.pq` and `genotypes/*.genotypes.pq`: run id, times, full-manifest path and sha256, code commit and dirty flag, parameters (SpliceAI/Pangolin distance, SPiP tier), input path and hash, `databases`, resource identities, tools (VEP image, bcftools, Nextflow, env fingerprints and key packages). Never overwrites an existing file. The table itself is not hashed; the file is tied to it by name, size, mtime and the manifest's sha256.
- `tests/python/test_write_run_manifest_releases.py` (13 tests).

## Checks

- Real Nextflow render of `VEP_ANNOTATE` (`.command.sh`, scratch wrapper, fake `vep` on PATH, `-profile local_dev`, scratch `-w` and `-log`): identical before and after the param refactor. Compared with the real cluster `.command.sh` of the `F2005_union` run (work dir `work/70/03db39...`), the only differences are the FASTA path (`local_dev`), the file names, and the new ClinVar `--custom` line. This also closes the "Nextflow-rendered command" gap listed in the ClinVar walkthrough.
- `launch` against the local reference mapping (not the cluster): 67 s, mostly hashing files up to 1 GiB. Releases read: Pangolin db "evidence-based annotation of the human genome (GRCh38), version 45 (Ensembl 111)"; ClinVar `2026-09-28`; dbNSFP "dbNSFP version 4.8a" (ClinVar 20240215 inside it); VEP cache 111 (sources listed); VEP software `111.0`. Not readable from content: gnomAD (header has no `fileDate`; only the name hint `v4.1`), SpliceVarDB, CADD, REVEL, AlphaMissense, SpliceAI scores, MasterMind (name hint `2024.07.03`), LoFtool, pLI, SpliceVault, UTRAnnotator: identity and name hint only. The SIF label `org.opencontainers.image.version` is `latest`, so it is recorded but says nothing.
- `completed` on a scratch run folder: provenance files written for the main table and the genotype table (14 KB each); a second call did not overwrite them.
- `pytest` (manifest, external provenance, releases, VEP column tests): 34 passed.

## Not verified / limits

- Not run on the cluster: the launcher end to end, the VEP singularity call for `vep --help` there, and the time cost of the extra hashing on the login node.
- The provenance file travels with the table only when the table is read by path. An uploaded `.pq` on its own loses it. Embedding the stamp in the Parquet file metadata at the pipeline's last write would fix that; it changes process inputs and the DAG, so it is not done (PROPOSAL).
- Existing run folders (S223, F2005_union) were not touched and have no provenance files. S223 would have `unknown` for every database anyway (no `external_provenance` in its manifest).
- If a run folder is resumed and a table is rewritten, the existing provenance file is kept (not overwritten) and no longer matches; delete it to regenerate.
- Release `null` means "not determinable from the file", not "unversioned". Module 3 should treat it as unknown.
