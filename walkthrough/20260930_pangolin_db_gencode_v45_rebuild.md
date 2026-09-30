# Pangolin annotation db rebuilt from GENCODE 45 (MISC-12 / PLT-139)

**2026-09-30T10:40+02:00**

## Context

`BUG_TRACKER.md` MISC-12: Pangolin's `annotation_file`
(`/data_lab_PGP/shared/utils/pangolin_db/pangolin_grch38.db`) is a mouse annotation. Pangolin uses
it for gene assignment (`get_genes()`, `src/external/Pangolin-main/pangolin/pangolin.py:80-104`)
and for exon masking, so every Pangolin score either pipeline has produced is affected.

## FACT: what the old db is

Read-only SQLite inspection of the `directives` table:
`description: evidence-based annotation of the mouse genome (GRCm38), version M23 (Ensembl 98)`,
`date: 2019-09-06`. It holds 32,285 genes, 118,153 transcripts and 805,200 exons, written by
gffutils 0.14. It was **not modified**.

Old baseline on the panel VCF (`nf_work/annotation_out/panel7_test/panel7_test.annPangolin.vcf.gz`):
22/72 records carry a Pangolin annotation, all 22 with `ENSMUSG` gene IDs. The other 50 were
treated as outside any gene.

## Decision: annotation release (user's choice)

Evidence presented before the choice. Counts are on primary contigs. "Curated" means the 217-panel
transcripts in `resources/gene_transcript_mapping.txt`.

| | GENCODE 38 | GENCODE 45 | GENCODE 48 | MANE v1.4 |
|---|---|---|---|---|
| Ensembl release | 104 | 111 (= VEP cache) | 114 | 114 |
| Genes | 60,649 | 63,187 | 78,686 | 19,276 |
| Selected transcript = curated ENST | 210 (4 differ: FHOD3, TAZ, BRAF, RAF1) | 212/212 | 212/212 | 212/212 |
| `Ensembl_canonical` == `MANE_Select` | 17,774/17,774 | 19,204/19,204 | 19,276/19,276 | n/a (0 canonical tags) |

Five curated symbols resolve in none of them (e.g. `TAZ` → `TAFAZZIN`).

- The MANE GTF has no `Ensembl_canonical` tags. `create_db.py`'s default filter would therefore
  drop every transcript and exon, and masking would silently get empty exon lists.
- `parse_pangolin` (`shared/utils/src/modules/splicing.py:263-300`) takes the max |score| across
  every gene Pangolin reports and never joins gene IDs to VEP. What matters downstream is which
  genes are in the db, not the ID format.

**Chosen:** GENCODE 45, `Ensembl_canonical` (create_db.py default), one transcript per gene. It
uses the same gene models as the VEP 111 cache.

## Build

In `~/pangolin_db_build/` (local disk; authorized). The env is new, created with approval:
`mamba create -n pangolin_dbbuild -c conda-forge -c bioconda python=3.12 gffutils=0.14`. That
gffutils version matches the cluster reader `pangolin_env`. The cluster env itself can't be run
from the `noexec` mount.

```bash
cd ~/pangolin_db_build
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/MD5SUMS
grep " gencode.v45.annotation.gtf.gz$" MD5SUMS | md5sum -c -        # OK
~/apps/miniforge3/envs/pangolin_dbbuild/bin/python \
    <repo>/src/external/Pangolin-main/scripts/create_db.py gencode.v45.annotation.gtf.gz
```

| | |
|---|---|
| GTF | `gencode.v45.annotation.gtf.gz`, GENCODE 45 (Ensembl 111), header date 2023-09-19 |
| GTF md5 | `b6eeb6c9791b7a43a5504a654ff09d9a` (matches GENCODE `MD5SUMS`) |
| GTF sha256 | `56a166e4a43c934a1910c0d09dec786ac048862fc65130731d41f3c7f451ae68` |
| Builder | `create_db.py` @ `6d366a3`, unmodified |
| db | `gencode.v45.annotation.db`, 392,273,920 bytes |
| db sha256 | `55831b7679132a7dd868f8aa65a5ae6cdfa298993a29460c846f276b3d5eeb8c` |
| Run | 2026-09-30T10:16:58+02:00 → 10:19:24+02:00, 112 MB peak RSS |

The full record is in `~/pangolin_db_build/BUILD_MANIFEST.txt`, next to the env export
(`pangolin_dbbuild_env.yml`) and `validation.log`.

## Validation (local, 0 failures)

`src/tools/validate_pangolin_db.py --db … --gtf … --pangolin-py src/external/Pangolin-main/pangolin/pangolin.py`

- **Header.** The directives declare human GRCh38, GENCODE 45.
- **Counts.** Gene, transcript and exon counts equal an independent parse of the GTF with the same
  filter: 63,187 / 63,187 / 305,766.
- **Contigs and IDs.** Contigs are exactly `chr1`–`chr22`, `chrX`, `chrY`, `chrM`. All gene IDs
  are `ENSG`, with 0 `ENSMUS`.
- **Transcripts.** Exactly one transcript per gene, and every one is tagged `Ensembl_canonical`.
- **Spot genes.** MYBPC3, TTN, MYH7, LMNA, SCN5A, EMD, DMD and SRY cover `+`/`-` strands, chrX
  and chrY. For each, gene coordinates, strand and the full canonical exon list equal the GTF
  (e.g. MYBPC3 35 exons, TTN 363).
- **`get_genes()`.** Pangolin's own function is compiled from `pangolin.py` by AST, so torch
  isn't imported. It finds each spot gene on the right strand bucket with 2 × exon-count
  boundaries, for both `chrN` and bare `N` contig names. It returns nothing at intergenic
  `chr1:5000`.

## Not done here

- **Deployment.** The file is not copied to `shared/utils` (user's manual step; commands below).
  `pangolin_grch38.db` stays untouched.
- **Repointed consumers.** `nextflow.config`'s `pangolin_db` and `src/hpc/annotate_pangolin_vars.sh`
  now name `gencode.v45.ensembl_canonical.grch38.db`. Both pipelines fail to find the db until it
  is deployed; they don't silently fall back to the mouse file.
- **No inference.** Pangolin was not run, locally or on the cluster. The panel7 re-test is
  pending on the cluster.
- **Earlier results.** Every Pangolin score produced before this fix, including the `S223` WGS
  run, stays invalid and needs re-running.
- **Other repos.** Git-tracked files of the sibling repos under `pipelines/` and `shared/utils`
  mention the old path only in documentation. Untracked consumers were not searched.

## Pending manual steps (user)

**1. Deploy** (adds a new file; the old db is not touched):

```bash
D=/home/mruizp/data_lab_PGP/shared/utils/pangolin_db
cp ~/pangolin_db_build/gencode.v45.annotation.db  $D/gencode.v45.ensembl_canonical.grch38.db
cp ~/pangolin_db_build/BUILD_MANIFEST.txt         $D/gencode.v45.ensembl_canonical.grch38.BUILD_MANIFEST.txt
sha256sum $D/gencode.v45.ensembl_canonical.grch38.db   # expect 55831b7679132a7dd868f8aa65a5ae6cdfa298993a29460c846f276b3d5eeb8c
chmod a-w $D/gencode.v45.ensembl_canonical.grch38.db
```

**2. Pangolin-only re-test** (on a cluster compute node):

```bash
cd /data_lab_PGP/pipelines/annotation_pipeline_new && mkdir -p RUNS/pangolin_v45_retest
zcat test_data/raw_vcfs/panel7_test.vcf.gz > RUNS/pangolin_v45_retest/panel7.vcf
PYTHONPATH=src/external/Pangolin-main /data_lab_PGP/shared/utils/conda_envs/pangolin_env/bin/python3 -m pangolin.pangolin \
  RUNS/pangolin_v45_retest/panel7.vcf \
  /references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta \
  /data_lab_PGP/shared/utils/pangolin_db/gencode.v45.ensembl_canonical.grch38.db \
  RUNS/pangolin_v45_retest/panel7.annPangolin.vcf -d 10000
grep -v '^#' RUNS/pangolin_v45_retest/panel7.annPangolin.vcf | grep -o 'ENS[A-Z]*G' | sort | uniq -c   # expect ENSG only
```

The old baseline was 22/72 records annotated, all `ENSMUSG`. Expect `ENSG` IDs on most records.

**3. Full pipeline re-test** (about 4 h). This also exercises the MISC-14 pre-flight check:

```bash
bash run_annotate_vcf.sh -profile standard --input_vcf test_data/raw_vcfs/panel7_test.vcf.gz --run_id panel7_v45_retest
# then: nf_work/annotation_out/panel7_v45_retest/panel7_v45_retest.assembly_check.tsv should say PASS
```

**4. Afterwards:**
- Re-run Pangolin for every earlier result (incl. `S223`).
- Update `carve-platform` `STATUS.md` / `PLT-139`.
- Set MISC-12 to 🟢 once the re-test passes.
