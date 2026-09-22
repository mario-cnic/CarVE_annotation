# Nextflow migration — whole-VCF redesign (course correction)

**2026-09-22**

## Context

Phases 1-5 (see `walkthrough/20260921_nextflow_migration_phase1_spliceai_spike.md`) built the
Nextflow port around a **legacy per-gene adapter**: one gene's pre-split VCF in, five predictors
run in parallel, merged, flattened to a table — mirroring `main.sh`'s existing one-gene-per-file
convention exactly.

Mid-session course correction: the user redirected that the 217-gene panel and its per-gene-file
convention are a **historical artifact** of how this pipeline was originally iterated, not a
requirement going forward. The primary goal is now: take an arbitrary input (VCF, parquet, etc.)
that may span many genes, annotate as much of it as possible, and produce **one output file** for
the whole run — not one file per gene. Per-gene splitting becomes optional, low-priority, future
work. Full reasoning and the revised architecture are in the plan:
`/home/mruizp/.claude/plans/scalable-wibbling-snowflake.md` (substantially rewritten this session —
not a new file, since Phases 1-5's actual process ports remain valid, reusable building blocks;
only the entry point, filtering, and transcript-selection layer around them changes).

## What was built: `resources/v5_genes_loc.bed`

First piece of the revised plan (item 6): the region-restriction BED the gene-restricted tier
(SpliceAI/Pangolin/SPiP) will subset against, replacing the 217-gene panel's `cardio_genes_loc.bed`
with coverage for `Complete_gene_list_V5` (5,877 unique ClinGen-curated genes after dedup).

**Key finding that shaped the approach**: `gene_coords.py`'s live Ensembl REST API (today's only
gene→coordinate resolution mechanism in this repo) does not scale to ~5,877 genes in reasonable
time — one HTTP round-trip per gene, with a `3`s sleep between genes in `gene_coords.sh`'s loop.
Investigation found two resources already present **locally**, closing this gap with zero external
API calls:
- `shared/utils/data/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz` — genome-wide MANE GTF, 19,334 genes
  with `chr`-prefixed GRCh38 gene-body coordinates.
- `shared/utils/data/ensembl_to_refseq.tsv.gz` — 533,741-row Ensembl↔RefSeq/MANE mapping (not
  needed for the BED itself in the end — see "Simplification" below — but confirmed as a viable
  cross-check resource).

**New file**: `src/python/build_gene_restriction_bed.py` — reads `Complete_gene_list_V5.csv`,
dedups by `GENE`, joins against the MANE GTF's `gene` feature lines by `gene_name`, pads each
interval by `10000`bp (matching SpliceAI/Pangolin's own `-d`/`-D` window so splice-prediction-range
variants at a gene's boundary aren't clipped), writes `resources/v5_genes_loc.bed` (same
6-column format as the existing `cardio_genes_loc.bed`). Purely local/offline — no cluster, no
network — logs unmatched genes rather than failing on them.

**Simplification found during design**: transcript *priority tagging* (tier 1: 217-panel curated,
tier 2: VEP's own `MANE_SELECT`, tier 3: everything else — the other half of this course
correction, item 10 in the plan, not yet built) turns out to need **neither** the MANE GTF **nor**
`ensembl_to_refseq.tsv.gz` at runtime: `vcf_parser_pysam.py` already emits a `MANE_SELECT` column
per row when `--add_vep` is used (it's a standard VEP/CSQ field), so tier-2 detection is a plain
column check on the already-existing output, not a new resource lookup. The MANE GTF's only actual
job in this design is building the region-restriction BED (a coordinates problem, not a
transcript-identity problem).

## Verification (fully local, matches the plan's stated verification steps)

Ran `build_gene_restriction_bed.py` with `--unmatched-output` to also capture the miss list:

```
Read 5877 unique genes from Complete_gene_list_V5.csv
Parsed 19334 gene coordinates from MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz
Matched 5742/5877 genes (97.7%)
Wrote 5742 regions to resources/v5_genes_loc.bed
Wrote 135 unmatched gene names
```

Exactly matches the plan's pre-computed coverage numbers (verified during planning, reproduced
here as the actual build artifact, not just a dry-run estimate).

1. **Row/coordinate sanity**: 5,742 lines, all `chr`-prefixed, no negative BED starts (the
   `max(0, ...)` guard in the margin-padding step matters for genes near a chromosome's start).
2. **Known-gene spot check**: `MYBPC3` → `chr11:47321405-47362702` — exactly
   `47331406 - 1 - 10000` to `47352702 + 10000`, matching the MANE GTF's raw gene-body span
   (`chr11:47331406-47352702`, confirmed directly during planning) with the margin applied
   correctly on both sides.
3. **217-panel overlap check**: of the 164 genes already resolved in `cardio_genes_loc.bed`, 163
   also appear in the new V5 BED — the one gap (`C10ORF71`) is the same gene-symbol-casing alias
   already flagged in the unmatched list (`C10orf71` in the MANE GTF's modern naming), not a new
   problem this script introduces.
4. **New finding, not predicted during planning**: 9 genes (`SLC37A4`, `ORAI1`, `POLR2A`, `NBPF1`,
   `MROH8`, `APOBEC3A`, `SHANK3`, `C4A`, `ASPN`) matched by name but their MANE GTF `gene` line sits
   on an **alt/fix contig** (e.g. `chr11_KN196481v1_fix`), not the primary chromosome assembly —
   confirmed no duplicate entry exists for any of these 9 on a primary chromosome. A `bcftools view
   -R` against a real clinical VCF (called on primary chromosomes only, as essentially all are)
   would never match these regions, silently excluding these 9 genes from the gene-restricted
   tier. Several of these (`NBPF1`, `SHANK3`, `POLR2A`) are known complex/segmentally-duplicated
   loci, plausibly why MANE represents them via an alt haplotype. **Not fixed here** — logged as a
   known gap alongside the 135 name-mismatch genes, same treatment: surfaced, not silently dropped,
   resolution (finding each gene's primary-chromosome coordinates) is out of this step's scope.

## Full unmatched-gene list (logged in `BUG_TRACKER.md` MISC-4, reproduced here in full)

Mostly mitochondrial genes (`MT-*` — the MANE GTF's `gene_name` field doesn't carry these the same
way as nuclear genes), older `C*ORFnn`-style aliases, a handful of literal Ensembl IDs present in
`Complete_gene_list_V5`'s `GENE` column instead of an HGNC symbol (`ENSG00000249624` etc.), and
small RNA genes (`SNORD*`, `MIR*`, `RNU*`):

```
C19ORF12, C9ORF72, COXFA4, DNAAF19, DRC2, DRC4, IGHM, LRTOMT, MT-ATP6, MT-CO1, MT-CO2, MT-CO3,
MT-CYB, MT-ND1, MT-ND3, MT-ND4, MT-ND5, MT-ND6, MT-RNR1, MT-TA, MT-TE, MT-TF, MT-TH, MT-TI, MT-TK,
MT-TL1, MT-TL2, MT-TM, MT-TN, MT-TP, MT-TS1, MT-TS2, MT-TV, MT-TW, MT-TY, MUC1, RNU2-2, RNU5B-1,
MIR96, MT-ND2, MT-TD, MT-TG, MT-TR, MT-TT, RNU7-1, TRAC, IGKC, MIR140, MT-ATP8, MT-ND4L, MT-RNR2,
MT-TC, MT-TQ, SNORA31, ARPIN-AP3S2, ASNA1, ATXN8, C10ORF71, C11ORF21, C11ORF58, C11ORF80, C11ORF91,
C12ORF57, C14ORF119, C14ORF180, C17ORF107, C18ORF32, C19ORF47, C1ORF105, C1ORF127, C1ORF131,
C20ORF27, C20ORF96, C21ORF140, C2ORF68, C2ORF69, C3ORF33, C3ORF49, C4ORF19, C4ORF46, C5ORF22,
C5ORF34, C6ORF163, C9ORF152, CEFIP, COMMD3-BMI1, CXORF38, CXORF58, DUX4L1, EEF1AKMT4-ECE2, ELAC,
ENSG00000249624, ENSG00000256349, ENSG00000260342, ENSG00000262660, ENSG00000271254,
ENSG00000278817, ENSG00000284906, ENSG00000285330, ENSG00000285827, ENSG00000286185,
ENSG00000289258, ENSG00000289697, FAM104A, FAM166B, FPGT-TNNI3K, FRA16E, GCOM1, HSPB11, HYMAI,
IL12A-AS1, IPW, LEXM, METTL7B, NBPF26, NDUFC2-KCTD14, NME1-NME2, NT5C1B-RDH14, NUTM2B-AS1, OAZ2,
PIFO, PKD1L2, PWAR1, PWRN1, RNF103-CHMP3, SAA2-SAA4, SLC9A3R2, SMIM4, SNORD115-1, SNORD116-1, TAZ,
TDGF1, TDGF1P3, TRU-TCA1-1, WISP1
```

Checked against `resources/gene_transcript_mapping.txt` (the 217-panel): only **`TAZ`** is an
actual panel gene in this unmatched list (confirmed present, twice — `NM_000116.4` cross-referenced
to two different ENST IDs, `ENST00000369776` and `ENST00000601016` — worth a separate look, unrelated
to this BED). `FPGT-TNNI3K` is a distinct readthrough-transcript gene symbol from `TNNI3K`, which
*is* a panel gene and matched fine on its own. Mitochondrial genes (`MT-*`) are lower concern for a
nuclear-splicing-predictor-restriction BED and aren't panel genes at all. `TAZ` is the one real gap
worth a manual coordinate override during item 7, not the others.

## What's next

Per the plan's revised migration sequencing: item 7 (whole-VCF entry point + generalized
`CHUNK_VCF`, reused for both tiers, + `GENE_SUBSET` using this new BED) is next.
