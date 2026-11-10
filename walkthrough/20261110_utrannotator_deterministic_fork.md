# UTRAnnotator: deterministic copy reporting every uORF class (PRED-15)

Kanvas SF-03. Mario chose option (a) on 2026-11-10: keep a copy of the plugin in this repo.

## Problem

FACT: the shared plugin (`/data_lab_PGP/resources/annotation/VEP/VEP_plugins/UTRAnnotator.pm`) loops over `keys %five_prime_flag` (line 262) and keeps the last flagged class. Perl randomises hash order per process, so a variant with several uORF effects gets one class picked at random. The order of the entries in `5UTR_annotation` (line 282) and of the key=value pairs inside each entry (line 283) also comes from hash order.

FACT, recomputed 2026-11-10 from the `F2005_union` and `F1830_union_avi` tables: 47 of 1,780 shared (`Locus`, `Feature`) rows, on 34 variants, report different classes in the two runs.

FACT, local, host Perl: the original selection loop applied to {gain, stop loss, frameshift} returned `stop_loss gain frameshift frameshift gain frameshift gain frameshift` under hash seeds 0–7.

## Change

- `src/external/vep_plugins/UTRAnnotatorPGP.pm` is a copy of the shared file. Source: Ensembl VEP_plugins, branch `release/116`, commit `867336c85c02a7842c87f484fd8d92cb44358992`, file md5 `fb9c8d5bc8ef2050d3b8303b76a10056`, unmodified in the shared checkout. The shared file is unchanged. The new `MODIFICATIONS` POD section records the changes, as the Apache-2.0 licence asks.
  - The package is renamed to `UTRAnnotatorPGP`. VEP puts `--dir_plugins` first in `@INC` (`Runner.pm:886` in the container), so a copy with the same name would be shadowed by the shared one.
  - `format_utr_effects` replaces the selection loop. `5UTR_consequence` lists every flagged class in a fixed order: uAUG gain, uAUG loss, uORF stop gain, uORF stop loss, uORF frameshift.
  - Each `5UTR_annotation` entry starts with `consequence=<SO term>`, followed by its key=value pairs sorted by key. Entries keep the plugin's numbering within each class.
  - VEP writes `,` as `&` in the VCF, so both fields come out `&`-separated in the table.
  - The JSON/REST branch uses the same order.
  - The field names (`5UTR_consequence`, `5UTR_annotation`, `Existing_*`) are unchanged.
- `modules/local/vep.nf` puts `params.vep_local_plugins_dir` on `PERL5LIB` inside the container and calls `--plugin UTRAnnotatorPGP`. `annotate_vcf.nf` adds that directory to `code_md5.vep`, so `-resume` reruns VEP when the plugin changes.
- `parse_utrannotator` (`src/python/parsing/modules/predictors.py`) grades an entry that carries `consequence=` only under its own class. Before, every entry was graded under every term. In a multi-class row that turned, for example, a uAUG-gain entry graded under "stop loss" into `unclassified`, because it has no `AltStop` field. Entries without `consequence=` (tables from the shared plugin) are graded as before.

## Effect on results

INFERENCE: on multi-effect rows `UTRAnnotator_status` now takes the worst class over all effects instead of the randomly chosen one. Rows can therefore become `high_impact`, which raises `NEW_IMPACT`. Single-effect rows only change in the annotation text: the `consequence=` prefix and the key order. Module 3 reads `5UTR_consequence` only as non-empty / empty (`clinical_variant_prioritization/src/modules/features.py:116`, `acmg.py:275`, `predictor_registry.py:65`), so a multi-valued value does not break it. It should still be told that the values changed.

## Checks (local, 2026-11-10)

- `tests/perl/run_in_vep_container.sh` (VEP container `vep.sif`): `perl -c` OK. `format_utr_effects` returns the same output on a synthetic gain + stop-loss + frameshift variant under `PERL_HASH_SEED` 0–7: 3 classes in order, 4 entries, the gain's 2 uORFs in plugin order.
- `tests/python/test_utrannotator_status.py`: the new linked-entry cases pass, and so does the old-format multi-term case.
- No VEP run locally: the VEP cache and FASTA are in `/references`, which is not mounted here.

## Cluster check (pending, Mario launches)

Input: `test_data/raw_vcfs/utr_multieffect_fixture.vcf.gz`, outside git, md5 `b71eb8b380fa500e6bb9608b6ff5d518`. It holds the 34 variants above, with synthetic trio genotypes (`0/1`, `0/0`, `0/1`).

How it was built, 2026-11-10:
- From `F2005_union` and `F1830_union_avi` (`tables/*.parsed.clean.pq`), read `Locus`, `Feature` and `5UTR_consequence` for the rows where `5UTR_consequence` is not empty or `-`.
- Join the two runs on (`Locus`, `Feature`) and keep the rows whose `5UTR_consequence` differs: 47 rows, 34 loci.
- Split each `Locus` (`chr:pos-ref-alt`) into a VCF record. The header is `clinvar_af_fixture.vcf.gz`'s, with the `##contig` lines of `F2005_union.annVEP.vcf.gz`.

FACT: all 34 records match `F2005_union.annVEP.vcf.gz` exactly on CHROM/POS/REF/ALT (tabix lookup), so the alleles are the ones that already passed the assembly check there.

Run both from the committed tree, without `-resume`, so that each VEP run is a new Perl process with its own hash seed:

```bash
bash run_annotate_vcf.sh -profile standard --input_vcf test_data/raw_vcfs/utr_multieffect_fixture.vcf.gz --annotation_cache_dir /mnt/scratch/mruizp/annotation_cache_test --chunk_size 5 --gene_restricted_chunk_size 2 --run_id utr5_multieffect_run1
bash run_annotate_vcf.sh -profile standard --input_vcf test_data/raw_vcfs/utr_multieffect_fixture.vcf.gz --annotation_cache_dir /mnt/scratch/mruizp/annotation_cache_test --chunk_size 5 --gene_restricted_chunk_size 2 --run_id utr5_multieffect_run2
```

Pass when all of these hold:

1. Both runs exit 0, and the VEP task log has no `Failed to compile plugin UTRAnnotatorPGP`. Without `--safe`, a plugin that fails only warns and its columns disappear.
2. On all 47 rows, `5UTR_consequence` is not empty or `-`, and `5UTR_annotation` starts with `consequence=`.
3. `5UTR_consequence`, `5UTR_annotation` and `UTRAnnotator_status` are identical between the two runs on every row.
4. On each of the 47 (`Locus`, `Feature`) rows, the new `5UTR_consequence` contains both classes the old runs reported (`F2005_union` and `F1830_union_avi`).

Criteria 2–4 in one command (local, after both runs; `annotation_tests` env):

```bash
O=nf_work/annotation_out
python src/tools/compare_utr5_runs.py \
  --run1 $O/utr5_multieffect_run1/tables/utr5_multieffect_run1.parsed.clean.pq \
  --run2 $O/utr5_multieffect_run2/tables/utr5_multieffect_run2.parsed.clean.pq \
  --old1 $O/F2005_union/tables/F2005_union.parsed.clean.pq \
  --old2 $O/F1830_union_avi/tables/F1830_union_avi.parsed.clean.pq
```

It prints `PASS` or `FAIL` and lists the rows that miss a class. It was tested on synthetic tables: a good case gives PASS, and a run that drops a class gives FAIL.
