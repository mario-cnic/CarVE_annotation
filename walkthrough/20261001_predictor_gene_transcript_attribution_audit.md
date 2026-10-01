# Predictor score attribution: wrong-gene / wrong-transcript scores on overlapping loci

**2026-10-01.** Trigger: the Pangolin re-test on `panel7` returned `MADD` scores for variants the fixture labels `MYBPC3`.
Scope: read-only audit of the real `S223` table + alias-resolution feasibility. No pipeline behaviour changed.

## FACTS

**Why the label said MYBPC3.** `panel7_test.vcf.gz` is synthetic; its `gene=` tag was set by the generator from the padded window
(`v5_genes_loc.bed` = gene body ±10 kb). GENCODE 45: MYBPC3 chr11:47,331,406-47,352,702 (− strand); the 8 variants at
47,328,406-47,330,000 are 1.4-3.0 kb **downstream** (3') of it and inside MADD (47,269,161-47,330,031, + strand). Pangolin correctly scored MADD.

**How scores are attached (shared parser, vendored copy `src/python/parsing/modules/splicing.py`):**
| Predictor | Unit | Parser behaviour |
|---|---|---|
| SpliceAI (`parse_spliceai_custom`) | gene | takes only the **first** entry; never matches row gene |
| Pangolin (`parse_pangolin`) | gene (ENSG) | max \|score\| over **all** genes; gene id never read |
| SPiP (`parse_spip`) | transcript (`NM_`) | first entry whose gene symbol == row SYMBOL, else **first entry of any gene**; transcript never matched |
Rows are variant x VEP transcript; predictor strings are variant-level, so every row of a variant gets the same string.

**Audit** (`src/tools/audit_gene_attribution.py`, S223, 8,993,048 populated rows; rows, not distinct variants; not filtered to curated transcripts):
- SpliceAI (6,233,537 rows): first entry = row gene 85.2%; wrong entry used though the right one was present 81,564 (270 cross the 0.20 floor); row gene not scored, neighbour score attached 576,599 (3,553 with score >= 0.20); row has no symbol 261,471.
- SPiP (7,756,696 rows): row gene absent -> other gene's first entry used 1,237,426 (4,042 >= 0.20); gene has several SPiP transcripts 4,878,242; of 470,073 multi-transcript rows with a MANE row, 301,855 take a non-MANE transcript (136 cross 0.20); MANE not among SPiP entries 34,281.
- Pangolin: S223 predates the db fix (all 2,791,464 rows `ENSMUSG`), so gene matching is not assessable; 128,206 rows multi-gene, 582 of them with the 0.20 outcome depending on which gene.

**Alias vs true neighbour** (`src/tools/audit_alias_resolution.py`; 576,599 SpliceAI "row gene not scored" rows). Resolve each SpliceAI symbol to ENSG
by exon overlap with GENCODE 45 (same chrom+strand), no alias list: alias 187,734 rows (32.6%, 372 distinct pairs), true neighbour 388,434 (67.4%), unresolved 431 (3 symbols).
Spot check of the alias class: PARK2->PRKN, SEPT9->SEPTIN9, IKBKAP->ELP1, FAM134B->RETREG1. Neighbours are mostly overlapping antisense genes (HRG-AS1/KNG1) and include PMS2 rows carrying AIMP2's score.
Cross-check with HGNC complete set (release file dated 2026-09-29, sha256 `91d0ad20...9d889e`; current/prev/alias symbol -> ENSG): coordinate and HGNC classes agree on 10,488 of 11,043 distinct pairs (95.0%; 94.1% of rows);
HGNC has no record for 535 pairs (5.7% of rows; coordinates decide); they **disagree on 19 pairs (1,178 rows, 0.2%)**; coordinates leave 3 pairs unresolved, HGNC resolves 1 of them. Combined: 2 pairs (3 rows) unresolved.

## INFERENCE / limits
- Part of the "neighbour" rows may be alias cases the two methods both miss (HGNC-less symbols judged only by coordinates).
- The 0.20 flips are the rows that can create a false report flag; how many reach tier-1 transcripts (Module 3's input) was not measured.
- One sample (S223). SPiP alias resolution not tested (no SPiP annotation locally; `NM_`->ENSG via HGNC `refseq_accession` + MANE GTF is the proposed route).

## Module 3 relationship (checked)
Module 3 has its **own** copies (`clinical_variant_prioritization/src/filter_variants.py`, `src/modules/splicing.py`), already diverged from the shared ones (shared has `parse_maxentscan`).
Earlier statement in this session that the shared module is "used by Module 3" was wrong. `DEC-0003`: Module 2 owns parsing to the table, Module 3 owns interpretation.

## Actions taken
- Vendored the parsing code into `src/python/parsing/` (commit d5a915c, verbatim, sha256-verified, provenance in its PROVENANCE.md); Nextflow config repointed (fa87d38). Legacy bash (`config/env.sh`, `vcf2parsed.sh`) not repointed.
- HGNC complete set downloaded locally to `~/hgnc_build/` (not deployed anywhere).
- Next: `TODO.md` item "Predictor score attribution".
