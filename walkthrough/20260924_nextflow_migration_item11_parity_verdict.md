# Nextflow migration — item 11 parity verdict (first real cluster run of `panel7_test.vcf.gz`)

**2026-09-24**

## Context

`bash run_annotate_vcf.sh -profile standard --input_vcf test_data/raw_vcfs/panel7_test.vcf.gz
--chunk_size 15` completed on the real SGE cluster (run name `suspicious_curie`, 2026-09-22
15:20-19:43, 4h22m). This entry records the parity verification against the known-good bash-`main.sh`
reference and against the plan's own item-11 checklist (`scalable-wibbling-snowflake.md`).

## 1. Record-count parity (all 7 genes)

All 5 predictor outputs + the merged VCF: **72/72** records, matching the constructed input exactly
(`ACTC1: 6, DSP: 9, KCNQ1: 10, MYBPC3: 22, PKP2: 10, SCN5A: 6, TNNT2: 9`). Table shape: 803 rows ×
577 columns (one row per variant×transcript, no `--gene_set` filter, as designed).

## 2. Content parity — MYBPC3 only (the only gene with a bash-`main.sh` baseline)

Raw `diff` of the 22 MYBPC3 records against `test_data/test_run/annotation/MYBPC3.annotated.vcf.gz`
showed all 22 records differing. Investigation (transcript count, gene symbols, sorted-CSQ-block
diff) isolated the entire remaining diff to **SPiP-field numeric formatting only** — same values,
different text representation (`10` vs `10.0000000`, leading whitespace padding, and, for small
probabilities, differing significant-figure counts in scientific notation, e.g. `7.8731384e-07` vs
`7.873138e-07`). Confirmed same command (`Rscript SPiPv2.1_main.r`, unchanged by the port) —
attributed to `spip_env`/R version drift between whenever the reference file was generated and now,
not a port defect.

**Normalization rule used** (records this so the claim is reproducible, not asserted as byte
parity): strip padding whitespace around `|` delimiters, then round every decimal number — including
scientific notation — to 6 significant figures on *both* sides before diffing:
```python
re.sub(r'-?\d+\.\d+(?:[eE][+-]?\d+)?', lambda m: f'{float(m.group(0)):.6g}' if float(m.group(0)) != int(float(m.group(0))) else str(int(float(m.group(0)))), line)
```
Result: **0/22 variants show any diff after normalization.**

**Verdict language, precisely**: MYBPC3's CSQ+INFO content is *identical to within 6 significant
figures* between the Nextflow cluster run and the bash reference — not byte-identical (unlike
Phases 3/4's Branchpointer and merge-subworkflow checks, which were byte-identical against their
references). This is a looser grade of evidence, and is recorded as such rather than overstated.

## 3. New-functionality checks (not covered by #1/#2 — these test the *new* code, not ported code)

**`GENE_SUBSET` behavior**: inspected the real work directory
(`work/a8/b9e960f14f1ef34ae122f1863a9903/subset.vcf.gz`) directly — 72/72 records, i.e. the subset
step removed **zero** variants. All 7 genes, including `KCNQ1` (deliberately built outside the
217-panel to test tier-2/3 fallback), are apparently within `Complete_gene_list_V5`'s ~5,877-gene
coverage. **This means `GENE_SUBSET`'s actual filtering behavior was not exercised by this test
run** — logged as a coverage gap, not a bug (see "What's not yet covered" below).

**`TRANSCRIPT_PRIORITY_TIER` tagging, on real cluster output** (not just the pre-run standalone
check against MYBPC3): read `panel7_test.parsed.clean.pq` directly.

| SYMBOL | Tier 1 | Tier 2 | Tier 3 |
|---|---|---|---|
| ACTC1 | 6 | 0 | 30 |
| DSP | 9 | 0 | 39 |
| GJD2-DT | 0 | 0 | 27 |
| KCNQ1 | **0** | 10 | 40 |
| KCNQ1OT1 | 0 | 0 | 8 |
| MADD | 0 | 9 | 167 |
| MYBPC3 | 17 | 0 | 28 |
| PKP2 | 10 | 0 | 78 |
| SCN5A | 6 | 0 | 47 |
| TNNT2 | 9 | 0 | 230 |

Exactly the expected discriminator: **every 217-panel gene** (ACTC1, DSP, MYBPC3, PKP2, SCN5A,
TNNT2) has nonzero tier-1 rows; **KCNQ1 — deliberately outside the panel — has zero**, falling back
to tier-2 (`MANE_SELECT`). `MADD` (a real overlapping neighbor gene near MYBPC3, not itself
synthetic-test-targeted) also correctly has no tier-1 rows and falls to tier-2/3. Confirms the
tiering logic is correct on real, non-synthetic-standalone cluster output, not just the earlier
isolated pre-run check.

**Provenance manifest**: `RUN_MANIFEST.json`'s `completed` stage has `exit_code: 0`;
`launch.resolved_paths.nextflow_binary` = `/opt/nextflow/nextflow` (the cluster-specific path found
on the first real-cluster attempt, correctly resolved and recorded — not the earlier bug's silent
`null`); `launch.nextflow_version` = `nextflow version 25.10.2.10555`; repo commit `24ca486`, clean
tree. All fields this session's two infrastructure-only cluster fixes were specifically about are
now confirmed populated from an actual real run, not just inferred safe.

## What's not yet covered (open, flagged for the user, not silently decided)

Per the plan's own item-11 text: "real multi-gene input (not just MYBPC3)... spot-checked content
against the equivalent `bash main.sh` outputs." Record-count and table-shape evidence now covers
all 7 genes; **content parity covers MYBPC3 only** — it is the only gene with an existing bash-run
baseline. The other 6 genes (`ACTC1`, `TNNT2`, `PKP2`, `DSP`, `SCN5A`, `KCNQ1`) are new synthetic
data with no prior bash-`main.sh` run to diff against. Closing this gap needs a real `bash main.sh`
run on the cluster for those 6 genes — real cluster time, not something to spend without asking.

`GENE_SUBSET`'s actual record-removal behavior is untested by this run (see above) — would need a
variant on a gene confirmed absent from `Complete_gene_list_V5` to exercise it.

## Verdict

Item 11 is **partially satisfied**: full record-count/table-shape parity across all 7 genes, full
content parity (within float-formatting tolerance) for MYBPC3, and both pieces of genuinely new
logic (`GENE_SUBSET`, `TRANSCRIPT_PRIORITY_TIER`) verified correct on real cluster output where
testable. Not yet closed: content parity for the other 6 genes, and a real test of `GENE_SUBSET`
actually excluding something. Left to the user to decide whether this is sufficient to proceed to
WGS-scale work now, or whether to spend more cluster time closing the remaining gap first.
