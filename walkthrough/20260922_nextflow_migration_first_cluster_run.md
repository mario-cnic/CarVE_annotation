# Nextflow migration — first real cluster run (item 11)

**2026-09-22, continuation of the same session**

## Context

First actual execution of `annotate_vcf.nf` on the real SGE cluster
(`bash run_annotate_vcf.sh -profile standard --input_vcf test_data/raw_vcfs/panel7_test.vcf.gz
--chunk_size 15`), run by the user (this sandbox cannot execute `qsub`/`qstat`). Two
infrastructure-only bugs (`nextflow` not on `PATH`, `write_run_manifest.py`'s Python 3.10+ syntax)
were already found and fixed on the first two attempts — see this file's neighbor,
`20260922_nextflow_migration_multi_gene_test_data.md`'s addendum. This entry covers the third
attempt's finding, the first one that's actually inside the pipeline's own logic.

## Real bug found: `VEP_ANNOTATE`'s `PERL5OPT=-X` (`BUG_TRACKER.md` `ORCH-11`)

Every `VEP_ANNOTATE` chunk failed identically:

```
Command error:
  Illegal switch in PERL5OPT: -X.
```

**Root cause, confirmed**: `-X` is a legal `perl` *command-line* switch (disables all warnings),
but `PERL5OPT` — the environment-variable route for passing switches — has its own, separate,
security-motivated allowed-switch list that does not include `-X`. Perl refuses it outright.

**Why this was never caught before, and why it's a real production finding, not just a porting
bug**: `src/hpc/annotate_vep_vars.sh:33` sets `export PERL5OPT=-X` as a plain host environment
variable, then calls `singularity exec -e ...` (`-e` is `--cleanenv`). `--cleanenv` strips the
host environment from the container — a plain `export` never reaches anything inside it; only
`SINGULARITYENV_`-prefixed variables do. So in the bash pipeline, this line has **always been a
silent no-op** — it has never actually suppressed the Mastermind plugin warnings it was written
for, in any real run, ever. Nextflow's own singularity integration propagates env vars exported
within a process's script into the container more faithfully than a manual bash `singularity exec`
call does — which is exactly why a bug that's been dormant in production surfaces as a hard
failure here instead.

**Fix**: `modules/local/vep.nf` simply drops the line — not replaced with a working equivalent,
since the original problem (disk-I/O noise from warning volume) is a performance/logging concern,
not a correctness one, and out of scope for this port. `main.sh` itself is untouched and still
carries the dormant, ineffective line (`ORCH-11`, left open — whether the original warning-flood
problem is actually happening in production, unnoticed, is now a real open question worth
someone checking, not something this session resolved).

## What's next

Fix pushed; next cluster attempt pending. This walkthrough will be updated (or a new dated entry
added) with the outcome once the run actually gets past `VEP_ANNOTATE`.
