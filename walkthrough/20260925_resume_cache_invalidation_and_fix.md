# `-resume` cache invalidation cascade on S223's relaunch, and the fix (`CHUNK_VCF_RESTRICTED`)

**2026-09-25**

## Correction of a wrong claim made earlier this session

Before this relaunch, I told the user `-resume` would likely still cache-hit the broad tier's
255 already-succeeded `VEP_ANNOTATE`/`BRANCHPOINT_ANNOTATE` tasks despite `CHUNK_VCF` itself
re-executing (its process definition changed for `MISC-11`), reasoning: "Nextflow caches
downstream tasks by file content, not by which upstream task produced them... the re-chunk should
produce byte-identical output... VEP/Branchpointer should then correctly cache-hit."

**This was wrong**, and it cost a real relaunch attempt to find out. Nextflow's *default* cache
mode (`cache 'standard'`) keys a `path` input's identity on **path + size + last-modified
timestamp**, not a content checksum — content-addressed hashing only happens under the explicit
`cache 'deep'` mode, which this pipeline does not use. When `CHUNK_VCF` re-executes (for any
reason — its own definition changed, its script text changed, whatever), it writes its output
chunks to a brand-new work directory, i.e. brand-new paths, every time. Every downstream task
consuming those chunks then sees a "new" input by Nextflow's own bookkeeping, regardless of
whether the bytes are identical to a previous run's chunks. Confirmed for real on this relaunch:
`executor > sge (103)` and climbing, with `VEP_ANNOTATE (S223) 0 of 255` — a full, unwanted
re-submission of the entire broad tier to SGE, not a cache hit.

## Root cause, precisely

`modules/local/chunk.nf`'s `CHUNK_VCF` process definition changed for `MISC-11` (a `val(chunk_size)`
input was added, and the script line changed from `${params.chunk_size}` to `${chunk_size}`).
Nextflow's task hash folds in each declared input's resolved value, not just the rendered script
text — so even though the *broad tier's* actual `chunk_size` value (`20000`) never changed, adding
a new declared input to the process changed its hash for **every** call site, including the one
whose behavior was otherwise 100% unchanged. This, independently, combined with the path/mtime-based
(not content-based) caching above to guarantee the entire broad tier would re-run regardless of
`-resume`.

Verified via advisor before acting further (given the cost of guessing wrong a second time):
`git diff HEAD -- modules/local/chunk.nf workflows/broad_pass_whole.nf
workflows/gene_restricted_whole.nf` against the exact commit (`eb5f13f`) `spontaneous_saha` ran
against, confirming the diff was minimal enough for a surgical revert to be viable, and confirmed
directly (not inferred) that the original broad-tier `CHUNK_VCF` work dir
(`work/bf/6aff974cd9e37edc1c362080eac87e/chunks/`) and its 255 downstream VEP/Branchpointer results,
plus the `spontaneous_saha` session's `.nextflow/cache/f2fbe04a-.../db`, all still physically exist
on disk.

## Fix: `CHUNK_VCF_RESTRICTED`, a deliberate duplicate, not a shared parameterized process

- `modules/local/chunk.nf`'s `CHUNK_VCF` reverted to **byte-identical** to the commit
  `spontaneous_saha` ran against (confirmed via `git diff HEAD` returning empty for
  `workflows/broad_pass_whole.nf`, and `CHUNK_VCF`'s own script/input block matching the original).
- New process `CHUNK_VCF_RESTRICTED` added to the same file — same script, but with the
  parameterized `val(chunk_size)` input `MISC-11` actually needs. Used **only** by
  `workflows/gene_restricted_whole.nf`.
- `broad_pass_whole.nf` calls `CHUNK_VCF(meta_vcf_ch)` exactly as before — should now correctly
  reuse the same task hash as `spontaneous_saha`'s successful run, letting `-resume` reuse both
  `CHUNK_VCF` itself and (assuming its own inputs are otherwise unchanged) the 255 already-paid-for
  VEP/Branchpointer results.
- `gene_restricted_whole.nf` calls `CHUNK_VCF_RESTRICTED(subset_ch, params.gene_restricted_chunk_size)`
  — this tier was never going to cache-hit regardless (its own upstream, `GENE_SUBSET`, has a
  genuinely-changed `v5_genes_loc.bed` input this run, `MISC-10`), so there's no resume cost to
  giving it a distinct, purpose-built process.

Deliberate duplication over a shared parameterized process, documented as such in
`modules/local/chunk.nf`'s own header comment on `CHUNK_VCF_RESTRICTED` — the DRY alternative is
exactly what broke resume the first time.

## Verification

- `git diff HEAD -- workflows/broad_pass_whole.nf` returns empty — confirmed byte-identical to the
  commit that produced the cached results we're trying to preserve.
- Full `-stub-run` (isolated scratch work dir, not the shared repo's own `work/`) completes cleanly
  end-to-end with both `CHUNK_VCF` and `CHUNK_VCF_RESTRICTED` correctly wired to their respective
  tiers.
- **Not yet verified**: whether the actual cache-hit happens on the real cluster relaunch. That's
  the whole point of trying again — this is a real, testable prediction, not a repeat of the same
  unverified claim as before. If it's wrong again, the next step is to stop reasoning about
  Nextflow's cache semantics from first principles and instead directly inspect
  `.nextflow/cache/<session>/db` or add `-with-trace`/`-dump-hashes`-style diagnostics before
  trying a third code change.

## Outcome (same day): the fix didn't work, reverted

Relaunched with this fix in place — still no cache-hit, same `executor > sge (N)` climbing pattern
for the broad tier as before. Root cause of *why the fix itself didn't work* was not chased further
(would need direct `.nextflow/cache` inspection or `-dump-hashes`, not more first-principles
reasoning about Nextflow's cache internals, which has now been wrong twice on this exact question).

Mario's call: revert the `CHUNK_VCF`/`CHUNK_VCF_RESTRICTED` split. It added a second,
near-duplicate process for a resume-preservation benefit that, empirically, wasn't materializing —
pure added complexity with no realized payoff. Given the user has stated repeatedly that a full
reprocess of both tiers on every relaunch is acceptable, there's no reason to keep chasing this.
Back to the single shared `CHUNK_VCF(vcf_tuple, chunk_size)` design from before this file started
(`modules/local/chunk.nf`, `workflows/{broad_pass_whole,gene_restricted_whole}.nf`) — verified with
a clean `-stub-run`. `-resume` on this pipeline should now be assumed **not** to preserve the broad
tier's results across a relaunch whenever `CHUNK_VCF`'s shared script changes for any reason, and
that's an accepted tradeoff, not an open bug to revisit.
