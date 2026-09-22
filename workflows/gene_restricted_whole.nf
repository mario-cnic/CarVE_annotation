// Whole-VCF redesign (plan item 7-8): gene-restricted tier for the new whole-VCF entry point.
// Subsets the FULL input to Complete_gene_list_V5's regions ONCE (GENE_SUBSET), chunks that
// subsetted VCF, then fans SpliceAI/Pangolin [+ SPiP if params.spip_tier=='restricted', DEFAULT]
// out over the chunks and reassembles each tool's chunks back into one whole-run file.
//
// SpliceAI here runs directly on CHUNK_VCF's output via SPLICEAI_ANNOTATE, same process the
// legacy path's SPLICEAI_SPIKE already uses — Phase 1's own inner chunk/concat (SPLICEAI_CHUNK/
// SPLICEAI_CONCAT) is NOT reused here, since it's now redundant with this subworkflow's outer
// chunking (see the plan's "Why one chunker, reused twice" section: chunk boundaries don't affect
// SpliceAI's correctness, each variant is scored independently against the FASTA).
//
// CONCAT_CHUNKS is imported under a distinct alias per tool — see
// workflows/broad_pass_whole.nf's header comment for why (Nextflow disallows the same process
// being invoked from more than one call site within one workflow scope without aliasing).
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

include { GENE_SUBSET } from '../modules/local/gene_subset.nf'
include { CHUNK_VCF } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_SPLICEAI } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_PANGOLIN } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_SPIP_RESTRICTED } from '../modules/local/chunk.nf'
include { SPLICEAI_ANNOTATE } from '../modules/local/spliceai.nf'
include { PANGOLIN_ANNOTATE } from '../modules/local/pangolin.nf'
include { SPIP_ANNOTATE } from '../modules/local/spip.nf'
include { FAN_CHUNKS } from './chunk_fan_concat.nf'

// Same inline collect-for-concat sequence as broad_pass_whole.nf — see that file's header
// comment for why this isn't a shared subworkflow, and for why there's no sort-by-filename here
// (removed per BUG_TRACKER.md MISC-8: it was a no-op for Pangolin/SPiP's identically-named
// per-chunk output, and unnecessary anyway now that CONCAT_CHUNKS's own `bcftools sort` guarantees
// correct final ordering regardless of input order).
def collectForConcat(annotated_chunk_ch, orig_vcf_ch, suffix) {
    annotated_chunk_ch
        .groupTuple()
        .join(orig_vcf_ch)
        .map { meta, chunk_files, vcf, tbi -> tuple(meta, chunk_files, vcf, tbi, suffix) }
}

workflow GENE_RESTRICTED_WHOLE_SUBWORKFLOW {
    take:
    meta_vcf_ch   // tuple(meta, vcf, tbi) — the FULL, unfiltered whole-run input

    main:
    GENE_SUBSET(meta_vcf_ch)
    subset_ch = GENE_SUBSET.out.vcf   // tuple(meta, subset.vcf.gz, subset.vcf.gz.tbi)

    CHUNK_VCF(subset_ch)
    chunk_ch = FAN_CHUNKS(CHUNK_VCF.out.chunks)

    // SPLICEAI_ANNOTATE takes/emits tuple(meta, chunk) — 2 elements, no tbi at all — unlike the
    // other four predictors, which all expect the legacy per-gene path's 3-tuple shape (it was
    // built chunk-shaped from Phase 1 onward, when chunks never carried a tbi either — FAN_CHUNKS
    // added tbi propagation later, for the other four's benefit, not SpliceAI's). Drop it here.
    SPLICEAI_ANNOTATE(chunk_ch.map { m, vcf, tbi -> tuple(m, vcf) })
    CONCAT_SPLICEAI(collectForConcat(SPLICEAI_ANNOTATE.out.annotated_chunk, subset_ch, 'annSpliceAI'))

    PANGOLIN_ANNOTATE(chunk_ch)
    CONCAT_PANGOLIN(collectForConcat(PANGOLIN_ANNOTATE.out.vcf.map { m, vcf, tbi -> tuple(m, vcf) }, subset_ch, 'annPangolin'))

    // Mirrors the legacy GENE_RESTRICTED_SUBWORKFLOW's gate: default tier, per the user's
    // explicit choice at plan approval — SPiP's broad-safety is unresolved.
    if (params.spip_tier == 'restricted') {
        SPIP_ANNOTATE(chunk_ch)
        CONCAT_SPIP_RESTRICTED(collectForConcat(SPIP_ANNOTATE.out.vcf.map { m, vcf, tbi -> tuple(m, vcf) }, subset_ch, 'annSPiP'))
        spip_out_ch = CONCAT_SPIP_RESTRICTED.out.vcf
    } else {
        spip_out_ch = Channel.empty()
    }

    emit:
    spliceai = CONCAT_SPLICEAI.out.vcf
    pangolin = CONCAT_PANGOLIN.out.vcf
    spip     = spip_out_ch
}
