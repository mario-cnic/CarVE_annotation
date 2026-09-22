// Whole-VCF redesign (plan item 7-8): broad-safe tier for the new whole-VCF entry point. Chunks
// the FULL, unfiltered input (VEP/Branchpointer are genome-wide-capable, no region restriction
// needed), fans VEP_ANNOTATE + BRANCHPOINT_ANNOTATE [+ SPIP_ANNOTATE if params.spip_tier=='broad']
// out over the chunks, then reassembles each tool's chunks back into one whole-run file.
// Reuses the exact same processes as the legacy per-gene BROAD_PASS_SUBWORKFLOW
// (workflows/broad_pass.nf, untouched) — only the wiring around them differs.
//
// CONCAT_CHUNKS is imported under a distinct alias per tool: Nextflow disallows the same process
// being invoked from more than one call site within one workflow scope unless each site aliases
// it (see workflows/chunk_fan_concat.nf's header comment for how this was found).
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

include { CHUNK_VCF } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_VEP } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_BRANCHPOINT } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_SPIP_BROAD } from '../modules/local/chunk.nf'
include { VEP_ANNOTATE } from '../modules/local/vep.nf'
include { BRANCHPOINT_ANNOTATE } from '../modules/local/branchpointer.nf'
include { SPIP_ANNOTATE } from '../modules/local/spip.nf'

// Groups a tool's per-chunk output back by run and pairs it with the pre-chunk input for
// CONCAT_CHUNKS's parity check. No sort-by-filename here (an earlier version had one, following
// Phase 1's spliceai_spike.nf pattern, on the theory that it fixed groupTuple()'s unspecified
// element order) — removed after BUG_TRACKER.md MISC-8 showed it was a no-op for every predictor
// except SpliceAI (whose per-chunk filename is uniquely named; VEP/Branchpointer/Pangolin/SPiP's
// aren't) and, more to the point, unnecessary regardless: CONCAT_CHUNKS's own `bcftools sort`
// (modules/local/chunk.nf) now guarantees correct final ordering no matter what order chunk_files
// arrives in here.
def collectForConcat(annotated_chunk_ch, orig_vcf_ch, suffix) {
    annotated_chunk_ch
        .groupTuple()
        .join(orig_vcf_ch)
        .map { meta, chunk_files, vcf, tbi -> tuple(meta, chunk_files, vcf, tbi, suffix) }
}

workflow BROAD_PASS_WHOLE_SUBWORKFLOW {
    take:
    meta_vcf_ch   // tuple(meta, vcf, tbi) — the FULL, unfiltered whole-run input

    main:
    CHUNK_VCF(meta_vcf_ch)

    // Fan out directly off CHUNK_VCF.out.chunks, same as Phase 1's spliceai_spike.nf — NOT
    // wrapped in a shared FAN_CHUNKS subworkflow (an earlier version of this file did; dropped
    // after empirically narrowing a hard-to-reproduce bug to that indirection, see the plan's
    // walkthrough for the diagnostic trail).
    chunk_ch = CHUNK_VCF.out.chunks
        .map { meta, vcf_list, tbi_list ->
            tuple(
                meta,
                vcf_list instanceof List ? vcf_list : [vcf_list],
                tbi_list instanceof List ? tbi_list : [tbi_list]
            )
        }
        .flatMap { meta, vcf_list, tbi_list ->
            [vcf_list, tbi_list].transpose().collect { vcf, tbi -> tuple(meta, vcf, tbi) }
        }

    VEP_ANNOTATE(chunk_ch)
    CONCAT_VEP(collectForConcat(VEP_ANNOTATE.out.vcf.map { m, vcf, tbi -> tuple(m, vcf) }, meta_vcf_ch, 'annVEP'))

    BRANCHPOINT_ANNOTATE(chunk_ch)
    CONCAT_BRANCHPOINT(collectForConcat(BRANCHPOINT_ANNOTATE.out.vcf.map { m, vcf, tbi -> tuple(m, vcf) }, meta_vcf_ch, 'annBranchpoint'))

    // Mirrors the legacy BROAD_PASS_SUBWORKFLOW's inverse gate: SPiP's tier is unresolved,
    // defaulted to 'restricted' (see the whole-VCF GENE_RESTRICTED_WHOLE_SUBWORKFLOW), so this
    // branch is a no-op today.
    if (params.spip_tier == 'broad') {
        SPIP_ANNOTATE(chunk_ch)
        CONCAT_SPIP_BROAD(collectForConcat(SPIP_ANNOTATE.out.vcf.map { m, vcf, tbi -> tuple(m, vcf) }, meta_vcf_ch, 'annSPiP'))
        spip_out_ch = CONCAT_SPIP_BROAD.out.vcf
    } else {
        spip_out_ch = Channel.empty()
    }

    emit:
    vep         = CONCAT_VEP.out.vcf
    branchpoint = CONCAT_BRANCHPOINT.out.vcf
    spip        = spip_out_ch
}
