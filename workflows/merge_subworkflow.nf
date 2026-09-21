// Phase 4: MERGE_SUBWORKFLOW — joins the two tiers' five predictor output channels by `meta` into
// one tuple and hands it to MERGE_ANNOTATIONS. Still the trivial 1:1 single-gene case (Phase 4
// scope per the plan); the overlapping-gene disambiguation the plan describes for this subworkflow
// is Phase 7 work, once multi-gene "regions" exist.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

include { MERGE_ANNOTATIONS } from '../modules/local/merge.nf'

workflow MERGE_SUBWORKFLOW {
    take:
    vep_ch         // tuple(meta, vcf, tbi)
    branchpoint_ch // tuple(meta, vcf, tbi)
    pangolin_ch    // tuple(meta, vcf, tbi)
    spliceai_ch    // tuple(meta, vcf, tbi)
    spip_ch        // tuple(meta, vcf, tbi) — already resolved to whichever tier is active

    main:
    // .join() keys on the first tuple element (meta) and folds it into a single copy per join,
    // same operator Phase 1's spliceai_spike.nf already relies on for meta-keyed joins.
    merge_input_ch = vep_ch
        .join(branchpoint_ch)
        .join(pangolin_ch)
        .join(spliceai_ch)
        .join(spip_ch)

    MERGE_ANNOTATIONS(merge_input_ch)

    emit:
    vcf = MERGE_ANNOTATIONS.out.vcf
}
