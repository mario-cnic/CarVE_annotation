// Joins the five predictor channels on meta and merges them into one annotated VCF.

include { MERGE_ANNOTATIONS } from '../modules/local/merge.nf'

workflow MERGE_SUBWORKFLOW {
    take:
    vep_ch         // tuple(meta, vcf, tbi)
    branchpoint_ch // tuple(meta, vcf, tbi)
    pangolin_ch    // tuple(meta, vcf, tbi)
    spliceai_ch    // tuple(meta, vcf, tbi)
    spip_ch        // tuple(meta, vcf, tbi) — already resolved to whichever tier is active

    main:
    merge_input_ch = vep_ch
        .join(branchpoint_ch)
        .join(pangolin_ch)
        .join(spliceai_ch)
        .join(spip_ch)

    MERGE_ANNOTATIONS(merge_input_ch)

    emit:
    vcf = MERGE_ANNOTATIONS.out.vcf
}
