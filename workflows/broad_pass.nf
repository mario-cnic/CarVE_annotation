// Phase 3: broad-safe tier — VEP (containerized) + Branchpointer, both run in parallel off the
// same pre-annotation per-gene VCF as GENE_RESTRICTED_SUBWORKFLOW's predictors, mirroring
// main.sh submitting VEP and Branchpointer as independent qsub jobs off the same $vcf_gz
// (main.sh:437-497). Still the legacy per-gene adapter only (Phase 7 adds the WGS/WES entry
// point) — single-gene-sized "regions" for now, per the plan's migration sequencing.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

include { VEP_ANNOTATE } from '../modules/local/vep.nf'
include { BRANCHPOINT_ANNOTATE } from '../modules/local/branchpointer.nf'
include { SPIP_ANNOTATE } from '../modules/local/spip.nf'

workflow BROAD_PASS_SUBWORKFLOW {
    take:
    meta_vcf_ch   // tuple(meta, vcf, tbi) — one element per gene

    main:
    VEP_ANNOTATE(meta_vcf_ch)
    BRANCHPOINT_ANNOTATE(meta_vcf_ch)

    // Mirrors GENE_RESTRICTED_SUBWORKFLOW's inverse gate: SPiP's tier is unresolved, defaulted to
    // 'restricted' (see that subworkflow), so this branch is a no-op today. Not invoked from both
    // subworkflows in the same run — each gate only fires for the tier params.spip_tier selects.
    if (params.spip_tier == 'broad') {
        SPIP_ANNOTATE(meta_vcf_ch)
        spip_out_ch = SPIP_ANNOTATE.out.vcf
    } else {
        spip_out_ch = Channel.empty()
    }

    emit:
    vep         = VEP_ANNOTATE.out.vcf
    branchpoint = BRANCHPOINT_ANNOTATE.out.vcf
    spip        = spip_out_ch
}
