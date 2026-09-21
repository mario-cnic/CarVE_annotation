// Phase 2: gene-restricted tier complete — SpliceAI (chunked, from Phase 1) + Pangolin + SPiP,
// all three fed from the same pre-annotation per-gene VCF, run in parallel (not chained) — mirrors
// main.sh submitting all three as independent qsub jobs off the same $vcf_gz (main.sh:420-467).
// Still the legacy per-gene adapter only (Phase 7 adds the WGS/WES entry point).
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

include { SPLICEAI_SPIKE } from './spliceai_spike.nf'
include { PANGOLIN_ANNOTATE } from '../modules/local/pangolin.nf'
include { SPIP_ANNOTATE } from '../modules/local/spip.nf'

workflow GENE_RESTRICTED_SUBWORKFLOW {
    take:
    meta_vcf_ch   // tuple(meta, vcf, tbi) — one element per gene

    main:
    SPLICEAI_SPIKE(meta_vcf_ch)
    PANGOLIN_ANNOTATE(meta_vcf_ch)

    // spip_tier default is 'restricted' (this subworkflow) per the user's explicit choice at plan
    // approval — SPiP's broad-safety is unresolved. A future 'broad' placement would instead wire
    // SPIP_ANNOTATE into BROAD_PASS_SUBWORKFLOW (not built yet) and this branch would go empty.
    if (params.spip_tier == 'restricted') {
        SPIP_ANNOTATE(meta_vcf_ch)
        spip_out_ch = SPIP_ANNOTATE.out.vcf
    } else {
        spip_out_ch = Channel.empty()
    }

    emit:
    spliceai = SPLICEAI_SPIKE.out.vcf
    pangolin = PANGOLIN_ANNOTATE.out.vcf
    spip     = spip_out_ch
}
