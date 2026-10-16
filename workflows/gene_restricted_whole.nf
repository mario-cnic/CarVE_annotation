// Gene-restricted tier: subsets the input to the gene BED once, chunks it, runs SpliceAI, Pangolin
// (and SPiP when spip_tier == 'restricted', the default) per chunk, and reassembles each predictor's
// chunks. CONCAT_CHUNKS needs one alias per call site.

include { GENE_SUBSET } from '../modules/local/gene_subset.nf'
include { CHUNK_VCF } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_SPLICEAI } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_PANGOLIN } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_SPIP_RESTRICTED } from '../modules/local/chunk.nf'
include { SPLICEAI_ANNOTATE } from '../modules/local/spliceai.nf'
include { PANGOLIN_ANNOTATE } from '../modules/local/pangolin.nf'
include { SPIP_ANNOTATE } from '../modules/local/spip.nf'
include { FAN_CHUNKS } from './chunk_fan_concat.nf'

// Groups a predictor's chunks by run and pairs them with the pre-chunk input for the parity check.
def collectForConcat(annotated_chunk_ch, orig_vcf_ch, suffix) {
    annotated_chunk_ch
        // Hash order, not completion order: the same chunk list on every launch, so -resume keeps CONCAT_* cached.
        .groupTuple(sort: 'hash')
        .join(orig_vcf_ch)
        .map { meta, chunk_files, vcf, tbi -> tuple(meta, chunk_files, vcf, tbi, suffix) }
}

workflow GENE_RESTRICTED_WHOLE_SUBWORKFLOW {
    take:
    meta_vcf_ch   // tuple(meta, vcf, tbi) — the FULL, unfiltered whole-run input

    main:
    GENE_SUBSET(meta_vcf_ch, file(params.gene_restriction_bed))
    subset_ch = GENE_SUBSET.out.vcf   // tuple(meta, subset.vcf.gz, subset.vcf.gz.tbi)

    CHUNK_VCF(subset_ch, params.gene_restricted_chunk_size)
    chunk_ch = FAN_CHUNKS(CHUNK_VCF.out.chunks)

    // SPLICEAI_ANNOTATE takes (meta, vcf) without the index.
    SPLICEAI_ANNOTATE(chunk_ch.map { m, vcf, tbi -> tuple(m, vcf) })
    CONCAT_SPLICEAI(collectForConcat(SPLICEAI_ANNOTATE.out.annotated_chunk, subset_ch, 'annSpliceAI'))

    PANGOLIN_ANNOTATE(chunk_ch)
    CONCAT_PANGOLIN(collectForConcat(PANGOLIN_ANNOTATE.out.vcf.map { m, vcf, tbi -> tuple(m, vcf) }, subset_ch, 'annPangolin'))

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
