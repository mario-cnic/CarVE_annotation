// Broad tier: chunks the full input, runs VEP, Branchpointer (and SPiP when spip_tier == 'broad')
// per chunk, and reassembles each predictor's chunks. CONCAT_CHUNKS needs one alias per call site.

include { CHUNK_VCF } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_VEP } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_BRANCHPOINT } from '../modules/local/chunk.nf'
include { CONCAT_CHUNKS as CONCAT_SPIP_BROAD } from '../modules/local/chunk.nf'
include { VEP_ANNOTATE } from '../modules/local/vep.nf'
include { BRANCHPOINT_ANNOTATE } from '../modules/local/branchpointer.nf'
include { SPIP_ANNOTATE } from '../modules/local/spip.nf'

// Groups a predictor's chunks by run and pairs them with the pre-chunk input for the parity check.
def collectForConcat(annotated_chunk_ch, orig_vcf_ch, suffix) {
    annotated_chunk_ch
        // Hash order, not completion order: the same chunk list on every launch, so -resume keeps CONCAT_* cached.
        .groupTuple(sort: 'hash')
        .join(orig_vcf_ch)
        .map { meta, chunk_files, vcf, tbi -> tuple(meta, chunk_files, vcf, tbi, suffix) }
}

workflow BROAD_PASS_WHOLE_SUBWORKFLOW {
    take:
    meta_vcf_ch   // tuple(meta, vcf, tbi) — the FULL, unfiltered whole-run input

    main:
    CHUNK_VCF(meta_vcf_ch, params.chunk_size)

    // One channel element per chunk; a single chunk arrives as a bare Path, so wrap it in a list.
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
