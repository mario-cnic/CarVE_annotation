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
include { CACHE_SPLIT as CACHE_SPLIT_SPLICEAI; CACHE_SPLIT as CACHE_SPLIT_PANGOLIN } from '../modules/local/annotation_cache.nf'
include { CACHE_WRITE as CACHE_WRITE_SPLICEAI; CACHE_WRITE as CACHE_WRITE_PANGOLIN } from '../modules/local/annotation_cache.nf'
include { CHUNK_VCF as CHUNK_SPLICEAI_MISSES; CHUNK_VCF as CHUNK_PANGOLIN_MISSES } from '../modules/local/chunk.nf'
include { FAN_CHUNKS as FAN_SPLICEAI_MISSES; FAN_CHUNKS as FAN_PANGOLIN_MISSES } from './chunk_fan_concat.nf'

// Groups a predictor's chunks by run and pairs them with the pre-chunk input for the parity check.
def collectForConcat(annotated_chunk_ch, orig_vcf_ch, suffix) {
    annotated_chunk_ch
        // Hash order, not completion order: the same chunk list on every launch, so -resume keeps CONCAT_* cached.
        .groupTuple(sort: 'hash')
        .join(orig_vcf_ch)
        .map { meta, chunk_files, vcf, tbi -> tuple(meta, chunk_files, vcf, tbi, suffix) }
}

// The cache segment files listed for this run by run_annotate_vcf.sh (frozen per run, so a
// -resume does not see the run's own new segment), or [] when the cache has none yet.
def cacheSegments(list_file) {
    Channel.value(file(list_file).readLines().findAll { it.trim() }.collect { file(it) })
}

// CACHE_SPLIT misses that still need the predictor, as (meta, vcf, tbi).
def withMisses(misses_ch) {
    misses_ch
        .filter { meta, vcf, tbi, n -> (n as long) > 0 }
        .map { meta, vcf, tbi, n -> tuple(meta, vcf, tbi) }
}

workflow GENE_RESTRICTED_WHOLE_SUBWORKFLOW {
    take:
    meta_vcf_ch   // tuple(meta, vcf, tbi) — the FULL, unfiltered whole-run input

    main:
    GENE_SUBSET(meta_vcf_ch, file(params.gene_restriction_bed))
    subset_ch = GENE_SUBSET.out.vcf   // tuple(meta, subset.vcf.gz, subset.vcf.gz.tbi)

    CHUNK_VCF(subset_ch, params.gene_restricted_chunk_size)
    chunk_ch = FAN_CHUNKS(CHUNK_VCF.out.chunks)

    // With the score cache, SpliceAI and Pangolin only see the sites their cache lacks; the cached
    // sites (hits) rejoin their chunks at CONCAT, whose parity check covers both.
    if (params.annotation_cache_dir) {
        CACHE_SPLIT_SPLICEAI(subset_ch, 'spliceai', cacheSegments(params.annotation_cache_segments_spliceai))
        CHUNK_SPLICEAI_MISSES(withMisses(CACHE_SPLIT_SPLICEAI.out.misses), params.gene_restricted_chunk_size)
        spliceai_chunk_ch = FAN_SPLICEAI_MISSES(CHUNK_SPLICEAI_MISSES.out.chunks)
        spliceai_hits_ch = CACHE_SPLIT_SPLICEAI.out.hits

        CACHE_SPLIT_PANGOLIN(subset_ch, 'pangolin', cacheSegments(params.annotation_cache_segments_pangolin))
        CHUNK_PANGOLIN_MISSES(withMisses(CACHE_SPLIT_PANGOLIN.out.misses), params.gene_restricted_chunk_size)
        pangolin_chunk_ch = FAN_PANGOLIN_MISSES(CHUNK_PANGOLIN_MISSES.out.chunks)
        pangolin_hits_ch = CACHE_SPLIT_PANGOLIN.out.hits
    } else {
        spliceai_chunk_ch = chunk_ch
        pangolin_chunk_ch = chunk_ch
        spliceai_hits_ch = Channel.empty()
        pangolin_hits_ch = Channel.empty()
    }

    // SPLICEAI_ANNOTATE takes (meta, vcf) without the index.
    SPLICEAI_ANNOTATE(spliceai_chunk_ch.map { m, vcf, tbi -> tuple(m, vcf) })
    CONCAT_SPLICEAI(collectForConcat(SPLICEAI_ANNOTATE.out.annotated_chunk.mix(spliceai_hits_ch), subset_ch, 'annSpliceAI'))

    PANGOLIN_ANNOTATE(pangolin_chunk_ch)
    CONCAT_PANGOLIN(collectForConcat(PANGOLIN_ANNOTATE.out.vcf.map { m, vcf, tbi -> tuple(m, vcf) }.mix(pangolin_hits_ch), subset_ch, 'annPangolin'))

    // New segments only from sites the predictor computed in this run; stub outputs carry no scores.
    if (params.annotation_cache_dir && params.annotation_cache_write && !workflow.stubRun) {
        CACHE_WRITE_SPLICEAI(CONCAT_SPLICEAI.out.vcf.join(withMisses(CACHE_SPLIT_SPLICEAI.out.misses)),
            'spliceai', params.annotation_cache_version_spliceai)
        CACHE_WRITE_PANGOLIN(CONCAT_PANGOLIN.out.vcf.join(withMisses(CACHE_SPLIT_PANGOLIN.out.misses)),
            'pangolin', params.annotation_cache_version_pangolin)
    }

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
