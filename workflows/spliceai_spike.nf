// Phase 1 spike workflow: wires SPLICEAI_CHUNK -> SPLICEAI_ANNOTATE -> SPLICEAI_CONCAT with
// native Nextflow fan-out/collection, replacing the hand-rolled bash worker-pool in
// src/hpc/annotate_spliceai_vars.sh:120-152. See the migration plan for scope.

include { SPLICEAI_CHUNK; SPLICEAI_ANNOTATE; SPLICEAI_CONCAT } from '../modules/local/spliceai.nf'

workflow SPLICEAI_SPIKE {
    take:
    meta_vcf_ch   // tuple(meta, vcf, tbi) — one element per gene for this spike

    main:
    SPLICEAI_CHUNK(meta_vcf_ch)

    // Fan out: one channel element per chunk, native replacement for the bash worker-pool.
    // NOTE: a `path(glob)` process output emits a scalar Path (not a List) when the glob
    // matches exactly one file — and java.nio.file.Path implements Iterable<Path>, iterating
    // its path *segments*. Without this normalization, a single-chunk gene silently fans out
    // into one bogus "chunk" per directory component of the real file's path instead of one
    // real chunk. Confirmed by hitting this for real against MYBPC3's small test VCF.
    chunks_ch = SPLICEAI_CHUNK.out.chunks
        .map { meta, chunk_list -> tuple(meta, chunk_list instanceof List ? chunk_list : [chunk_list]) }
        .flatMap { meta, chunk_list -> chunk_list.collect { c -> tuple(meta, c) } }

    SPLICEAI_ANNOTATE(chunks_ch)

    // Collect annotated chunks back per gene, sorted by filename to preserve original record
    // order (matches annotate_spliceai_vars.sh:156-159's `sorted by filename` step) —
    // groupTuple() does not guarantee completion order, so this sort is not optional.
    grouped_ch = SPLICEAI_ANNOTATE.out.annotated_chunk
        .groupTuple()
        .map { meta, chunk_files -> tuple(meta, chunk_files.sort { it.name }) }

    // Rejoin the original pre-chunk VCF for SPLICEAI_CONCAT's record-count parity check.
    concat_input_ch = grouped_ch
        .join(meta_vcf_ch)   // (meta, chunk_files) join (meta, vcf, tbi) on meta -> (meta, chunk_files, vcf, tbi)

    SPLICEAI_CONCAT(concat_input_ch)

    emit:
    vcf = SPLICEAI_CONCAT.out.vcf
}
