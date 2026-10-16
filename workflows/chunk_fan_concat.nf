// Turns CHUNK_VCF's (meta, vcf_list, tbi_list) into one (meta, vcf, tbi) element per chunk.

workflow FAN_CHUNKS {
    take:
    chunk_ch   // tuple(meta, vcf_list, tbi_list)

    main:
    // A single chunk arrives as a bare Path, so wrap it in a list.
    out = chunk_ch
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

    emit:
    out   // tuple(meta, vcf, tbi)
}
