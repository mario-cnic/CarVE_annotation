// Whole-VCF redesign (plan item 7): FAN_CHUNKS is the reusable half of the chunk-fan-out /
// collect-and-concat pattern Phase 1's spliceai_spike.nf first hand-wired — pure channel
// operators, no process call inside, so it's safe to call from multiple tier subworkflows without
// Nextflow's "process already used in this context" restriction (confirmed by hitting that error
// when CONCAT_CHUNKS was wrapped in a shared collect+concat subworkflow instead: Nextflow
// disallows the same process being reused across call sites in one workflow scope unless each
// site imports it under its own alias). The collect+concat half is NOT extracted into a shared
// subworkflow for that reason — each tier workflow (broad_pass_whole.nf,
// gene_restricted_whole.nf) inlines its own groupTuple/sort/join sequence and imports
// CONCAT_CHUNKS under its own per-tool alias.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

workflow FAN_CHUNKS {
    take:
    chunk_ch   // tuple(meta, vcf_list, tbi_list) — CHUNK_VCF's raw output (both globs sort into
               // matching per-index order, confirmed: chunk_NNNN.vcf.gz / chunk_NNNN.vcf.gz.tbi
               // share the same zero-padded index prefix)

    main:
    // Same List-normalization fix as Phase 1's spliceai_spike.nf: a single-chunk run's
    // path(glob) output is a scalar Path (which java.nio.file.Path iterates as path segments,
    // not a one-element list) unless explicitly wrapped first — applied to both globs here.
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
    out   // tuple(meta, vcf, tbi) — one channel element per chunk, matching the 3-tuple shape
          // VEP_ANNOTATE/BRANCHPOINT_ANNOTATE/PANGOLIN_ANNOTATE/SPIP_ANNOTATE/SPLICEAI_ANNOTATE
          // already expect
}
