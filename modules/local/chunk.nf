// Whole-VCF redesign (plan item 7): generic record-count chunk/concat, reused by BOTH tiers
// instead of each tool having its own chunking layer. Wraps the same, unmodified
// src/python/split_vcf_chunks.py that Phase 1's SPLICEAI_CHUNK already wraps — pure record-count
// chunking with no gene-awareness, so one implementation suffices for the broad tier's full-input
// chunks and the gene-restricted tier's subsetted-input chunks alike.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

process CHUNK_VCF {
    tag "${meta.partition_id}"
    label 'process_low'

    // split_vcf_chunks.py already tabix-indexes every chunk it writes (confirmed: real .tbi files
    // land alongside each chunk_NNNN.vcf.gz) — captured here (both globs sort into matching
    // per-index order) so downstream per-chunk calls to VEP_ANNOTATE/BRANCHPOINT_ANNOTATE/
    // PANGOLIN_ANNOTATE/SPIP_ANNOTATE (modules/local/{vep,branchpointer,pangolin,spip}.nf) get the
    // 3-tuple `(meta, vcf, tbi)` shape those processes already expect from the legacy per-gene
    // path, with no changes needed to those processes themselves.
    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("chunks/chunk_*.vcf.gz"), path("chunks/chunk_*.vcf.gz.tbi"), emit: chunks

    // No stub: chunking needs no FASTA/model and runs for real even under -stub-run (same
    // rationale as SPLICEAI_CHUNK, which this process supersedes for the whole-VCF path).
    script:
    """
    mkdir -p chunks
    ${params.python_spliceai} ${params.split_vcf_chunks_script} \\
        --input ${vcf} \\
        --output-dir chunks \\
        --chunk-size ${params.chunk_size}
    """
}

// Reassembles one tool's per-chunk outputs back into a single whole-run VCF, generalizing
// SPLICEAI_CONCAT (same bcftools concat + index + record-count parity check against the
// pre-chunk input) with an explicit `suffix` so one process serves all five predictors instead of
// one per tool. Always publishes — its output is always a real whole-run deliverable, unlike the
// per-chunk ANNOTATE calls upstream (see modules/local/{vep,branchpointer,pangolin,spip}.nf's
// conditional publishDir, which suppresses publishing those chunk-level intermediates).
process CONCAT_CHUNKS {
    tag "${meta.partition_id}/${suffix}"
    label 'process_single'
    publishDir "${params.outdir}/${meta.partition_id}", mode: 'copy'

    // Real bug caught in review: every per-chunk ANNOTATE output is named
    // "${meta.partition_id}.${suffix}.vcf.gz" (same convention the legacy per-gene path's
    // processes already use, see e.g. modules/local/vep.nf) — identical to what THIS process
    // declares as its own output. With a single-chunk run, the staged input and the declared
    // output resolved to the exact same filename, so `bcftools concat -o <name>` truncated its
    // own input to 0 bytes before reading it (`bcf_hdr_read` then fails on the now-empty file).
    // Staging each chunk into its own numbered subdirectory avoids both this collision and any
    // collision between chunks sharing a name — same mechanism already used for this exact class
    // of bug in modules/local/tag_transcript_priority.nf's `stageAs`.
    input:
    tuple val(meta), path(annotated_chunks, stageAs: 'chunk_??/*'), path(orig_vcf), path(orig_tbi), val(suffix)

    output:
    tuple val(meta), path("${meta.partition_id}.${suffix}.vcf.gz"), path("${meta.partition_id}.${suffix}.vcf.gz.tbi"), emit: vcf

    // No stub: concat + parity check run for real regardless of -stub-run, same rationale as
    // SPLICEAI_CONCAT (chunk reassembly is genuinely verified even when the per-chunk ANNOTATE
    // step itself is stubbed).
    script:
    def final_vcf = "${meta.partition_id}.${suffix}.vcf.gz"
    // Same scalar-vs-List Path gotcha Phase 1 first caught (java.nio.file.Path implements
    // Iterable<Path>, so `.join()` on a bare single Path iterates its path *segments* instead of
    // treating it as one element) — resurfaces here because the `chunk_??/*` stageAs directory
    // pattern above means a single-chunk run's `annotated_chunks` binds to one Path two segments
    // deep (`chunk_01/<file>`), and `.join(' ')` on that bare Path split it into two bogus
    // arguments ("chunk_01" and the real filename) rather than one.
    def chunk_list = annotated_chunks instanceof List ? annotated_chunks : [annotated_chunks]
    // No -a/--allow-overlaps: chunks are non-overlapping by construction (sequential
    // record-count split), same rationale as SPLICEAI_CONCAT.
    """
    ${params.bcftools} concat ${chunk_list.join(' ')} -Oz -o ${final_vcf}
    ${params.tabix} -p vcf ${final_vcf}

    orig_count=\$(${params.bcftools} view -H ${orig_vcf} | wc -l)
    final_count=\$(${params.bcftools} view -H ${final_vcf} | wc -l)
    if [ "\$orig_count" != "\$final_count" ]; then
        echo "CONCAT_CHUNKS (${suffix}): record count mismatch — input=\$orig_count output=\$final_count" >&2
        exit 1
    fi
    echo "CONCAT_CHUNKS (${suffix}): record count parity confirmed (\$final_count records)"
    """
}
