// Record-count chunking and reassembly shared by both predictor tiers.

process CHUNK_VCF {
    tag "${meta.partition_id}"
    label 'process_low'

    // split_vcf_chunks.py writes a .tbi next to each chunk; both globs sort into matching order.
    input:
    tuple val(meta), path(vcf), path(tbi)
    val(chunk_size)

    output:
    tuple val(meta), path("chunks/chunk_*.vcf.gz"), path("chunks/chunk_*.vcf.gz.tbi"), emit: chunks

    // No stub: chunking also runs under -stub-run. chunk_size is set per tier by the caller.
    script:
    """
    mkdir -p chunks
    ${params.python_spliceai} ${params.split_vcf_chunks_script} \\
        --input ${vcf} \\
        --output-dir chunks \\
        --chunk-size ${chunk_size}
    """
}

// Reassembles one predictor's chunks into a whole-run VCF and checks record-count parity
// against the pre-chunk input.
process CONCAT_CHUNKS {
    tag "${meta.partition_id}/${suffix}"
    label 'process_single'
    publishDir "${params.outdir}/${meta.partition_id}/${params.out_subdir.predictors}", mode: 'copy'

    // Chunks share the output's file name, so each is staged in its own numbered directory.
    input:
    tuple val(meta), path(annotated_chunks, stageAs: 'chunk_??/*'), path(orig_vcf), path(orig_tbi), val(suffix)

    output:
    tuple val(meta), path("${meta.partition_id}.${suffix}.vcf.gz"), path("${meta.partition_id}.${suffix}.vcf.gz.tbi"), emit: vcf

    // No stub: concat and the parity check also run under -stub-run.
    // Chunks may arrive in any order: index each, `concat -a`, then `sort` the result.
    script:
    def final_vcf = "${meta.partition_id}.${suffix}.vcf.gz"
    // A single staged file is a bare Path, which `.join` would split into path segments.
    def chunk_list = annotated_chunks instanceof List ? annotated_chunks : [annotated_chunks]
    """
    set -o pipefail
    for f in ${chunk_list.join(' ')}; do
        ${params.tabix} -f -p vcf "\$f"
    done
    # Concat streams uncompressed BCF into sort; sort temp files go in the task directory, not the node's TMPDIR.
    ${params.bcftools} concat -a ${chunk_list.join(' ')} -Ou \\
        | ${params.bcftools} sort -T ./bcftools-sort.XXXXXX -Oz -o ${final_vcf} -
    ${params.tabix} -p vcf ${final_vcf}

    # Record counts from the indexes; full decompression only if an index has no counts.
    count_records() { ${params.bcftools} index -n "\$1" 2>/dev/null || ${params.bcftools} view -H "\$1" | wc -l; }
    orig_count=\$(count_records ${orig_vcf})
    final_count=\$(count_records ${final_vcf})
    if [ "\$orig_count" != "\$final_count" ]; then
        echo "CONCAT_CHUNKS (${suffix}): record count mismatch — input=\$orig_count output=\$final_count" >&2
        exit 1
    fi
    echo "CONCAT_CHUNKS (${suffix}): record count parity confirmed (\$final_count records)"
    """
}
