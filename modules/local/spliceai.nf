// Phase 1 spike: SpliceAI chunking + annotation + concat, ported from
// src/hpc/annotate_spliceai_vars.sh + src/python/split_vcf_chunks.py.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

process SPLICEAI_CHUNK {
    tag "${meta.partition_id}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("chunks/chunk_*.vcf.gz"), emit: chunks

    // No stub: chunking needs no FASTA/model and runs for real even under -stub-run,
    // so this stage's real behavior is genuinely verified in this sandbox.
    script:
    """
    mkdir -p chunks
    ${params.python_spliceai} ${params.split_vcf_chunks_script} \\
        --input ${vcf} \\
        --output-dir chunks \\
        --chunk-size ${params.chunk_size}
    """
}

process SPLICEAI_ANNOTATE {
    tag "${meta.partition_id}/${chunk.baseName}"
    // process_long, not process_medium: real WGS-scale wall-clock kill found 2026-09-25
    // (BUG_TRACKER.md MISC-11) -- see nextflow.config's process_long comment for the qacct
    // evidence (time problem, not memory).
    label 'process_long'

    input:
    tuple val(meta), path(chunk)

    output:
    tuple val(meta), path("${chunk.baseName.replace('.vcf','')}.annSpliceAI.vcf.gz"), emit: annotated_chunk

    script:
    def out_vcf = "${chunk.baseName.replace('.vcf','')}.annSpliceAI.vcf"
    """
    zcat ${chunk} > raw_input.vcf
    ${params.python_spliceai} ${params.annotate_spliceai_script} \\
        raw_input.vcf ${out_vcf} ${params.fasta} -d ${params.spliceai_distance}
    ${params.bgzip} -c ${out_vcf} > ${out_vcf}.gz
    """

    // Real model inference needs the GRCh38 FASTA + SpliceAI weights, unavailable in this
    // sandbox. Stub proves the DAG/channel wiring without them (run with -stub-run).
    stub:
    def out_vcf = "${chunk.baseName.replace('.vcf','')}.annSpliceAI.vcf"
    """
    cp ${chunk} ${out_vcf}.gz
    """
}

process SPLICEAI_CONCAT {
    tag "${meta.partition_id}"
    label 'process_single'
    publishDir "${params.outdir}/${meta.partition_id}", mode: 'copy'

    input:
    tuple val(meta), path(annotated_chunks), path(orig_vcf), path(orig_tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annSpliceAI.vcf.gz"), path("${meta.partition_id}.annSpliceAI.vcf.gz.tbi"), emit: vcf

    // No stub: concat + the record-count parity check (ported from
    // annotate_spliceai_vars.sh:170-171) run for real regardless of -stub-run, so chunk
    // reassembly is genuinely verified even when SPLICEAI_ANNOTATE itself is stubbed.
    script:
    def final_vcf = "${meta.partition_id}.annSpliceAI.vcf.gz"
    // No -a/--allow-overlaps: chunks are non-overlapping by construction (sequential
    // record-count split of one file, not a region-based split), and -a requires each input
    // to be independently tabix-indexed to check for overlaps — confirmed by hitting
    // "Could not retrieve index file" when -a was used without per-chunk indices.
    """
    ${params.bcftools} concat ${annotated_chunks.join(' ')} -Oz -o ${final_vcf}
    ${params.tabix} -p vcf ${final_vcf}

    orig_count=\$(${params.bcftools} view -H ${orig_vcf} | wc -l)
    final_count=\$(${params.bcftools} view -H ${final_vcf} | wc -l)
    if [ "\$orig_count" != "\$final_count" ]; then
        echo "SPLICEAI_CONCAT: record count mismatch — input=\$orig_count output=\$final_count" >&2
        exit 1
    fi
    echo "SPLICEAI_CONCAT: record count parity confirmed (\$final_count records)"
    """
}
