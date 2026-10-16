// SpliceAI annotation, one task per chunk (gene-restricted tier).

process SPLICEAI_ANNOTATE {
    tag "${meta.partition_id}/${chunk.baseName}"
    // Wall-clock, not memory-bound: per-chunk runtime is highly variable (minutes to many hours).
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
