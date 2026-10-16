// Adds TRANSCRIPT_PRIORITY_TIER to every row (1 = curated 217-panel transcript, 2 = MANE Select,
// 3 = other; no row dropped) and publishes the final table.

process TAG_TRANSCRIPT_PRIORITY {
    tag "${meta.partition_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.partition_id}/${params.out_subdir.tables}", mode: 'copy'

    // Input and output share a file name, so the input is staged under another.
    input:
    tuple val(meta), path(table, stageAs: "pretag_input.${params.output_format}")

    output:
    tuple val(meta), path("${meta.partition_id}.parsed.clean.${params.output_format}"), emit: table

    script:
    def stamp_flag = params.provenance_stamp ? "--stamp ${params.provenance_stamp}" : ""
    """
    ${params.datasci_python} ${params.gene_transcript_priority_script} \\
        --input ${table} \\
        --output ${meta.partition_id}.parsed.clean.${params.output_format} \\
        --gene-transcript-mapping ${projectDir}/resources/gene_transcript_mapping.txt \\
        ${stamp_flag}
    """

    stub:
    """
    touch ${meta.partition_id}.parsed.clean.${params.output_format}
    """
}
