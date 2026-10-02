// Genotypes of the annotated VCF as a long table (one row per variant x sample), and an optional
// wide GT_<sample> view of the main table. Logic: src/python/extract_genotypes.py and
// src/python/genotypes_to_wide.py.
//
// EXTRACT_GENOTYPES writes nothing when the input has no sample columns or no GT field
// (params.genotypes = 'off' skips the process entirely).

process EXTRACT_GENOTYPES {
    tag "${meta.partition_id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.partition_id}", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.genotypes.pq"), emit: genotypes, optional: true

    script:
    """
    ${params.datasci_python} ${params.extract_genotypes_script} \\
        --vcf ${vcf} \\
        --output ${meta.partition_id}.genotypes.pq \\
        --bcftools ${params.bcftools}
    """

    stub:
    """
    touch ${meta.partition_id}.genotypes.pq
    """
}

process ADD_WIDE_GENOTYPES {
    tag "${meta.partition_id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.partition_id}", mode: 'copy'

    input:
    tuple val(meta), path(table, stageAs: "main_table_input.${params.output_format}"), path(genotypes)

    output:
    tuple val(meta), path("${meta.partition_id}.parsed.clean.withGT.${params.output_format}"), emit: table

    script:
    """
    ${params.datasci_python} ${params.genotypes_wide_script} \\
        --table ${table} \\
        --genotypes ${genotypes} \\
        --output ${meta.partition_id}.parsed.clean.withGT.${params.output_format} \\
        --max-samples ${params.genotype_wide_max_samples}
    """

    stub:
    """
    touch ${meta.partition_id}.parsed.clean.withGT.${params.output_format}
    """
}
