// Branchpointer + LaBranchoR: lookup against a precomputed genome-wide branch-point BED (broad
// tier).

process BRANCHPOINT_ANNOTATE {
    tag "${meta.partition_id}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annBranchpoint.vcf.gz"), path("${meta.partition_id}.annBranchpoint.vcf.gz.tbi"), emit: vcf

    script:
    """
    zcat ${vcf} > raw_input.vcf
    ${params.branchpointer_python} ${params.branchpointer_script} \\
        raw_input.vcf ${meta.partition_id}.annBranchpoint.vcf \\
        --labranchor-bed ${params.labranchor_bed}
    ${params.bgzip} -f ${meta.partition_id}.annBranchpoint.vcf
    ${params.tabix} -f ${meta.partition_id}.annBranchpoint.vcf.gz
    """

    stub:
    """
    cp ${vcf} ${meta.partition_id}.annBranchpoint.vcf.gz
    cp ${tbi} ${meta.partition_id}.annBranchpoint.vcf.gz.tbi
    """
}
