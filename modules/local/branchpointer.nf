// LaBranchoR (broad tier): positional lookup against the precomputed top-branchpoint BED (PRED-18) and,
// for SNVs, the allele-aware in silico mutagenesis scores.

process BRANCHPOINT_ANNOTATE {
    tag "${meta.partition_id}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annBranchpoint.vcf.gz"), path("${meta.partition_id}.annBranchpoint.vcf.gz.tbi"), emit: vcf

    script:
    """
    # code ${params.code_md5.branchpoint}
    # resources ${params.resource_id.branchpoint}
    zcat ${vcf} > raw_input.vcf
    ${params.branchpointer_python} ${params.branchpointer_script} \\
        raw_input.vcf ${meta.partition_id}.annBranchpoint.vcf \\
        --labranchor-bed ${params.labranchor_bed} \\
        --labranchor-ism ${params.labranchor_ism}
    ${params.bgzip} -f ${meta.partition_id}.annBranchpoint.vcf
    ${params.tabix} -f ${meta.partition_id}.annBranchpoint.vcf.gz
    """

    stub:
    """
    cp ${vcf} ${meta.partition_id}.annBranchpoint.vcf.gz
    cp ${tbi} ${meta.partition_id}.annBranchpoint.vcf.gz.tbi
    """
}
