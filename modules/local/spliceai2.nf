// SpliceAI2 (broad tier): exact-allele lookup in the precomputed GRCh38 scores (informational).

process SPLICEAI2_ANNOTATE {
    tag "${meta.partition_id}"
    label 'process_single'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annSpliceAI2.vcf.gz"), path("${meta.partition_id}.annSpliceAI2.vcf.gz.tbi"), emit: vcf

    script:
    """
    # code ${params.code_md5.spliceai2}
    # resources ${params.resource_id.spliceai2}
    ${params.spliceai2_python} ${params.spliceai2_script} \\
        ${vcf} ${meta.partition_id}.annSpliceAI2.vcf \\
        --scores-dir ${params.spliceai2_scores_dir}
    ${params.bgzip} -f ${meta.partition_id}.annSpliceAI2.vcf
    ${params.tabix} -f ${meta.partition_id}.annSpliceAI2.vcf.gz
    """

    stub:
    """
    cp ${vcf} ${meta.partition_id}.annSpliceAI2.vcf.gz
    cp ${tbi} ${meta.partition_id}.annSpliceAI2.vcf.gz.tbi
    """
}
