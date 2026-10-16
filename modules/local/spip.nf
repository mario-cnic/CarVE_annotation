// SPiP splice predictor (gene-restricted tier by default, see params.spip_tier). BLAS/OpenMP
// threads are pinned to 1 because SPiP parallelises with --threads.

process SPIP_ANNOTATE {
    tag "${meta.partition_id}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annSPiP.vcf.gz"), path("${meta.partition_id}.annSPiP.vcf.gz.tbi"), emit: vcf

    script:
    def raw_out = "${meta.partition_id}.annSPiP.vcf"
    """
    zcat ${vcf} > raw_input.vcf
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    ${params.spip_rscript} ${params.spip_script} \\
        --input raw_input.vcf \\
        --output ${raw_out} \\
        -g hg38 \\
        -t ${task.cpus} \\
        --maxLines ${params.spip_max_lines} \\
        --VCF

    # Ported from annotate_spip_vars.sh's post-processing: SPiP emits a non-standard header
    # line that breaks downstream VCF parsers.
    grep -v "^##SPiP output v2.1\$" ${raw_out} > ${raw_out}.clean
    mv ${raw_out}.clean ${raw_out}

    ${params.bgzip} -f ${raw_out}
    ${params.tabix} -f ${raw_out}.gz
    """

    stub:
    """
    cp ${vcf} ${meta.partition_id}.annSPiP.vcf.gz
    cp ${tbi} ${meta.partition_id}.annSPiP.vcf.gz.tbi
    """
}
