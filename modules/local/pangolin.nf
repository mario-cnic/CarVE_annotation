// Pangolin splice predictor. Runs per-chunk in the gene-restricted tier (chunk boundaries don't
// affect correctness — each variant is scored independently against the FASTA).

process PANGOLIN_ANNOTATE {
    // Includes the chunk filename so Nextflow's progress display distinguishes concurrent chunks.
    tag "${meta.partition_id}/${vcf.baseName}"
    // Wall-clock, not memory-bound: per-chunk runtime is highly variable (minutes to many hours),
    // so this needs a much larger time ceiling than the other predictors' process_medium.
    label 'process_long'
    // Chunk-level runs are unpublished intermediates; only CONCAT_CHUNKS's reassembled output is.
    publishDir "${params.outdir}/${meta.partition_id}/${params.out_subdir.predictors}", mode: 'copy',
        saveAs: { filename -> meta.partition_type == 'chunk' ? null : filename }

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annPangolin.vcf.gz"), path("${meta.partition_id}.annPangolin.vcf.gz.tbi"), emit: vcf

    // Needs the GRCh38 FASTA + the model DB; stub proves DAG wiring only.
    script:
    def raw_out = "${meta.partition_id}.annPangolin.vcf"
    """
    zcat ${vcf} > raw_input.vcf
    # Without these, torch/numpy oversubscribe threads past the SGE-allocated slot count.
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    PYTHONPATH="${params.pangolin_repo}:\${PYTHONPATH:-}" \\
        ${params.pangolin_python} -m pangolin.pangolin \\
        raw_input.vcf ${params.fasta} ${params.pangolin_db} ${raw_out} \\
        -d ${params.pangolin_distance}
    ${params.bgzip} -f -@ ${task.cpus} ${raw_out}
    ${params.tabix} -f ${raw_out}.gz
    """

    stub:
    """
    cp ${vcf} ${meta.partition_id}.annPangolin.vcf.gz
    cp ${tbi} ${meta.partition_id}.annPangolin.vcf.gz.tbi
    """
}
