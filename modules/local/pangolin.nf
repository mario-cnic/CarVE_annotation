// Phase 2: Pangolin splice predictor, ported from src/hpc/annotate_pangolin_vars.sh.
// Unlike SpliceAI, main.sh runs Pangolin directly on the whole per-gene VCF.gz (main.sh:458-463)
// — no chunking stage exists for it, so this is a single process, not a chunk/fan-out subworkflow.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

process PANGOLIN_ANNOTATE {
    tag "${meta.partition_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.partition_id}", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annPangolin.vcf.gz"), path("${meta.partition_id}.annPangolin.vcf.gz.tbi"), emit: vcf

    // Real inference needs the GRCh38 FASTA + the ~1GB pangolin_grch38.db, unavailable in this
    // sandbox (mirrors SPLICEAI_ANNOTATE's stub rationale). Stub proves the DAG wiring only.
    script:
    def raw_out = "${meta.partition_id}.annPangolin.vcf"
    """
    zcat ${vcf} > raw_input.vcf
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
