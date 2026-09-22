// Phase 2: Pangolin splice predictor, ported from src/hpc/annotate_pangolin_vars.sh.
// Unlike SpliceAI, main.sh runs Pangolin directly on the whole per-gene VCF.gz (main.sh:458-463)
// — no chunking stage exists for it, so this is a single process, not a chunk/fan-out subworkflow.
// OMP/MKL/OPENBLAS_NUM_THREADS exports (added 2026-09-21, Phase 3 review) port
// annotate_pangolin_vars.sh's own `export ..._NUM_THREADS=${THREADS:-4}` — a real gap in the
// original Phase 2 commit, caught late: without it, torch/numpy oversubscribe threads past the
// SGE-allocated slot count on the cluster (unreproducible in this local sandbox).
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.
//
// Whole-VCF redesign (plan item 7-8): also reused per-CHUNK in workflows/gene_restricted_whole.nf —
// publishDir skips chunk-level calls via `saveAs` (see modules/local/vep.nf's comment for why
// `enabled:` doesn't work for this and `saveAs` does); only CONCAT_CHUNKS's reassembled file
// (modules/local/chunk.nf) is a real deliverable there.

process PANGOLIN_ANNOTATE {
    tag "${meta.partition_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.partition_id}", mode: 'copy',
        saveAs: { filename -> meta.partition_type == 'chunk' ? null : filename }

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
