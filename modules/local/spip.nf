// Phase 2: SPiP splice predictor, ported from src/hpc/annotate_spip_vars.sh.
// Like Pangolin, main.sh runs SPiP directly on the whole per-gene VCF.gz (main.sh:420-429) — no
// chunking stage exists for it either; SPiPv2.1_main.r's own --maxLines is an internal detail of
// that script, kept verbatim as a non-goal (same treatment as SpliceAI's -d value).
// Tier: gene-restricted by default (params.spip_tier == 'restricted') per the user's explicit
// choice during plan approval — SPiP's broad-safety is unresolved, no profiling data exists yet.
// OMP/MKL/OPENBLAS_NUM_THREADS=1 exports (added 2026-09-21, Phase 3 review) port
// annotate_spip_vars.sh's own explicit `=1` pins (unconditional, unlike Pangolin's ${THREADS:-4}
// — SPiP does its own internal parallelism via -t/--threads, so these are deliberately pinned to
// 1 to prevent double-parallelizing) — a real gap in the original Phase 2 commit, caught late.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.
//
// Whole-VCF redesign (plan item 7-8): also reused per-CHUNK in both whole-VCF tier subworkflows
// (whichever one params.spip_tier selects) — publishDir skips chunk-level calls via `saveAs` (see
// modules/local/vep.nf's comment for why `enabled:` doesn't work for this and `saveAs` does); only
// CONCAT_CHUNKS's reassembled file (modules/local/chunk.nf) is a real deliverable there.

process SPIP_ANNOTATE {
    tag "${meta.partition_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.partition_id}/${params.out_subdir.predictors}", mode: 'copy',
        saveAs: { filename -> meta.partition_type == 'chunk' ? null : filename }

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annSPiP.vcf.gz"), path("${meta.partition_id}.annSPiP.vcf.gz.tbi"), emit: vcf

    // Real inference needs the spip_env R stack + the external SPiPv2.1_main.r script, neither
    // executable in this sandbox (mirrors SPLICEAI_ANNOTATE's stub rationale). Stub proves the
    // DAG wiring only; real output must be checked on the cluster.
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
