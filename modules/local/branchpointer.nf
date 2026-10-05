// Phase 3: Branchpointer + LaBranchoR predictor pair, ported from
// src/hpc/annotate_branchpointer_vars.sh. Broad-safe by design: it's a lookup against a
// precomputed genome-wide BED catalog (labranchor_grch38_top.bed.gz), not live model inference —
// confirmed during this migration's research (see plan Context). Not containerized in this phase,
// same treatment as Pangolin/SPiP — the bash source reuses the plain spliceai_env python, not a
// container, and porting that choice is a non-goal here.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.
//
// Whole-VCF redesign (plan item 7-8): also reused per-CHUNK in workflows/broad_pass_whole.nf —
// publishDir skips chunk-level calls via `saveAs` (see modules/local/vep.nf's comment for why
// `enabled:` doesn't work for this and `saveAs` does); only CONCAT_CHUNKS's reassembled file
// (modules/local/chunk.nf) is a real deliverable there.

process BRANCHPOINT_ANNOTATE {
    tag "${meta.partition_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.partition_id}/${params.out_subdir.predictors}", mode: 'copy',
        saveAs: { filename -> meta.partition_type == 'chunk' ? null : filename }

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annBranchpoint.vcf.gz"), path("${meta.partition_id}.annBranchpoint.vcf.gz.tbi"), emit: vcf

    // Unlike every other predictor ported so far, this one's real command was actually run and
    // verified in this sandbox (outside Nextflow, same script/args this process invokes): against
    // test_data/test_run/_tmp/MYBPC3.vcf.gz, output was byte-identical (22/22 records) to the
    // known-good test_data/test_run/annotation/MYBPC3.annBranchpoint.vcf.gz reference — its only
    // real dependency is `pysam` against the LaBranchoR BED, both present locally. Still exercised
    // via -stub-run when run through the full main.nf pipeline (see gene_restricted/broad_pass
    // wiring), since Nextflow's -stub-run is an all-or-nothing switch for the whole run.
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
