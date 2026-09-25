// Phase 4: merge all five predictors' output into one final annotated VCF per gene, ported from
// src/hpc/merge_vep_spip.sh. VEP is the authoritative record set (matches the plan's target
// architecture); SPiP, Pangolin, SpliceAI, and Branchpointer are layered onto it in that order via
// `bcftools annotate -a <file> -c <tag>`, each transferring exactly the INFO tag(s) that predictor
// writes (confirmed by reading known-good reference headers under test_data/test_run/annotation/:
// SPiP, Pangolin, SpliceAI each write one INFO tag of that literal name; Branchpointer writes the
// five explicit fields listed below). Still the trivial single-gene case (Phase 4 scope per the
// plan) — no cross-gene overlap can occur yet, so the AMBIGUOUS_GENE_SOURCE dedup logic the plan
// describes for Phase 7 isn't needed here.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.
//
// Two deliberate departures from the bash original, not silent ports:
// 1. `large_sv_vcf` handling is dropped entirely. Confirmed dead code (BUG_TRACKER.md SV-1): no
//    script anywhere in this repo ever produces a `*.large_svs.vcf.gz` file, so the bash branch
//    that concatenates it back in is permanently a no-op today. Reviving SV/CNV handling is
//    PLT-102 (future, separate), not something to silently recreate as an always-skipped branch
//    here.
// 2. No per-predictor soft-continue. The bash original wraps each `bcftools annotate` in
//    `if ... ; then ... else echo Warning ...; fi`, so the merge "succeeds" with a predictor
//    silently missing from the output. That, plus ORCH-3 (the bash script's `CMD_EXIT_CODE=$?` is
//    captured after its cleanup `rm` loop, not after the actual merge), together mean a real merge
//    failure today can go completely unnoticed by main.sh's dependency chain. This port drops the
//    soft-continue entirely — any `bcftools annotate` failure fails the process. Note this is a
//    channel/logic decision, not a shell-flag one: Nextflow already runs every script under
//    `#!/bin/bash -ue` by default (confirmed by inspecting a real .command.sh in this sandbox), so
//    a missing/malformed predictor file failing the whole task is Nextflow's native behavior, not
//    something this script has to opt into — ORCH-3 is closed because exit status is now
//    determined natively rather than from a manually captured `$?` after unrelated cleanup.

process MERGE_ANNOTATIONS {
    tag "${meta.partition_id}"
    label 'process_single'
    publishDir "${params.outdir}/${meta.partition_id}", mode: 'copy'

    input:
    tuple val(meta), path(vep_vcf), path(vep_tbi), path(bp_vcf), path(bp_tbi), path(pango_vcf), path(pango_tbi), path(sai_vcf), path(sai_tbi), path(spip_vcf), path(spip_tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annotated.vcf.gz"), path("${meta.partition_id}.annotated.vcf.gz.tbi"), emit: vcf

    // `pipefail` isn't part of Nextflow's default `-ue` wrapper (unlike `-e`/`-u`, which are
    // already applied and would make a redundant `set -e`/`-u` here misleading, not just inert),
    // so it's the one flag actually worth adding explicitly — without it, a failure inside a
    // pipe's left-hand side (none exist in this script today, but future edits could add one)
    // would go unnoticed.
    script:
    """
    set -o pipefail

    ${params.bcftools} annotate -a ${spip_vcf} -c SPiP ${vep_vcf} -Oz -o tmp_spip.vcf.gz
    ${params.tabix} -f tmp_spip.vcf.gz

    ${params.bcftools} annotate -a ${pango_vcf} -c Pangolin,PangolinTissue tmp_spip.vcf.gz -Oz -o tmp_pangolin.vcf.gz
    ${params.tabix} -f tmp_pangolin.vcf.gz

    ${params.bcftools} annotate -a ${sai_vcf} -c SpliceAI tmp_pangolin.vcf.gz -Oz -o tmp_spliceai.vcf.gz
    ${params.tabix} -f tmp_spliceai.vcf.gz

    ${params.bcftools} annotate -a ${bp_vcf} -c Branchpointer_prob,Branchpointer_U2_energy,Branchpoint_disrupted,LaBranchoR_score,LaBranchoR_acc_dist tmp_spliceai.vcf.gz -Oz -o ${meta.partition_id}.annotated.vcf.gz
    ${params.tabix} -f ${meta.partition_id}.annotated.vcf.gz
    """

    // Real run needs all five predictors' actual INFO tags to exist in their headers (stub outputs
    // upstream are plain `cp` of the pre-annotation VCF, with no such tags), so this process can't
    // be exercised end-to-end for real through Nextflow itself in this sandbox. Its exact command
    // sequence WAS verified for real, though (outside Nextflow, same commands this script runs):
    // ran against the five known-good MYBPC3 reference predictor outputs under
    // test_data/test_run/annotation/ and diffed the result against that same directory's
    // MYBPC3.annotated.vcf.gz — all 22 data records identical. This is the strongest verification
    // in the migration so far: it confirms the merge order, `-c` tag list, and bcftools invocation
    // are correct, not just that the DAG wires up. Stub (below) proves the 5-way channel join only.
    stub:
    """
    cp ${vep_vcf} ${meta.partition_id}.annotated.vcf.gz
    cp ${vep_tbi} ${meta.partition_id}.annotated.vcf.gz.tbi
    """
}
