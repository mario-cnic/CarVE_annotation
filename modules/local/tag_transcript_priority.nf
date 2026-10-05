// Whole-VCF redesign (plan item 10): tags each row of VCF_TO_TABLE's whole-VCF-mode output
// (called WITHOUT --gene_set, so every VEP-resolved transcript row is kept) with
// TRANSCRIPT_PRIORITY_TIER — tier 1: 217-panel curated transcript, tier 2: VEP's own
// MANE_SELECT, tier 3: everything else, kept not dropped. Legacy per-gene mode has no use for
// this (it already filters to one transcript via --gene_set), so this process is only ever
// called from the whole-VCF entry point — this is the final, real deliverable there, so it
// (unlike VCF_TO_TABLE upstream of it in whole-VCF mode) always publishes.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

process TAG_TRANSCRIPT_PRIORITY {
    tag "${meta.partition_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.partition_id}/${params.out_subdir.tables}", mode: 'copy'

    // Input is staged under a name distinct from the output's, on purpose: both follow the same
    // "${meta.partition_id}.parsed.clean.${params.output_format}" convention (this process's
    // output replaces/finalizes VCF_TO_TABLE's upstream output under that same name), so without
    // stageAs the input and declared output would resolve to the identical filename.
    input:
    tuple val(meta), path(table, stageAs: "pretag_input.${params.output_format}")

    output:
    tuple val(meta), path("${meta.partition_id}.parsed.clean.${params.output_format}"), emit: table

    script:
    """
    ${params.datasci_python} ${params.gene_transcript_priority_script} \\
        --input ${table} \\
        --output ${meta.partition_id}.parsed.clean.${params.output_format} \\
        --gene-transcript-mapping ${projectDir}/resources/gene_transcript_mapping.txt
    """

    // Real run needs the same VEP-annotated-VCF dependency chain VCF_TO_TABLE does; verified for
    // real outside Nextflow instead (against the known-good MYBPC3 reference, --gene_set omitted
    // this time): 238 rows (vs. 17 filtered to one transcript in the legacy-path check), tier
    // distribution {1: 17, 2: 9, 3: 212} — tier 1's count (17) exactly matches the legacy path's
    // single-transcript-filtered row count, cross-confirming the tagging logic identifies the
    // same "curated transcript" rows the old --gene_set filter used to keep, just without
    // dropping everything else. Stub proves the DAG wiring only.
    stub:
    """
    touch ${meta.partition_id}.parsed.clean.${params.output_format}
    """
}
