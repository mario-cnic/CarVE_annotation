#!/usr/bin/env nextflow
// Whole-VCF redesign (plan items 7-10): the new PRIMARY entry point — one arbitrary-size VCF in,
// spanning any number of genes, annotated as fully as possible, one final table out. Replaces
// the legacy per-gene adapter (main.nf, kept working untouched as a deprioritized opt-in — see
// the plan's "Per-gene mode" section) as the path getting further development attention.
// Does NOT replace main.sh. See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md.
//
// Usage:
//   nextflow run annotate_vcf.nf --input_vcf <any>.vcf.gz [-profile standard|local_dev] [-stub-run]
//
// --input_vcf must already be tabix-indexed alongside it (<input_vcf>.tbi).

include { CHECK_INPUT_ASSEMBLY } from './modules/local/check_assembly.nf'
include { BROAD_PASS_WHOLE_SUBWORKFLOW } from './workflows/broad_pass_whole.nf'
include { GENE_RESTRICTED_WHOLE_SUBWORKFLOW } from './workflows/gene_restricted_whole.nf'
include { MERGE_SUBWORKFLOW } from './workflows/merge_subworkflow.nf'
include { VCF_TO_TABLE } from './modules/local/vcf_to_table.nf'
include { TAG_TRANSCRIPT_PRIORITY } from './modules/local/tag_transcript_priority.nf'
include { EXTRACT_GENOTYPES; ADD_WIDE_GENOTYPES } from './modules/local/genotypes.nf'

workflow {
    if (!params.input_vcf) {
        error "Missing required --input_vcf <file>.vcf.gz (must have a sibling .tbi index)"
    }
    if (!(params.spip_tier in ['restricted', 'broad'])) {
        error "params.spip_tier must be 'restricted' or 'broad', got: ${params.spip_tier}"
    }

    input_vcf_file = file(params.input_vcf)
    input_tbi_file = file("${params.input_vcf}.tbi")
    if (!input_vcf_file.exists()) { error "Input VCF not found: ${input_vcf_file}" }
    if (!input_tbi_file.exists()) { error "Missing tabix index: ${input_tbi_file}" }

    // `partition_type: 'chunk'` marks this whole run as the whole-VCF path everywhere it matters
    // downstream: it's what VEP_ANNOTATE/BRANCHPOINT_ANNOTATE/PANGOLIN_ANNOTATE/SPIP_ANNOTATE
    // check to suppress publishing their per-chunk intermediates (modules/local/{vep,
    // branchpointer,pangolin,spip}.nf), and what VCF_TO_TABLE checks to suppress publishing its
    // pre-tagging intermediate (modules/local/vcf_to_table.nf). No `gene`/`transcript` fields —
    // transcript selection happens per-row in TAG_TRANSCRIPT_PRIORITY, not per-partition (see the
    // plan's "Transcript priority tiering" section).
    //
    // run_id (closes BUG_TRACKER.md MISC-7): derived from the input filename by default, same
    // .baseName.replaceAll(/\.vcf$/, '') pattern main.nf already uses for gene_name — NOT a
    // Nextflow session UUID (workflow.sessionId), which is fresh and meaningless on every
    // invocation regardless of input, so two runs of the same sample produced two unrelated,
    // unfindable output directories. Deterministic naming also means a re-run of the same input
    // lands in the same output directory (matching main.sh's own is_valid_file-checkpoint
    // philosophy) rather than scattering a new one every time. --run_id overrides this when the
    // filename itself isn't a meaningful identifier (e.g. a generic "input.vcf.gz").
    run_id = params.run_id ?: input_vcf_file.baseName.replaceAll(/\.vcf$/, '')
    meta = [
        partition_type : 'chunk',
        partition_id   : run_id,
        build           : 'GRCh38',
        run_id          : run_id
    ]

    // Every predictor consumes CHECK_INPUT_ASSEMBLY's output, so nothing is scheduled until the
    // input's contig lengths and REF alleles are confirmed against params.fasta (MISC-14).
    CHECK_INPUT_ASSEMBLY(Channel.of(tuple(meta, input_vcf_file, input_tbi_file)))
    meta_vcf_ch = CHECK_INPUT_ASSEMBLY.out.vcf

    BROAD_PASS_WHOLE_SUBWORKFLOW(meta_vcf_ch)
    GENE_RESTRICTED_WHOLE_SUBWORKFLOW(meta_vcf_ch)

    // spip_tier picks exactly one of these two channels to actually emit — same .mix() pattern
    // as the legacy path's main.nf.
    spip_merged_ch = GENE_RESTRICTED_WHOLE_SUBWORKFLOW.out.spip.mix(BROAD_PASS_WHOLE_SUBWORKFLOW.out.spip)

    // MERGE_SUBWORKFLOW is reused verbatim from the legacy path (workflows/merge_subworkflow.nf)
    // — its bcftools-annotate chain doesn't care whether its five inputs are one-gene files or
    // whole-run files, only that they're each one file with the right INFO tags (see the plan's
    // "Why the merge logic itself doesn't need to change" section).
    MERGE_SUBWORKFLOW(
        BROAD_PASS_WHOLE_SUBWORKFLOW.out.vep,
        BROAD_PASS_WHOLE_SUBWORKFLOW.out.branchpoint,
        GENE_RESTRICTED_WHOLE_SUBWORKFLOW.out.pangolin,
        GENE_RESTRICTED_WHOLE_SUBWORKFLOW.out.spliceai,
        spip_merged_ch
    )

    MERGE_SUBWORKFLOW.out.vcf.view { m, vcf, tbi -> "MERGE done: ${m.partition_id} -> ${vcf}" }

    // VCF_TO_TABLE runs with --gene_set omitted (no meta.transcript in this path's meta) — keeps
    // every VEP-resolved transcript row instead of filtering to one. TAG_TRANSCRIPT_PRIORITY then
    // adds TRANSCRIPT_PRIORITY_TIER and is the real final deliverable (VCF_TO_TABLE's own publish
    // is suppressed here, see modules/local/vcf_to_table.nf).
    VCF_TO_TABLE(MERGE_SUBWORKFLOW.out.vcf)
    TAG_TRANSCRIPT_PRIORITY(VCF_TO_TABLE.out.table)

    TAG_TRANSCRIPT_PRIORITY.out.table.view { m, table -> "Final table: ${m.partition_id} -> ${table}" }

    // Long genotype table from the merged VCF; no output when the input has no genotypes.
    if (params.genotypes != 'off') {
        EXTRACT_GENOTYPES(MERGE_SUBWORKFLOW.out.vcf)
        if (params.genotype_wide) {
            ADD_WIDE_GENOTYPES(TAG_TRANSCRIPT_PRIORITY.out.table.join(EXTRACT_GENOTYPES.out.genotypes))
            EXTRACT_GENOTYPES.out.genotypes.count().subscribe { n ->
                if (n == 0) log.warn "genotype_wide: the input has no genotypes (no samples or no GT), so no withGT table is written"
            }
        }
    }
}
