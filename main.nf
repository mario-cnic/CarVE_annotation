#!/usr/bin/env nextflow
// Phase 5: both tiers + merge + VCF_TO_TABLE — SpliceAI + Pangolin + SPiP (gene-restricted),
// VEP + Branchpointer (broad), merged into one final annotated VCF per gene, then flattened into
// the final parsed table. Legacy per-gene adapter. Does NOT replace main.sh. See
// /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md.
//
// Usage:
//   nextflow run main.nf --input_vcf <gene>.vcf.gz [-profile standard|local_dev] [-stub-run]
//
// --input_vcf must already be tabix-indexed alongside it (<input_vcf>.tbi), matching today's
// Step-1 (variant_converter) output shape — this spike starts downstream of that step.

include { GENE_RESTRICTED_SUBWORKFLOW } from './workflows/gene_restricted.nf'
include { BROAD_PASS_SUBWORKFLOW } from './workflows/broad_pass.nf'
include { MERGE_SUBWORKFLOW } from './workflows/merge_subworkflow.nf'
include { VCF_TO_TABLE } from './modules/local/vcf_to_table.nf'

workflow {
    if (!params.input_vcf) {
        error "Missing required --input_vcf <gene>.vcf.gz (must have a sibling .tbi index)"
    }

    // gene_restricted.nf and broad_pass.nf each gate SPIP_ANNOTATE behind an if/else that falls
    // back to Channel.empty() for any params.spip_tier value other than their own tier name — so
    // an invalid value silently emits from NEITHER, MERGE_SUBWORKFLOW's spip_ch never fires, and
    // the pipeline exits 0 with no final VCF for the gene, no error. Caught here instead.
    if (!(params.spip_tier in ['restricted', 'broad'])) {
        error "params.spip_tier must be 'restricted' or 'broad', got: ${params.spip_tier}"
    }

    input_vcf_file = file(params.input_vcf)
    input_tbi_file = file("${params.input_vcf}.tbi")
    if (!input_vcf_file.exists())  { error "Input VCF not found: ${input_vcf_file}" }
    if (!input_tbi_file.exists())  { error "Missing tabix index: ${input_tbi_file}" }

    // Legacy adapter: gene identity derived from the input filename, same convention as
    // main.sh:363-364, since this spike only exercises the gene-restricted-panel special case.
    gene_name = input_vcf_file.baseName.replaceAll(/\.vcf$/, '').split('_')[0].toUpperCase()

    // Transcript resolution for VCF_TO_TABLE, ported from main.sh:378's exact lookup: whole-word
    // match against resources/gene_transcript_mapping.txt (Gen,NM,ENST columns), 3rd field, or
    // 'UNKNOWN' if the gene has no mapping row. Deliberately NOT the live
    // query_new_transcripts.py --auto-append write path (ORCH-9 — unlocked write, non-MANE
    // fallback) — read-only lookup against the existing file, same as main.sh's read side.
    transcript = 'UNKNOWN'
    mapping_file = file("${projectDir}/resources/gene_transcript_mapping.txt")
    mapping_file.eachLine { line ->
        if (line =~ /\b${java.util.regex.Pattern.quote(gene_name)}\b/) {
            def fields = line.split(',')
            if (fields.size() >= 3) { transcript = fields[2].trim() }
        }
    }

    // Behavior change from main.sh, not a verbatim port: 'UNKNOWN' reaching VCF_TO_TABLE's
    // --gene_set isn't harmless there — vcf_parser_pysam.py's read_gene_set() treats it as a
    // literal one-element gene list, no VEP CSQ block ever has Feature == 'UNKNOWN', every variant
    // gets skipped, and filter_variants.py writes an empty-but-valid table — exit 0, published,
    // no warning. In main.sh that was one hand-driven gene's blast radius; here it's this DAG's
    // terminal stage, and at Phase 7 scale a missing mapping row would silently ship an empty
    // deliverable to Module 3. Failing loudly here, alongside the existing spip_tier guard.
    if (transcript == 'UNKNOWN') {
        error "No transcript mapping for ${gene_name} in resources/gene_transcript_mapping.txt — VCF_TO_TABLE's --gene_set would match nothing and silently produce an empty table."
    }

    meta = [
        partition_type : 'gene',
        partition_id   : gene_name,
        gene           : gene_name,
        transcript     : transcript,
        build          : 'GRCh38',
        run_id         : workflow.sessionId.toString()
    ]

    meta_vcf_ch = Channel.of(tuple(meta, input_vcf_file, input_tbi_file))

    GENE_RESTRICTED_SUBWORKFLOW(meta_vcf_ch)
    BROAD_PASS_SUBWORKFLOW(meta_vcf_ch)

    GENE_RESTRICTED_SUBWORKFLOW.out.spliceai.view { m, vcf, tbi -> "SpliceAI done: ${m.partition_id} -> ${vcf}" }
    GENE_RESTRICTED_SUBWORKFLOW.out.pangolin.view { m, vcf, tbi -> "Pangolin done: ${m.partition_id} -> ${vcf}" }
    GENE_RESTRICTED_SUBWORKFLOW.out.spip.view    { m, vcf, tbi -> "SPiP (restricted) done: ${m.partition_id} -> ${vcf}" }
    BROAD_PASS_SUBWORKFLOW.out.vep.view          { m, vcf, tbi -> "VEP done: ${m.partition_id} -> ${vcf}" }
    BROAD_PASS_SUBWORKFLOW.out.branchpoint.view  { m, vcf, tbi -> "Branchpoint done: ${m.partition_id} -> ${vcf}" }
    BROAD_PASS_SUBWORKFLOW.out.spip.view         { m, vcf, tbi -> "SPiP (broad) done: ${m.partition_id} -> ${vcf}" }

    // spip_tier picks exactly one of these two channels to actually emit (the other is
    // Channel.empty(), see gene_restricted.nf/broad_pass.nf) — .mix() combines them into the one
    // real stream MERGE_SUBWORKFLOW needs, without caring which tier was active.
    spip_merged_ch = GENE_RESTRICTED_SUBWORKFLOW.out.spip.mix(BROAD_PASS_SUBWORKFLOW.out.spip)

    MERGE_SUBWORKFLOW(
        BROAD_PASS_SUBWORKFLOW.out.vep,
        BROAD_PASS_SUBWORKFLOW.out.branchpoint,
        GENE_RESTRICTED_SUBWORKFLOW.out.pangolin,
        GENE_RESTRICTED_SUBWORKFLOW.out.spliceai,
        spip_merged_ch
    )

    MERGE_SUBWORKFLOW.out.vcf.view { m, vcf, tbi -> "MERGE done: ${m.partition_id} -> ${vcf}" }

    VCF_TO_TABLE(MERGE_SUBWORKFLOW.out.vcf)

    VCF_TO_TABLE.out.table.view { m, table -> "VCF_TO_TABLE done: ${m.partition_id} -> ${table}" }
}
