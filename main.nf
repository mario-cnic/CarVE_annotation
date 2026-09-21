#!/usr/bin/env nextflow
// Phase 2: gene-restricted tier complete — SpliceAI + Pangolin + SPiP, legacy per-gene adapter.
// Does NOT replace main.sh. See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md.
//
// Usage:
//   nextflow run main.nf --input_vcf <gene>.vcf.gz [-profile standard|local_dev] [-stub-run]
//
// --input_vcf must already be tabix-indexed alongside it (<input_vcf>.tbi), matching today's
// Step-1 (variant_converter) output shape — this spike starts downstream of that step.

include { GENE_RESTRICTED_SUBWORKFLOW } from './workflows/gene_restricted.nf'

workflow {
    if (!params.input_vcf) {
        error "Missing required --input_vcf <gene>.vcf.gz (must have a sibling .tbi index)"
    }

    input_vcf_file = file(params.input_vcf)
    input_tbi_file = file("${params.input_vcf}.tbi")
    if (!input_vcf_file.exists())  { error "Input VCF not found: ${input_vcf_file}" }
    if (!input_tbi_file.exists())  { error "Missing tabix index: ${input_tbi_file}" }

    // Legacy adapter: gene identity derived from the input filename, same convention as
    // main.sh:363-364, since this spike only exercises the gene-restricted-panel special case.
    gene_name = input_vcf_file.baseName.replaceAll(/\.vcf$/, '').split('_')[0].toUpperCase()

    meta = [
        partition_type : 'gene',
        partition_id   : gene_name,
        gene           : gene_name,
        transcript     : null,
        build          : 'GRCh38',
        run_id         : workflow.sessionId.toString()
    ]

    meta_vcf_ch = Channel.of(tuple(meta, input_vcf_file, input_tbi_file))

    GENE_RESTRICTED_SUBWORKFLOW(meta_vcf_ch)

    GENE_RESTRICTED_SUBWORKFLOW.out.spliceai.view { m, vcf, tbi -> "SpliceAI done: ${m.partition_id} -> ${vcf}" }
    GENE_RESTRICTED_SUBWORKFLOW.out.pangolin.view { m, vcf, tbi -> "Pangolin done: ${m.partition_id} -> ${vcf}" }
    GENE_RESTRICTED_SUBWORKFLOW.out.spip.view    { m, vcf, tbi -> "SPiP done: ${m.partition_id} -> ${vcf}" }
}
