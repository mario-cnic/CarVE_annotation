#!/usr/bin/env nextflow
// Phase 1 spike entry point — SpliceAI only, legacy per-gene adapter.
// Does NOT replace main.sh. See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md.
//
// Usage:
//   nextflow run main.nf --input_vcf <gene>.vcf.gz [-profile standard|local_dev] [-stub-run]
//
// --input_vcf must already be tabix-indexed alongside it (<input_vcf>.tbi), matching today's
// Step-1 (variant_converter) output shape — this spike starts downstream of that step.

include { SPLICEAI_SPIKE } from './workflows/spliceai_spike.nf'

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

    SPLICEAI_SPIKE(meta_vcf_ch)

    SPLICEAI_SPIKE.out.vcf.view { m, vcf, tbi -> "SPLICEAI_SPIKE done: ${m.partition_id} -> ${vcf}" }
}
