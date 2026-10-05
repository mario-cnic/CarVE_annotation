// Pre-flight input check for the whole-VCF entry point: fails the run before any predictor is
// scheduled if the input VCF's ##contig lengths or REF alleles disagree with params.fasta, if it
// has multiallelic or non-normalised records, is a gVCF, or has duplicate sample names. The report
// also records the genotype mode and sample names. Logic and thresholds:
// src/python/check_vcf_assembly.py.
//
// No stub block, so the check also runs for real under -stub-run: it reads the whole VCF twice
// (bcftools norm and a pysam pass; minutes on a WGS input) and the FASTA.

process CHECK_INPUT_ASSEMBLY {
    tag "${meta.partition_id}"
    label 'process_single'
    // A failed check is a property of the input, not a transient error; retrying can't fix it.
    errorStrategy 'terminate'
    // Only the report; the VCF passes through unchanged and must not be copied.
    publishDir "${params.outdir}/${meta.partition_id}/${params.out_subdir.reports}", mode: 'copy', pattern: '*.assembly_check.tsv'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path(vcf), path(tbi), emit: vcf
    path "${meta.partition_id}.assembly_check.tsv", emit: report

    script:
    """
    ${params.check_assembly_python} ${params.check_assembly_script} \\
        --vcf ${vcf} \\
        --fasta ${params.fasta} \\
        --report ${meta.partition_id}.assembly_check.tsv \\
        --ref-check-records ${params.assembly_ref_check_records} \\
        --max-ref-mismatch-frac ${params.assembly_max_ref_mismatch_frac} \\
        --min-ref-checked ${params.assembly_min_ref_checked} \\
        --max-multiallelic ${params.input_max_multiallelic} \\
        --max-unnormalised ${params.input_max_unnormalised} \\
        --bcftools ${params.bcftools}
    """
}
