// Merges the five predictor VCFs into one annotated VCF. VEP is the record set; SPiP, Pangolin,
// SpliceAI and Branchpointer INFO tags are transferred onto it in that order. Any failure fails the task.

process MERGE_ANNOTATIONS {
    tag "${meta.partition_id}"
    label 'process_single'
    publishDir "${params.outdir}/${meta.partition_id}/${params.out_subdir.annotation}", mode: 'copy'

    input:
    tuple val(meta), path(vep_vcf), path(vep_tbi), path(bp_vcf), path(bp_tbi), path(pango_vcf), path(pango_tbi), path(sai_vcf), path(sai_tbi), path(spip_vcf), path(spip_tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annotated.vcf.gz"), path("${meta.partition_id}.annotated.vcf.gz.tbi"), emit: vcf

    // Nextflow already runs scripts with -ue; pipefail is added here.
    script:
    """
    set -o pipefail

    # Intermediates are fast-compressed BCF indexed on write (annotate -a needs an indexed target).
    ${params.bcftools} annotate -a ${spip_vcf} -c SPiP ${vep_vcf} -Ob1 -W -o tmp_spip.bcf

    ${params.bcftools} annotate -a ${pango_vcf} -c Pangolin,PangolinTissue tmp_spip.bcf -Ob1 -W -o tmp_pangolin.bcf

    ${params.bcftools} annotate -a ${sai_vcf} -c SpliceAI tmp_pangolin.bcf -Ob1 -W -o tmp_spliceai.bcf

    ${params.bcftools} annotate -a ${bp_vcf} -c Branchpointer_prob,Branchpointer_U2_energy,Branchpoint_disrupted,LaBranchoR_score,LaBranchoR_acc_dist tmp_spliceai.bcf -Oz -o ${meta.partition_id}.annotated.vcf.gz
    ${params.tabix} -f ${meta.partition_id}.annotated.vcf.gz
    """

    stub:
    """
    cp ${vep_vcf} ${meta.partition_id}.annotated.vcf.gz
    cp ${vep_tbi} ${meta.partition_id}.annotated.vcf.gz.tbi
    """
}
