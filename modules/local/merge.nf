// Merges the six predictor VCFs into one annotated VCF. VEP is the record set; SPiP, Pangolin,
// SpliceAI, LaBranchoR and SpliceAI2 INFO tags are transferred onto it in that order. Any failure fails the task.

process MERGE_ANNOTATIONS {
    tag "${meta.partition_id}"
    label 'process_single'
    publishDir "${params.outdir}/${meta.partition_id}/${params.out_subdir.annotation}", mode: 'copy'

    input:
    tuple val(meta), path(vep_vcf), path(vep_tbi), path(bp_vcf), path(bp_tbi), path(sai2_vcf), path(sai2_tbi), path(pango_vcf), path(pango_tbi), path(sai_vcf), path(sai_tbi), path(spip_vcf), path(spip_tbi)

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

    ${params.bcftools} annotate -a ${bp_vcf} -c LaBranchoR_score,LaBranchoR_acc_dist,LaBranchoR_ISM_d_best,LaBranchoR_ISM_d_max,LaBranchoR_ISM_dist,LaBranchoR_ISM_n_3pss tmp_spliceai.bcf -Ob1 -W -o tmp_branchpoint.bcf

    ${params.bcftools} annotate -a ${sai2_vcf} -c SpliceAI2_score,SpliceAI2_DS_AG,SpliceAI2_DS_AL,SpliceAI2_DS_DG,SpliceAI2_DS_DL,SpliceAI2_DP_AG,SpliceAI2_DP_AL,SpliceAI2_DP_DG,SpliceAI2_DP_DL,SpliceAI2_strand tmp_branchpoint.bcf -Oz -o ${meta.partition_id}.annotated.vcf.gz
    ${params.tabix} -f ${meta.partition_id}.annotated.vcf.gz

    # The merged VCF keeps exactly VEP's CHROM/POS/REF/ALT records (MISC-6).
    site_keys() { ${params.bcftools} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' "\$1" | LC_ALL=C sort -S 1G -T . | sha256sum | cut -d' ' -f1; }
    if [ "\$(site_keys ${vep_vcf})" != "\$(site_keys ${meta.partition_id}.annotated.vcf.gz)" ]; then
        echo "MERGE_ANNOTATIONS: CHROM/POS/REF/ALT of the merged VCF differ from VEP's" >&2
        exit 1
    fi
    """

    stub:
    """
    cp ${vep_vcf} ${meta.partition_id}.annotated.vcf.gz
    cp ${vep_tbi} ${meta.partition_id}.annotated.vcf.gz.tbi
    """
}
