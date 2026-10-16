// VEP 111 with its plugins and custom annotations, run in the VEP container (broad tier).

process VEP_ANNOTATE {
    tag "${meta.partition_id}"
    label 'process_medium'
    container "${params.vep_sif}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annVEP.vcf.gz"), path("${meta.partition_id}.annVEP.vcf.gz.tbi"), emit: vcf

    // `--plugin SpliceAI` is VEP's precomputed-score lookup, separate from SPLICEAI_ANNOTATE.
    // `tabix` is the container's own binary. Do not set PERL5OPT=-X: Perl rejects it.
    script:
    def vep_files = params.vep_dir
    def dbnsfp_fields = file(params.vep_dbnsfp_fields).readLines()
        .findAll { it.trim() && !it.trim().startsWith('#') }
        .collect { it.trim() }.join(',')
    """
    vep --fork ${task.cpus} -species homo_sapiens \\
        --dir_plugins ${params.data_lab_pgp}/resources/annotation/VEP/VEP_plugins \\
        --offline --cache --cache_version ${params.vep_cache_version} --dir ${vep_files}/cache \\
        --synonyms ${vep_files}/cache/homo_sapiens/${params.vep_cache_version}_GRCh38/chr_synonyms.txt \\
        --fasta ${params.fasta} --assembly GRCh38 \\
        -i ${vcf} \\
        -o ${meta.partition_id}.annVEP.vcf.gz \\
        --vcf --compress_output bgzip --force_overwrite \\
        --check_existing \\
        --sift b --polyphen b --ccds --hgvs --symbol --numbers --domains \\
        --regulatory --canonical --protein --biotype --af --af_1kg --pubmed \\
        --uniprot --mane --tsl --appris --variant_class --gene_phenotype --mirna \\
        --no_stats \\
        --plugin CADD,"${params.vep_cadd}" \\
        --plugin UTRAnnotator,file=${params.vep_utrannotator} \\
        --plugin SpliceAI,snv=${params.vep_spliceai_snv},indel=${params.vep_spliceai_indel} \\
        --plugin SpliceRegion \\
        --plugin REVEL,${params.vep_revel} \\
        --plugin Mastermind,${params.vep_mastermind} \\
        --plugin LoFtool,${params.vep_loftool} \\
        --plugin NMD \\
        --plugin pLI,${params.vep_pli} \\
        --plugin dbNSFP,"${params.vep_dbnsfp},${dbnsfp_fields}" \\
        --plugin AlphaMissense,file=${params.vep_alphamissense} \\
        --plugin MaxEntScan,${params.vep_maxentscan_dir} \\
        --custom file=${params.splicevardb_vcf},short_name=splicevardb,format=vcf,type=exact,coords=0,fields=gene%hgvs%method%classification%location%doi \\
        --custom file=${params.clinvar_vcf},short_name=ClinVarVCF,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNSIGCONF%CLNREVSTAT%CLNDN%CLNVC%ALLELEID \\
        --plugin SpliceVault,file=${params.splicevault_file} \\
        --custom file=${params.vep_remm},short_name=ReMM,format=bigwig,type=exact \\
        --custom file=${params.gnomad_vcf},short_name=gnomADv4,format=vcf,type=exact,coords=0,fields=AF_joint%AF_grpmax_joint%grpmax_joint%fafmax_faf95_max_joint%faf95_joint%AF_exomes%AF_grpmax_exomes%fafmax_faf95_max_exomes%AF_genomes%AF_grpmax_genomes%fafmax_faf95_max_genomes%AF_joint_afr%AF_joint_amr%AF_joint_asj%AF_joint_eas%AF_joint_fin%AF_joint_mid%AF_joint_nfe%AF_joint_sas%AF_joint_ami%AF_joint_remaining%faf95_joint_afr%faf95_joint_amr%faf95_joint_eas%faf95_joint_mid%faf95_joint_nfe%faf95_joint_sas%nhomalt_joint%exomes_filters%genomes_filters%AN_joint%AN_joint_afr%AN_joint_ami%AN_joint_amr%AN_joint_asj%AN_joint_eas%AN_joint_fin%AN_joint_mid%AN_joint_nfe%AN_joint_raw%AN_joint_remaining%AN_joint_sas%AN_grpmax_joint%AN_genomes%AN_exomes%AC_joint%AC_joint_afr%AC_joint_ami%AC_joint_amr%AC_joint_asj%AC_joint_eas%AC_joint_fin%AC_joint_mid%AC_joint_nfe%AC_joint_raw%AC_joint_remaining%AC_joint_sas%AC_grpmax_joint%AC_genomes%AC_exomes%nhomalt_joint_afr%nhomalt_joint_ami%nhomalt_joint_amr%nhomalt_joint_asj%nhomalt_joint_eas%nhomalt_joint_fin%nhomalt_joint_mid%nhomalt_joint_nfe%nhomalt_joint_raw%nhomalt_joint_remaining%nhomalt_joint_sas%nhomalt_grpmax_joint%nhomalt_exomes%nhomalt_genomes

    tabix ${meta.partition_id}.annVEP.vcf.gz
    """

    stub:
    """
    cp ${vcf} ${meta.partition_id}.annVEP.vcf.gz
    cp ${tbi} ${meta.partition_id}.annVEP.vcf.gz.tbi
    """
}
