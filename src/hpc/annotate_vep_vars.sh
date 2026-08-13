#!/bin/bash
#$ -P PGP
#$ -N VEP_ANNOTATION
#$ -A PGP
#$ -l thread=2
#$ -l h_vmem=10G
#$ -o _log/vep_annotation.stdout
#$ -e _log/vep_annotation.stderr

in_file=$1
out_file=$2

vep_files=/references/genomes/Homo_sapiens/annotations/VEP/GRCh38
plugins=${vep_files}/plugins
gatk_bundle=/references/genomes/Homo_sapiens/GATK_bundle
# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "vep" "$in_file"
    TIME_LOG="$LOG_DIR/vep.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
if ! command -v run_command_timed &>/dev/null; then
    run_command_timed() { "$@"; }
fi
# ------------------------------------------------

# Silence all Perl interpreter warnings to prevent millions of Mastermind plugin
# non-numeric warnings from flooding stderr and causing disk I/O bottlenecks.
export PERL5OPT=-X

run_command_timed singularity exec \
    -e --bind /home/${USER}:/home/${USER} -B /data_tmp:/tmp -B /references:/references -B /data_lab_PGP:/data_lab_PGP \
    /data_lab_PGP/resources/sif_images/vep.sif \
    vep --fork ${THREADS:-1} -species homo_sapiens \
    --dir_plugins /data_lab_PGP/resources/annotation/VEP/VEP_plugins \
    --offline --cache --cache_version 111 --dir ${vep_files}/cache \
    --fasta ${gatk_bundle}/v0/Homo_sapiens_assembly38.fasta --assembly GRCh38 \
    -i ${in_file} \
    -o ${out_file} \
    --vcf --compress_output bgzip --force_overwrite \
    --check_existing \
    --everything \
    --no_stats \
    --plugin CADD,"${plugins}/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz" \
    --plugin UTRAnnotator,file=${plugins}/UTRAnnotator/uORF_5UTR_GRCh38_PUBLIC.txt \
    --plugin SpliceAI,snv=${plugins}/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=${plugins}/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
    --plugin SpliceRegion \
    --plugin REVEL,${plugins}/REVEL/GRCh38/new_tabbed_revel_grch38.tsv.gz \
    --plugin Mastermind,${plugins}/MasterMind/mastermind_cited_variants_reference-2024.07.03-grch38.vcf.gz \
    --plugin LoFtool,${plugins}/LoFtool/LoFtool_scores.txt \
    --plugin pLI,${plugins}/pLI/plI_gene.txt \
    --plugin dbNSFP,"${plugins}/dbNSFP/dbNSFP4.8a_grch38.gz,ALL" \
    --plugin AlphaMissense,file=${plugins}/AlphaMissense/AlphaMissense_hg38.tsv.gz \
    --plugin MaxEntScan,${plugins}/MaxEntScan/download/fordownload \
    --custom file=/references/genomes/Homo_sapiens/annotations/splicevardb/splicevardb.to.annotate.tsv.vcf.gz,short_name=splicevardb,format=vcf,type=exact,coords=0,fields=gene%hgvs%method%classification%location%doi \
    --plugin SpliceVault,file=/data_lab_PGP/resources/annotation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz \
    --custom file=${vep_files}/../../gnomAD/GRCh38/v4.1/gnomAD.v4.1.vcf.gz,short_name=gnomADv4,format=vcf,type=exact,coords=0,fields=AF_joint%AF_grpmax_joint%grpmax_joint%fafmax_faf95_max_joint%faf95_joint%AF_exomes%AF_grpmax_exomes%fafmax_faf95_max_exomes%AF_genomes%AF_grpmax_genomes%fafmax_faf95_max_genomes%AF_joint_afr%AF_joint_amr%AF_joint_asj%AF_joint_eas%AF_joint_fin%AF_joint_mid%AF_joint_nfe%AF_joint_sas%AF_joint_ami%AF_joint_remaining%faf95_joint_afr%faf95_joint_amr%faf95_joint_eas%faf95_joint_mid%faf95_joint_nfe%faf95_joint_sas%nhomalt_joint%exomes_filters%genomes_filters%AN_joint%AN_joint_afr%AN_joint_ami%AN_joint_amr%AN_joint_asj%AN_joint_eas%AN_joint_fin%AN_joint_mid%AN_joint_nfe%AN_joint_raw%AN_joint_remaining%AN_joint_sas%AN_grpmax_joint%AN_genomes%AN_exomes%AC_joint%AC_joint_afr%AC_joint_ami%AC_joint_amr%AC_joint_asj%AC_joint_eas%AC_joint_fin%AC_joint_mid%AC_joint_nfe%AC_joint_raw%AC_joint_remaining%AC_joint_sas%AC_grpmax_joint%AC_genomes%AC_exomes%nhomalt_joint%nhomalt_joint_afr%nhomalt_joint_ami%nhomalt_joint_amr%nhomalt_joint_asj%nhomalt_joint_eas%nhomalt_joint_fin%nhomalt_joint_mid%nhomalt_joint_nfe%nhomalt_joint_raw%nhomalt_joint_remaining%nhomalt_joint_sas%nhomalt_grpmax_joint%age_hist_hom_bin_freq_joint%age_hist_hom_n_smaller_joint%age_hist_hom_n_larger_joint%nhomalt_exomes%nhomalt_genomes%exomes_filters%genomes_filters \
    # --per_gene \
    #--plugin GeneSplicer,/data3/genomes/homo_sapiens/annotations/GeneSplicer/bin/linux/genesplicer,/data3/genomes/homo_sapiens/annotations/GeneSplicer/human \
    #--plugin dbscSNV,/data3/genomes/Homo_sapiens/annotations/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
    #--custom file=/data3/genomes/homo_sapiens/annotations/selv/selv.scores.vcf.gz,short_name=selv,format=vcf,type=exact,coords=0,fields=%SELV_score%SELV_weight_score \
    # --gene_phenotype \ # ya en --everything
    # --plugin FATHMM,"python /data3/genomes/Homo_sapiens/annotations/FATHMM/fathmm.py" \ # lento y no informativo
    
CMD_EXIT_CODE=$?

# index vep file after annotation is finished
if [ $CMD_EXIT_CODE -eq 0 ]; then
    /data_lab_PGP/shared/utils/conda_envs/genomics/bin/tabix ${out_file}
fi

echo "Variants saved at ${out_file} at $(date)"

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    log_step_end "vep" "$out_file" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------
