// Phase 3: VEP annotation, ported from src/hpc/annotate_vep_vars.sh — the first containerized
// process in this migration. main.sh already ran this inside vep.sif via a manual
// `singularity exec` call; here the container is a native Nextflow `container` directive instead,
// so the process script is a plain `vep ...` invocation with no singularity wrapper of its own —
// Nextflow's singularity integration (enabled in nextflow.config's `standard` profile) injects it.
// Runs on the same un-chunked per-gene VCF as every other predictor (main.sh:437-446), no chunking
// stage exists for VEP either. Broad-safe tier — every plugin here is one variant-level pass.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.
//
// Whole-VCF redesign (plan item 7-8): this process is also reused per-CHUNK in the new whole-VCF
// path (workflows/broad_pass_whole.nf), where every chunk of the same run shares one `meta`
// (chunk identity lives only in the task's own work directory, not in the output filename) — so
// `publishDir` below skips publishing for those chunk-level calls (`meta.partition_type ==
// 'chunk'`) via `saveAs`, since only CONCAT_CHUNKS's reassembled whole-run file
// (modules/local/chunk.nf) is a real deliverable; real per-gene calls publish as before. NOTE:
// Nextflow's `publishDir enabled:` option does NOT support a per-task closure (confirmed by
// testing directly — `enabled: { ... }` silently disables publishing unconditionally, even for a
// closure that always returns `true`) — `saveAs` returning `null` to skip a file is the correct,
// tested idiom for this. Without it, every chunk of the same run would try to publish to the
// identical filename `${meta.partition_id}.annVEP.vcf.gz` and silently overwrite each other.

process VEP_ANNOTATE {
    tag "${meta.partition_id}"
    label 'process_medium'
    container "${params.vep_sif}"
    publishDir "${params.outdir}/${meta.partition_id}/${params.out_subdir.predictors}", mode: 'copy',
        saveAs: { filename -> meta.partition_type == 'chunk' ? null : filename }

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.annVEP.vcf.gz"), path("${meta.partition_id}.annVEP.vcf.gz.tbi"), emit: vcf

    // Plugins/custom files run in the broad tier. The commented-out plugins in the bash
    // source (GeneSplicer, dbscSNV, selv, FATHMM) were already dead there and aren't ported.
    // Note: this --plugin SpliceAI is VEP's own precomputed-lookup plugin (hg38 raw score VCFs
    // under vep_dir/plugins/spliceAI), unrelated to and not a duplicate of this pipeline's own
    // live SPLICEAI_ANNOTATE process (gene-restricted tier) — main.sh already ran both together.
    //
    // Behavior change from the bash original, noted rather than hidden: main.sh submits this job
    // with `-pe smp 4` but never exports $THREADS, so the bash version's `--fork ${THREADS:-1}`
    // actually ran with --fork 1 — a latent bug, not an intentional throttle. This port uses
    // --fork ${task.cpus} (real parallelism matching the allocated slots), a deliberate behavior
    // change, not a verbatim port.
    //
    // tabix below is the container's OWN /usr/local/bin/tabix (htslib 1.9, confirmed present in
    // vep.sif), not params.tabix — the bash original indexes on the HOST after singularity exec
    // returns; this script runs entirely inside the container, so a host conda tabix binary would
    // be visible via the bind mount but not guaranteed to load (built against host libs). Using
    // the container's own tabix avoids that mismatch entirely.
    //
    // NOT ported: annotate_vep_vars.sh's `export PERL5OPT=-X` (meant to silence Mastermind
    // plugin warnings flooding stderr). Real bug found on the first live cluster run of this
    // process: `-X`, while a legal `perl` command-line switch, is NOT on PERL5OPT's own allowed-
    // switch list — Perl itself refuses it ("Illegal switch in PERL5OPT: -X."), which hard-failed
    // every chunk. Root cause of why this was never caught in the bash pipeline: `annotate_vep_
    // vars.sh` sets it as a plain host `export` before calling `singularity exec -e ...` (`-e`
    // = `--cleanenv`), and `--cleanenv` does NOT pass plain host env vars into the container —
    // only `SINGULARITYENV_`-prefixed ones are. So that `export` never actually reached Perl
    // inside the container in production either; it was always a silent no-op. Nextflow's own
    // singularity integration propagates script-exported env vars into the container more
    // faithfully than a manual bash `singularity exec` call does, which is exactly why a
    // dormant, always-broken setting only surfaces as a hard failure here. Not replaced with a
    // working equivalent — the original problem was disk-I/O noise, not correctness, and
    // reintroducing warning suppression is out of scope for this port.
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
        --plugin pLI,${params.vep_pli} \\
        --plugin dbNSFP,"${params.vep_dbnsfp},${dbnsfp_fields}" \\
        --plugin AlphaMissense,file=${params.vep_alphamissense} \\
        --plugin MaxEntScan,${params.vep_maxentscan_dir} \\
        --custom file=${params.splicevardb_vcf},short_name=splicevardb,format=vcf,type=exact,coords=0,fields=gene%hgvs%method%classification%location%doi \\
        --custom file=${params.clinvar_vcf},short_name=ClinVarVCF,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNSIGCONF%CLNREVSTAT%CLNDN%CLNVC%ALLELEID \\
        --plugin SpliceVault,file=${params.splicevault_file} \\
        --custom file=${params.gnomad_vcf},short_name=gnomADv4,format=vcf,type=exact,coords=0,fields=AF_joint%AF_grpmax_joint%grpmax_joint%fafmax_faf95_max_joint%faf95_joint%AF_exomes%AF_grpmax_exomes%fafmax_faf95_max_exomes%AF_genomes%AF_grpmax_genomes%fafmax_faf95_max_genomes%AF_joint_afr%AF_joint_amr%AF_joint_asj%AF_joint_eas%AF_joint_fin%AF_joint_mid%AF_joint_nfe%AF_joint_sas%AF_joint_ami%AF_joint_remaining%faf95_joint_afr%faf95_joint_amr%faf95_joint_eas%faf95_joint_mid%faf95_joint_nfe%faf95_joint_sas%nhomalt_joint%exomes_filters%genomes_filters%AN_joint%AN_joint_afr%AN_joint_ami%AN_joint_amr%AN_joint_asj%AN_joint_eas%AN_joint_fin%AN_joint_mid%AN_joint_nfe%AN_joint_raw%AN_joint_remaining%AN_joint_sas%AN_grpmax_joint%AN_genomes%AN_exomes%AC_joint%AC_joint_afr%AC_joint_ami%AC_joint_amr%AC_joint_asj%AC_joint_eas%AC_joint_fin%AC_joint_mid%AC_joint_nfe%AC_joint_raw%AC_joint_remaining%AC_joint_sas%AC_grpmax_joint%AC_genomes%AC_exomes%nhomalt_joint_afr%nhomalt_joint_ami%nhomalt_joint_amr%nhomalt_joint_asj%nhomalt_joint_eas%nhomalt_joint_fin%nhomalt_joint_mid%nhomalt_joint_nfe%nhomalt_joint_raw%nhomalt_joint_remaining%nhomalt_joint_sas%nhomalt_grpmax_joint%nhomalt_exomes%nhomalt_genomes

    tabix ${meta.partition_id}.annVEP.vcf.gz
    """

    // Real inference needs the full VEP cache/plugin data tree under /references (many GB,
    // several plugin files not mirrored to this sandbox's local_references at all — unlike
    // SpliceAI's single FASTA), and local_dev leaves singularity.enabled=false, so `container`
    // is never exercised here either. -stub-run skips this process's `script:` block entirely
    // (Nextflow does not execute it in stub mode) — the stub proves channel/DAG wiring only,
    // nothing about the command line above or the container directive; both are unverified in
    // this sandbox and must be checked on the cluster.
    stub:
    """
    cp ${vcf} ${meta.partition_id}.annVEP.vcf.gz
    cp ${tbi} ${meta.partition_id}.annVEP.vcf.gz.tbi
    """
}
