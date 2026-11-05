#!/usr/bin/env nextflow
// Whole-VCF entry point: one GRCh38 VCF of any size in, all predictors run, one final table out.
//
// Usage:
//   nextflow run annotate_vcf.nf --input_vcf <any>.vcf.gz [-profile standard|local_dev] [-stub-run]
//
// --input_vcf must already be tabix-indexed alongside it (<input_vcf>.tbi).

import java.nio.file.Files
import java.security.MessageDigest

// Nextflow does not see the content of scripts passed by path, so -resume would reuse a task whose code
// changed. params.code_md5.<key> is the md5 of the code (and small files read by path) of each process;
// the process names it on a `# code` line of its script, which puts that value in the task hash. A
// directory is walked recursively (no __pycache__); a missing path counts as "missing:<path>".
def codeMd5(List paths) {
    def md = MessageDigest.getInstance('MD5')
    paths.each { p ->
        def root = file(p.toString())
        if (!root.exists()) {
            md.update("missing:${p}".toString().bytes)
            return
        }
        def files = root.isDirectory()
            ? Files.walk(root).toList().findAll { Files.isRegularFile(it) && !it.toString().contains('__pycache__') && !it.toString().endsWith('.pyc') }.sort { it.toString() }
            : [root]
        files.each { f ->
            md.update(root.relativize(f).toString().bytes)
            md.update(Files.readAllBytes(f))
        }
    }
    md.digest().encodeHex().toString()
}
def stampPy = "${projectDir}/src/python/provenance_stamp.py"
params.code_md5 = [
    check       : codeMd5([params.check_assembly_script]),
    chunk       : codeMd5([params.split_vcf_chunks_script]),
    spliceai    : codeMd5([params.annotate_spliceai_script]),
    pangolin    : codeMd5(["${params.pangolin_repo}/pangolin"]),
    spip        : codeMd5([params.spip_script]),
    branchpoint : codeMd5([params.branchpointer_script]),
    vep         : codeMd5([params.vep_dbnsfp_fields]),
    vcf_to_tsv  : codeMd5([params.vcf_parser_script, "${file(params.vcf_parser_script).parent}/logger.py", params.vep_cols_file]),
    tsv_to_table: codeMd5([file(params.filter_variants_script).parent, params.spliceai_symbol_map, params.enst_spip_map,
                           params.dbnsfp_aligned_columns, params.gene_constraint_table, params.disease_gene_list, params.hgnc_table]),
    tag         : codeMd5([params.gene_transcript_priority_script, stampPy, "${projectDir}/resources/gene_transcript_mapping.txt"]),
    genotypes   : codeMd5([params.extract_genotypes_script, stampPy]),
    wide        : codeMd5([params.genotypes_wide_script, stampPy]),
    cache       : codeMd5([params.annotation_cache_script]),
]

include { CHECK_INPUT_ASSEMBLY } from './modules/local/check_assembly.nf'
include { BROAD_PASS_WHOLE_SUBWORKFLOW } from './workflows/broad_pass_whole.nf'
include { GENE_RESTRICTED_WHOLE_SUBWORKFLOW } from './workflows/gene_restricted_whole.nf'
include { MERGE_SUBWORKFLOW } from './workflows/merge_subworkflow.nf'
include { VCF_TO_TABLE } from './modules/local/vcf_to_table.nf'
include { TAG_TRANSCRIPT_PRIORITY } from './modules/local/tag_transcript_priority.nf'
include { EXTRACT_GENOTYPES; ADD_WIDE_GENOTYPES } from './modules/local/genotypes.nf'

workflow {
    // run_annotate_vcf.sh reads the params this exact command line resolves to (profile, -c files, CLI
    // values and the config values derived from them) from here, then launches the real run.
    if (params.dump_params) {
        file(params.dump_params).text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))
        return
    }
    if (!params.input_vcf) {
        error "Missing required --input_vcf <file>.vcf.gz (must have a sibling .tbi index)"
    }
    if (!(params.spip_tier in ['restricted', 'broad'])) {
        error "params.spip_tier must be 'restricted' or 'broad', got: ${params.spip_tier}"
    }
    if (params.annotation_cache_dir && !(params.annotation_cache_version_spliceai && params.annotation_cache_version_pangolin
            && params.annotation_cache_segments_spliceai && params.annotation_cache_segments_pangolin)) {
        error "--annotation_cache_dir needs the cache version ids and segment lists; launch with run_annotate_vcf.sh, which derives them from the run manifest"
    }

    input_vcf_file = file(params.input_vcf)
    input_tbi_file = file("${params.input_vcf}.tbi")
    if (!input_vcf_file.exists()) { error "Input VCF not found: ${input_vcf_file}" }
    if (!input_tbi_file.exists()) { error "Missing tabix index: ${input_tbi_file}" }

    // run_id defaults to the input file name, so re-runs of one input share an output folder.
    run_id = params.run_id ?: input_vcf_file.baseName.replaceAll(/\.vcf$/, '')
    meta = [
        partition_type : 'chunk',
        partition_id   : run_id,
        build           : 'GRCh38',
        run_id          : run_id
    ]

    // Every predictor consumes CHECK_INPUT_ASSEMBLY's output, so nothing is scheduled until the
    // input's contig lengths and REF alleles are confirmed against params.fasta.
    CHECK_INPUT_ASSEMBLY(Channel.of(tuple(meta, input_vcf_file, input_tbi_file)))
    meta_vcf_ch = CHECK_INPUT_ASSEMBLY.out.vcf

    BROAD_PASS_WHOLE_SUBWORKFLOW(meta_vcf_ch)
    GENE_RESTRICTED_WHOLE_SUBWORKFLOW(meta_vcf_ch)

    // Only the tier selected by spip_tier emits SPiP output.
    spip_merged_ch = GENE_RESTRICTED_WHOLE_SUBWORKFLOW.out.spip.mix(BROAD_PASS_WHOLE_SUBWORKFLOW.out.spip)

    MERGE_SUBWORKFLOW(
        BROAD_PASS_WHOLE_SUBWORKFLOW.out.vep,
        BROAD_PASS_WHOLE_SUBWORKFLOW.out.branchpoint,
        GENE_RESTRICTED_WHOLE_SUBWORKFLOW.out.pangolin,
        GENE_RESTRICTED_WHOLE_SUBWORKFLOW.out.spliceai,
        spip_merged_ch
    )

    MERGE_SUBWORKFLOW.out.vcf.view { m, vcf, tbi -> "MERGE done: ${m.partition_id} -> ${vcf}" }

    // Every VEP transcript row is kept; TAG_TRANSCRIPT_PRIORITY tags it and publishes the final table.
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
