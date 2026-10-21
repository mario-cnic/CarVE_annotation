// Annotated VCF -> final table. VCF_TO_TSV flattens INFO/CSQ into a TSV; TSV_TO_TABLE applies status
// contracts and writes the table. A retry of TSV_TO_TABLE reuses the cached TSV. Neither publishes:
// TAG_TRANSCRIPT_PRIORITY publishes the final table.

process VCF_TO_TSV {
    tag "${meta.partition_id}"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("scratch.tsv"), emit: tsv

    script:
    """
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export ARROW_IO_THREADS=1
    export PYTHONUNBUFFERED=1

    PYTHONPATH="${params.vcf_parser_pythonpath}" ${params.vcf_parser_python} ${params.vcf_parser_script} \\
        --input ${vcf} \\
        --output scratch.tsv \\
        --vep_columns ${params.vep_cols_file} \\
        --add_info --add_vep --overwrite \\
        --logging_level INFO
    """

    stub:
    """
    touch scratch.tsv
    """
}

process TSV_TO_TABLE {
    tag "${meta.partition_id}"
    label 'process_medium'

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("${meta.partition_id}.parsed.clean.${params.output_format}"), emit: table

    script:
    """
    export OMP_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export ARROW_IO_THREADS=1
    export PYTHONUNBUFFERED=1

    ${params.datasci_python} ${params.filter_variants_script} \\
        --input ${tsv} \\
        --output ${meta.partition_id}.parsed.clean.${params.output_format} \\
        --logging_level INFO \\
        --hgnc-table ${params.hgnc_table} \\
        --disease-gene-list ${params.disease_gene_list} \\
        --spliceai-symbol-map ${params.spliceai_symbol_map} \\
        --enst-spip-map ${params.enst_spip_map} \\
        --dbnsfp-aligned-columns ${params.dbnsfp_aligned_columns} \\
        ${params.table_chunk_rows ? "--chunk-rows ${params.table_chunk_rows}" : ''} \\
        --skip_quality_filter
    """

    stub:
    """
    touch ${meta.partition_id}.parsed.clean.${params.output_format}
    """
}

// Annotated VCF -> table; emits `table` as tuple(meta, table file).
workflow VCF_TO_TABLE {
    take:
    vcf_ch   // tuple(meta, vcf, tbi)

    main:
    VCF_TO_TSV(vcf_ch)
    TSV_TO_TABLE(VCF_TO_TSV.out.tsv)

    emit:
    table = TSV_TO_TABLE.out.table
}
