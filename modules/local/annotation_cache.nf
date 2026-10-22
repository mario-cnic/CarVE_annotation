// Predictor score cache (src/python/annotation_cache.py): splits the gene-restricted input into
// sites already scored for this predictor version and sites still to compute, and adds the newly
// computed sites to the cache as one segment.

process CACHE_SPLIT {
    tag "${meta.partition_id}/${predictor}"
    label 'process_low'
    publishDir "${params.pipeline_info_dir}", mode: 'copy', pattern: '*.cache_split.json'

    input:
    tuple val(meta), path(vcf), path(tbi)
    val(predictor)
    // Every file of the version's segments/ folder; a path input, so new segments invalidate -resume.
    path(segments, stageAs: 'segments/*')

    output:
    tuple val(meta), path("${predictor}.hits.vcf.gz"), emit: hits
    tuple val(meta), path("${predictor}.misses.vcf.gz"), path("${predictor}.misses.vcf.gz.tbi"), env(N_MISSES), emit: misses
    path("${meta.partition_id}.${predictor}.cache_split.json"), emit: report

    // No stub: the split also runs under -stub-run.
    script:
    """
    ${params.python_spliceai} ${params.annotation_cache_script} split \\
        --predictor ${predictor} --input ${vcf} --segments-dir segments \\
        --hits ${predictor}.hits.vcf.gz --misses ${predictor}.misses.vcf.gz \\
        --report ${meta.partition_id}.${predictor}.cache_split.json \\
        --bcftools ${params.bcftools} --tabix ${params.tabix} --bgzip ${params.bgzip}
    N_MISSES=\$(${params.bcftools} index -n ${predictor}.misses.vcf.gz)
    """
}

process CACHE_WRITE {
    tag "${meta.partition_id}/${predictor}"
    label 'process_low'
    // Segments are append-only; an existing file of the same name is never replaced.
    publishDir "${params.annotation_cache_dir}/${predictor}/${version_id}/segments", mode: 'copy', overwrite: false,
        saveAs: { fn -> fn.tokenize('/')[-1] }

    input:
    tuple val(meta), path(annotated), path(annotated_tbi), path(misses), path(misses_tbi)
    val(predictor)
    val(version_id)

    output:
    path("segment/*")

    // The UTC stamp is fixed when the task runs, so a resumed (cached) task republishes the same names.
    script:
    """
    mkdir segment
    ${params.python_spliceai} ${params.annotation_cache_script} make-segment \\
        --predictor ${predictor} --annotated ${annotated} --misses ${misses} \\
        --output-prefix segment/${meta.run_id}.\$(date -u +%Y%m%dT%H%M%SZ) \\
        --run-id ${meta.run_id} --version-id ${version_id} \\
        --bcftools ${params.bcftools} --tabix ${params.tabix} --bgzip ${params.bgzip}
    """
}
