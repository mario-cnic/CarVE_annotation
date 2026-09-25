// Subsets the full input VCF to the gene-restricted tier's region-of-interest list, once, before
// chunking (not per-gene — doesn't scale to a ~5,742-region BED).

process GENE_SUBSET {
    tag "${meta.partition_id}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)
    // Declared `path` input, not a bare params interpolation: this BED is actively-maintained
    // pipeline data, not an immutable reference, and only a declared path input keeps future
    // revisions visible to `-resume`.
    path(gene_bed)

    output:
    tuple val(meta), path("subset.vcf.gz"), path("subset.vcf.gz.tbi"), emit: vcf

    script:
    """
    ${params.bcftools} view -R ${gene_bed} ${vcf} -Oz -o subset.vcf.gz
    ${params.tabix} -p vcf subset.vcf.gz
    """
}
