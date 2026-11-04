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

    // `view -R` writes records in the BED's contig order; `sort` restores the header's order.
    // Sort temp files go in the task directory, not the node's TMPDIR.
    script:
    """
    set -o pipefail
    ${params.bcftools} view -R ${gene_bed} ${vcf} -Ou \\
        | ${params.bcftools} sort -T ./bcftools-sort.XXXXXX -Oz -o subset.vcf.gz -
    ${params.tabix} -p vcf subset.vcf.gz
    """
}
