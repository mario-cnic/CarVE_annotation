// Whole-VCF redesign (plan item 7): subsets the full input VCF to the gene-restricted tier's
// region-of-interest list, ONCE, before chunking — not per-gene (see the plan's Context section,
// "Course correction": an earlier draft of this design ran bcftools per gene, which does not
// scale to Complete_gene_list_V5's ~5,742-region BED). No stub: this is pure bcftools, no model
// dependency, genuinely verifiable in this sandbox.
// See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md for scope/non-goals.

process GENE_SUBSET {
    tag "${meta.partition_id}"
    label 'process_low'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("subset.vcf.gz"), path("subset.vcf.gz.tbi"), emit: vcf

    script:
    """
    ${params.bcftools} view -R ${params.gene_restriction_bed} ${vcf} -Oz -o subset.vcf.gz
    ${params.tabix} -p vcf subset.vcf.gz
    """
}
