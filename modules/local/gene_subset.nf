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
    path(gene_bed)

    output:
    tuple val(meta), path("subset.vcf.gz"), path("subset.vcf.gz.tbi"), emit: vcf

    script:
    // gene_bed is a real Nextflow `path` input (not a bare params.gene_restriction_bed string
    // interpolation, unlike this codebase's other static-reference-file processes) deliberately:
    // this file is actively-maintained pipeline data, not an immutable external reference (it was
    // revised mid-migration, BUG_TRACKER.md MISC-10, after silently matching ~0 records against a
    // real chr-prefixed WGS input) — a bare params interpolation only tracks the file's PATH in
    // the task's resume-cache signature, not its CONTENT, so a future BED fix would otherwise be
    // silently invisible to `-resume` and a stale subset.vcf.gz would be reused from cache.
    """
    ${params.bcftools} view -R ${gene_bed} ${vcf} -Oz -o subset.vcf.gz
    ${params.tabix} -p vcf subset.vcf.gz
    """
}
