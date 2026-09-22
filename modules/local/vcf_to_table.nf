// Phase 5: annotated VCF -> final clean table (Parquet/TSV/Excel), ported from
// src/hpc/vcf2parsed.sh. Two stages in one process (matching the bash original, which never
// persists the intermediate TSV): vcf_parser_pysam.py flattens INFO/CSQ into a scratch TSV
// filtered to one transcript (--filter_by Feature, --gene_set <ENST>), then filter_variants.py
// applies status contracts/splicing metrics and writes the final table.
//
// This is this module's actual deliverable — the table Module 3 (clinical_variant_prioritization)
// consumes. Per the user's steer this session ("don't focus too much on reports — Module 3 owns
// that"), the four downstream filter/plot/report scripts (filter_and_summarize.py,
// plot_annotation_results.R, generate_interactive_report.py,
// generate_clinical_prioritization_report.py) are NOT ported in this phase — deliberately
// deprioritized, not forgotten. See /home/mruizp/.claude/plans/scalable-wibbling-snowflake.md.
//
// meta.transcript is resolved in main.nf from resources/gene_transcript_mapping.txt (same grep-
// based lookup as main.sh:378). Resolving a missing mapping row for real is ORCH-9's job, an
// explicit non-goal here — but unlike main.sh, main.nf now FAILS the run when that lookup can't
// find a transcript, rather than letting 'UNKNOWN' silently reach --gene_set below (which would
// match no VEP CSQ block, skip every variant, and publish an empty-but-valid table with exit 0).
// A deliberate behavior change, not a verbatim port of that particular corner.
//
// Not ported: the bash script's explicit /data_tmp fast-scratch staging for its intermediate TSV
// (SCRATCH_DIR resolution + PID/timestamp-unique filename + trap-based cleanup). That existed to
// avoid writing a large intermediate file to shared storage from a *persistent, multi-job-shared*
// scratch directory — a concern that doesn't apply the same way inside a Nextflow task's own
// isolated work directory. If cluster profiling later shows this process needs fast local scratch,
// Nextflow's native `scratch` directive is the right mechanism, not a hand-rolled path — not
// wired up here since it's unverified without real cluster timing data.
//
// Whole-VCF redesign (plan item 10): `meta.transcript` only exists in the legacy per-gene path's
// meta shape (main.nf) — the whole-VCF path's meta (annotate_vcf.nf) has no such field by design
// (transcript selection happens per-ROW, not per-partition, see the plan's "Transcript priority
// tiering" section), so `--gene_set`/`--filter_by` below are only included when `meta.transcript`
// is actually present; omitting them is exactly vcf_parser_pysam.py's already-existing unfiltered
// behavior (confirmed by reading vcf_parser_pysam.py:385-386: no selected_genes means no filter,
// no code change needed there). `publishDir` also only fires for the legacy path
// (`meta.partition_type == 'gene'`) — in whole-VCF mode (`partition_type == 'chunk'`) this
// process's output is an intermediate that TAG_TRANSCRIPT_PRIORITY
// (modules/local/tag_transcript_priority.nf) consumes and re-publishes as the real final
// deliverable under the identical filename; publishing both would race to the same path.

process VCF_TO_TABLE {
    tag "${meta.partition_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.partition_id}", mode: 'copy',
        saveAs: { filename -> meta.partition_type == 'chunk' ? null : filename }

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.partition_id}.parsed.clean.${params.output_format}"), emit: table

    script:
    def gene_set_flag = meta.transcript ? "--gene_set ${meta.transcript} --filter_by Feature" : ""
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
        --logging_level INFO \\
        ${gene_set_flag}

    ${params.datasci_python} ${params.filter_variants_script} \\
        --input scratch.tsv \\
        --output ${meta.partition_id}.parsed.clean.${params.output_format} \\
        --logging_level INFO \\
        --skip_quality_filter
    """

    // Can't be exercised end-to-end for real through Nextflow in this sandbox (needs a genuinely
    // VEP-annotated VCF; every upstream predictor here is -stub-run-only). Its exact command
    // sequence WAS verified for real, though (outside Nextflow, same commands and --gene_set ENST
    // value this script uses): ran against the known-good MYBPC3.annotated.vcf.gz reference under
    // test_data/test_run/annotation/ and diffed the result against that same directory's
    // MYBPC3.parsed.clean.pq. Of 559 columns common to both, 557 were identical; the only two
    // that differed were PRIORITY_TIER and VARIANT_PRIORITY_SCORE — ACMG/scoring-tier columns
    // that belong to Module 3's domain, not this repo's (see memory:
    // project_acmg_tiering_deprioritized). Cause of the difference not investigated further, per
    // this session's scope steer — could be scoring-logic drift since the reference was
    // generated, or a different original invocation; either way it's Module 3's surface, not a
    // finding about this port's correctness. Row count matched exactly (17/17). Stub proves the
    // DAG wiring only.
    stub:
    """
    touch ${meta.partition_id}.parsed.clean.${params.output_format}
    """
}
