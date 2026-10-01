"""
Clinical Variant Prioritization Engine.
Top-Level Pipeline Orchestrator and API Entrypoint.
Re-exports all modular sub-packages from `src.modules` for 100% backward compatibility.
"""
import argparse
import pandas as pd

try:
    from .modules.gene_identity import GeneIdentity
    from .modules.io_qc import (
        logger, LOGGING_MAP, GENOTYPE_QUAL_COLUMNS, GENE_PRIO_CHOICES, file_logging,
        _safe_numeric_series, _safe_string_series, parse_genotype, fix_dup_columns,
        build_cDNA_MANE, build_locus, build_qc_status, filter_variant_quality,
        filter_genotype_quality, select_columns, select_samples, merge_sample_cols,
        remove_columns, sort_by, order_columns, save_data,
    )
    from .modules.splicing import (
        SPLICING_COLUMNS, NEW_SPIP_COLUMNS, NEW_PANGOLIN_COLUMNS, NEW_BRANCHPOINT_COLUMNS,
        NEW_MAXENTSCAN_COLUMNS,
        parse_spip, decode_splicevault_prediction, parse_splicevault, parse_pangolin,
        parse_branchpointer, parse_maxentscan, build_signed_intron_offset, parse_spliceai_custom,
        build_spliceMAX, _extract_spliceai_details, _extract_spliceai_score,
    )
    from .modules.predictors import (
        DBNSFP_COLUMNS, MISSENSE_COLUMNS, parse_missense, parse_dbnsfp,
        parse_dbnsfp_by_row_transcript, load_aligned_columns,
    )
    from .modules.disease_hpo import (
        ACMG_SF_V3_2_CARDIAC_GENES, HPO_GENE_PANELS, load_cardiac_disease_gene_curation,
        apply_disease_phenotype_curation, apply_hpo_phenotype_weights, get_dynamic_priority,
        prioritize_genes, validate_gene_categories, filter_by_gene_priority,
    )
    from .modules.acmg import (
        _extract_clinvar_status, build_newImpact, build_acmg_criteria,
    )
    from .modules.pedigree import (
        read_pedigree_file, filter_phenotype, filter_by_dominant, mark_compound_het,
        _mark_pairwise_trans, _mark_unphased_candidates, _mark_with_pandas, _mark_with_dict,
        analyze_pedigree_inheritance,
    )
    from .modules.scoring import (
        build_priority_tier, order_clinical_columns, audit_annotation_availability,
        filter_freq, classification_dominant,
    )
except ImportError:
    from modules.gene_identity import GeneIdentity
    from modules.io_qc import (
        logger, LOGGING_MAP, GENOTYPE_QUAL_COLUMNS, GENE_PRIO_CHOICES, file_logging,
        _safe_numeric_series, _safe_string_series, parse_genotype, fix_dup_columns,
        build_cDNA_MANE, build_locus, build_qc_status, filter_variant_quality,
        filter_genotype_quality, select_columns, select_samples, merge_sample_cols,
        remove_columns, sort_by, order_columns, save_data,
    )
    from modules.splicing import (
        SPLICING_COLUMNS, NEW_SPIP_COLUMNS, NEW_PANGOLIN_COLUMNS, NEW_BRANCHPOINT_COLUMNS,
        NEW_MAXENTSCAN_COLUMNS,
        parse_spip, decode_splicevault_prediction, parse_splicevault, parse_pangolin,
        parse_branchpointer, parse_maxentscan, build_signed_intron_offset, parse_spliceai_custom,
        build_spliceMAX, _extract_spliceai_details, _extract_spliceai_score,
    )
    from modules.predictors import (
        DBNSFP_COLUMNS, MISSENSE_COLUMNS, parse_missense, parse_dbnsfp,
        parse_dbnsfp_by_row_transcript, load_aligned_columns,
    )
    from modules.disease_hpo import (
        ACMG_SF_V3_2_CARDIAC_GENES, HPO_GENE_PANELS, load_cardiac_disease_gene_curation,
        apply_disease_phenotype_curation, apply_hpo_phenotype_weights, get_dynamic_priority,
        prioritize_genes, validate_gene_categories, filter_by_gene_priority,
    )
    from modules.acmg import (
        _extract_clinvar_status, build_newImpact, build_acmg_criteria,
    )
    from modules.pedigree import (
        read_pedigree_file, filter_phenotype, filter_by_dominant, mark_compound_het,
        _mark_pairwise_trans, _mark_unphased_candidates, _mark_with_pandas, _mark_with_dict,
        analyze_pedigree_inheritance,
    )
    from modules.scoring import (
        build_priority_tier, order_clinical_columns, audit_annotation_availability,
        filter_freq, classification_dominant,
    )


def argparsing() -> argparse.ArgumentParser:
    """Parse args either from cli or passed as arguments."""
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", type=str, help="path to intermediate tsv file")
    parser.add_argument("--output", type=str, help="file path and name")
    parser.add_argument(
        "--columns",
        type=str,
        help="columns to select or path to csv file with columns to select",
        default=None
    )
    parser.add_argument(
        "--genes_file", type=str, help="path to csv file with genes of interest"
    )
    parser.add_argument(
        "--gene_priority",
        type=str,
        help="Filter genes by priority [Prioritary, Secondary, Candidate] or a comma separated list or all",
        default=None,
    )
    parser.add_argument(
        "--genes_category",
        type=str,
        help="[DEPRECATED] Filter genes by priority",
        default=None,
    )
    parser.add_argument(
        "--disease_context",
        type=str,
        help="Disease context to filter genes, if --genes_category is not 'all' and --genes_file is provided",
        default='all',
    )
    parser.add_argument(
        "--genotype_file",
        type=str,
        help="Path to PED/JSON file with genotypes to filter or a dictionary-like string.",
        default=None,
    )
    parser.add_argument(
        "--genotype_filter",
        action="store_false",
        help="Filter by genotype instead of just labeling variants with matching genotypes",
    )
    parser.add_argument(
        "--genotype_group",
        type=str,
        help="Group of samples to filter by genotype.",
        default=None,
    )
    parser.add_argument(
        "--freq",
        type=float,
        help="frequency to filter variants",
        default=None,
    )
    parser.add_argument(
        "--freq_col",
        type=str,
        help="column to filter variants",
        default="gnomADv4_AF_grpmax_joint",
    )
    parser.add_argument(
        "--classification",
        type=str,
        help="Add classification to filter variants",
        default="",
    )
    parser.add_argument(
        "--merge_samples",
        action="store_true",
        help="merge samples columns and remove original columns",
    )
    parser.add_argument(
        "--logging_level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Log level to show [DEBUG,INFO,WARNING,ERROR]",
        default="INFO",
    )
    parser.add_argument(
        "--hgnc-table",
        type=str,
        default=None,
        help="HGNC complete set TSV (symbol history -> Ensembl gene ID); required with --spliceai-symbol-map",
    )
    parser.add_argument(
        "--spliceai-symbol-map",
        type=str,
        default=None,
        help="resources/spliceai_symbol_to_ensg.*.tsv (SpliceAI symbol -> Ensembl gene ID by coordinates)",
    )
    parser.add_argument(
        "--dbnsfp-aligned-columns",
        type=str,
        default=None,
        help="resources/dbnsfp_transcript_aligned_columns.txt (dbNSFP columns whose '&' lists align with "
             "Ensembl_transcriptid); required when dbNSFP columns are present",
    )
    parser.add_argument(
        "--enst-spip-map",
        type=str,
        default=None,
        help="resources/enst_to_spip_nm.*.tsv (Ensembl transcript -> SPiP RefSeq transcript by exon structure)",
    )
    parser.add_argument(
        "--skip_quality_filter",
        action="store_true",
        help="Skip quality filter, this may lead to false positives in the results",
    )
    parser.add_argument(
        "--log_file",
        type=str,
        help="path to log file",
        default=None,
    )
    parser.add_argument(
        "--add_format",
        action="store_true",
        help="Add FORMAT column and it's samples the output file (needs 'GT_' columns)",
    )
    parser.add_argument(
        "--low_gq_to_missing",
        action="store_true",
        help="Convert low GQ (<10) values to MISSING instead of removing the row",
    )
    parser.add_argument(
        "--dbnsfp_cols",
        type=str,
        help="path to file with dbNSFP columns to parse and add to the dataframe",
        default="data/dbNSFP_parsing.tsv",
    )

    return parser


def main(args: list[str] | None):
    args = argparsing().parse_args(args)

    file_logging(path="_log", filename=args.log_file if args.log_file else args.output)
    logger.setLevel(LOGGING_MAP[args.logging_level])

    logger.info("Starting script to filter variants and add new columns...")

    if args.genes_category:
        logger.error("--genes_category is deprecated, please use --disease_context and/or --gene_priority")
        raise ValueError("--genes_category is deprecated, please use --disease_context and/or --gene_priority")

    if args.gene_priority:
        priorities = args.gene_priority.split(",")
        if any(x not in GENE_PRIO_CHOICES for x in priorities):
            logger.error(f"Invalid gene priority {priorities}, should be one of {GENE_PRIO_CHOICES}")
            raise ValueError(f"Invalid gene priority {priorities}, should be one of {GENE_PRIO_CHOICES}")

    first_cols = [
        "Locus",
        "SYMBOL",
        "QC_STATUS",
        "QC_DETAIL",
        "HGVSc",
        "HGVSp",
        "CDNA_NAME",
        "PRIORITY_TIER",
        "VARIANT_PRIORITY_SCORE",
        "NEW_IMPACT",
        "IMPACT",
        "intron_offset_signed",
        "splice_side",
        f"{args.freq_col}",
        "spliceAI_MAX",
        "SpliceAI_status",
        "spliceai_custom_MAX",
        "spliceai_custom_DS_AG",
        "spliceai_custom_DS_AL",
        "spliceai_custom_DS_DG",
        "spliceai_custom_DS_DL",
        "spliceai_custom_DP_AG",
        "spliceai_custom_DP_AL",
        "spliceai_custom_DP_DG",
        "spliceai_custom_DP_DL",
        "spliceai_custom_SYMBOL",
        "SPiP_interpretation",
        "SPiP_prediction",
        "SPiP_score",
        "SPiP_mechanism",
        "SPiP_status",
        "SpliceVault_top_events",
        "SpliceVault_Predictions_Decoded",
        "SpliceVault_site_sample_count",
        "SpliceVault_site_max_depth",
        "SpliceVault_out_of_frame_events",
        "SpliceVault_site_pos",
        "SpliceVault_site_type",
        "SpliceVault_SpliceAI_delta",
        "SpliceVault_status",
        "Pangolin_max_score",
        "Pangolin_heart_lv_score",
        "Pangolin_heart_aa_score",
        "Pangolin_status",
        "Branchpointer_prob",
        "Branchpointer_U2_energy",
        "Branchpoint_disrupted",
        "Branchpoint_status",
        "LaBranchoR_score",
        "LaBranchoR_acc_dist",
        "LaBranchoR_status",
        "MaxEntScan_ref",
        "MaxEntScan_alt",
        "MaxEntScan_diff",
        "MaxEntScan_pct_decrease",
        "MaxEntScan_disrupted",
        "MaxEntScan_status",
        "5UTR_annotation",
        "5UTR_consequence",
        "gene_priority",
        "SOURCE",
        "EVIDENCE",
        "pLI_gene_value"
    ]

    if args.input.endswith(".parquet") or args.input.endswith(".pq"):
        logger.info(f"Loading input file {args.input} as parquet")
        try:
            data = pd.read_parquet(args.input)
        except ImportError as e:
            logger.error("Failed to read Parquet file. Please ensure 'pyarrow' or 'fastparquet' is installed.")
            raise e
    elif args.input.endswith(".feather") or args.input.endswith(".arrow"):
        logger.info(f"Loading input file {args.input} as feather")
        try:
            data = pd.read_feather(args.input)
        except ImportError as e:
            logger.error("Failed to read Feather file. Please ensure 'pyarrow' is installed.")
            raise e
    else:
        compression = "gzip" if args.input.endswith(".gz") else "infer"
        data = pd.read_csv(
            args.input,
            na_values=".",
            sep="\t",
            compression=compression,
            header=0,
            index_col=False,
            low_memory=False,
        )
    logger.debug(f"Dataframe columns: {data.columns}")
    freq_cols = [x for x in data.columns if args.freq_col in x]
    logger.debug(f"Frequency columns: {freq_cols}")
    logger.info(f"Input file {args.input} loaded with shape {data.shape}")

    data = select_columns(data, column_file=args.columns, formato=args.add_format)
    data = build_cDNA_MANE(data)
    data = build_signed_intron_offset(data)
    data = build_locus(data)
    data = fix_dup_columns(data)
    data = parse_genotype(data)

    if args.genotype_file:
        pedigree = read_pedigree_file(args.genotype_file)
        data = filter_phenotype(data, pedigree)
        first_cols += ["compound_het"]
    else:
        pedigree = None
        logger.warning("No genotype file provided, skipping filtering by genotype")
    
    if args.freq:
        data = filter_freq(data, freq=args.freq, freq_column=args.freq_col)
    else:
        logger.warning("No frequency provided, skipping filtering by frequency")

    if args.skip_quality_filter:
        logger.warning("Skipping quality filter")
    else:
        logger.info("Filtering data by variant quality score")
        data = filter_variant_quality(data, qual_column="QUAL", variant_quality_score=20)
        logger.info("Filtering data by genotype quality score")
        data = filter_genotype_quality(data, columns=GENOTYPE_QUAL_COLUMNS, genotype_quality_score=30, append_lowqual=args.low_gq_to_missing)

    if args.genes_file:
        data = prioritize_genes(data, args.genes_file, args.disease_context)
        if args.gene_priority:
            data = filter_by_gene_priority(data, args.gene_priority)

    needs_identity = any(c in data.columns for c in ("SpliceAI", "SPiP"))
    identity = None
    if needs_identity:
        if not (args.hgnc_table and args.spliceai_symbol_map and args.enst_spip_map):
            raise ValueError(
                "SpliceAI/SPiP columns are present: --hgnc-table, --spliceai-symbol-map and --enst-spip-map are required so "
                "scores are attached to the row's own gene/transcript. Refusing to fall back to unmatched scores."
            )
        identity = GeneIdentity(args.spliceai_symbol_map, args.hgnc_table, args.enst_spip_map)

    data = parse_spliceai_custom(data, identity=identity)
    data = build_spliceMAX(data, SPLICING_COLUMNS)

    data = parse_spip(data, identity=identity)
    if "Ensembl_transcriptid" in data.columns:
        if not args.dbnsfp_aligned_columns:
            raise ValueError(
                "dbNSFP columns are present: --dbnsfp-aligned-columns is required so each row keeps the value "
                "of its own transcript. Refusing to leave per-transcript lists on every transcript row."
            )
        data = parse_dbnsfp_by_row_transcript(data, load_aligned_columns(args.dbnsfp_aligned_columns))
    data = parse_dbnsfp(data, transcript=None, dbnsfp_cols=args.dbnsfp_cols)
    data = parse_missense(data, MISSENSE_COLUMNS)
    data = parse_splicevault(data)
    data = parse_pangolin(data)
    data = parse_branchpointer(data)
    data = parse_maxentscan(data)

    data = build_qc_status(data)
    data = build_newImpact(data)
    data = analyze_pedigree_inheritance(data, pedigree)
    data = build_priority_tier(data, freq_col=args.freq_col)

    if args.classification == 'dominant':
        data = classification_dominant(data=data, freq_col=args.freq_col, freq=args.freq)
        first_cols += ["Variant_rank"]
    elif args.classification:
        logger.warning(f"Invalid classification {args.classification}, should be one of [dominant]")
    
    if args.merge_samples:
        data = merge_sample_cols(data, pattern="GT_")
        data = merge_sample_cols(data, pattern="GQ_")
        sample_cols = ["merged_GT", "count_GT", "merged_GQ"]
        if any(col.startswith("AD_") for col in data.columns):
            data = merge_sample_cols(data, pattern="AD_")
            sample_cols += ["merged_AD"]
        if any(col.startswith("DP_") for col in data.columns):
            data = merge_sample_cols(data, pattern="DP_")
            sample_cols += ["merged_DP"]
        logger.info("Merging sample columns")
        first_cols += sample_cols

    data = remove_columns(data, columns=["CHROM", "POS", "REF", "ALT"])

    if pedigree is not None and isinstance(pedigree, (dict, pd.DataFrame)):
        data = mark_compound_het(data, pedigree)
    else:
        logger.warning("Pedigree is not properly initialized or valid, skipping marking compound heterozygous variants")

    data = order_columns(data, columns=first_cols)

    sorting_dict = {"NEW_IMPACT": True, "SPiP_prediction": False, args.freq_col: True}
    if args.classification:
        sorting_dict["Variant_rank"] = True

    data = sort_by(data, by=sorting_dict)
    save_data(data, args.output, args.log_file)


if __name__ == "__main__":
    main(None)
