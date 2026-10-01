###################################
# vcf_parser                      #
# read,write VCFs and transforms  #
#                                 #
###################################

# import pandas as pd
from logger import logger, LOGGING_MAP, file_logging
import argparse
from vcf_parser import VCFParser  # requires vcf_parser package from pip
from typing import TextIO
import gzip
import os
import re


def get_vep_header(vcf, arg_vep_columns) -> list[str]:
    """Retrieves VEP columns specified by the user.

    Params:
        - vcf: vcf file
        - selected_vep_columns:CSQ(VEP) columns specified by the user

    Returns:
        VEP columns as list of strings
    """
    final_vep_columns = []
    if not arg_vep_columns:
        logger.info("No VEP columns specified, using all VEP columns")
        arg_vep_columns = []
    elif arg_vep_columns.endswith(".txt"):
        with open(arg_vep_columns, "r") as f:
            lines = f.readlines()
        selected_vep_columns = [x.strip() for x in lines]
        selected_vep_columns = [x.replace(",", "") for x in selected_vep_columns]
        logger.info(f"Selected VEP columns from file {arg_vep_columns}")
    else:
        selected_vep_columns = arg_vep_columns
        logger.info(f"Selected VEP columns {selected_vep_columns}")
        selected_vep_columns = selected_vep_columns.split(",")

    if selected_vep_columns and vcf.metadata.vep_columns:
        # For each user-specified column, include all columns in vep_columns that match base or suffixed
        for user_col in selected_vep_columns:
            # Match exact or suffixed with .[digits] or without leading underscore
            clean_col = user_col.lstrip("_")
            pattern = re.compile(rf"^(_)?({re.escape(clean_col)}|{re.escape(user_col)})(\.\d+)?$")
            matched = [col for col in vcf.metadata.vep_columns if pattern.match(col)]
            if not matched:
                logger.warning(f"{user_col} not found in vep_columns")
            else:
                final_vep_columns.extend(matched)
        logger.info(f"VEP columns selected and present {', '.join(final_vep_columns)}")
        logger.info("VEP columns added to header")
    elif vcf.metadata.vep_columns:
        logger.info("No VEP columns specified, using all VEP columns")
        final_vep_columns = vcf.metadata.vep_columns
    else:
        logger.error("No VEP columns found in VCF")
        raise ValueError("No VEP columns found in VCF, please remove `--add_vep` flag")
    return final_vep_columns


def get_samples(all_individuals: list, selected_samples: str) -> list[str]:
    """Check and return samples in vcf specified by the user

    Params:
        - all_individuals: list of samples in VCF
        - selected_samples: list of samples specified by the user

    Returns:
        - List of samples
    """

    # if selected_samples is a path, retrieve samples from file function
    def get_samples_from_file(path: str) -> list[str]:
        """Retrieve samples from file"""
        with open(path, "r") as f:
            lines = f.readlines()
        samples = [x.strip() for x in lines]
        return samples

    # check if selected_samples is a path
    if selected_samples.endswith(".txt"):
        logger.info(f"Selected samples from file {selected_samples}")
        final_sample_columns = get_samples_from_file(selected_samples)
    elif selected_samples == "all":
        logger.info("Including all samples")
        final_sample_columns = all_individuals
    else:
        selected_samples = selected_samples.split(",")
        final_sample_columns = [x for x in all_individuals if x in selected_samples]
        if not final_sample_columns:
            logger.warning(f"{selected_samples} not found in VCF, looking for exact matches failed, trying partial matches")
            final_sample_columns = [x for x in all_individuals if any(s in x for s in selected_samples)]
            if not final_sample_columns:
                logger.error(f"{selected_samples} not found in VCF, using all samples")
                final_sample_columns = all_individuals
            else:
                logger.warning(f"Selected samples with partial matches {', '.join(final_sample_columns)}")

    logger.info(f"Sample columns selected {', '.join(final_sample_columns)}")

    return final_sample_columns


def add_sample_header(
    user_samples: list[str], vcf_samples: list[str], header: list[str]
) -> list[str]:
    """Adds sample names to header."""
    # for sample in user_samples:
    #     if sample in vcf_samples:
    #         header.append(sample)
    # return header
    pass


# def parse_sample_calling() -> None:
#     """Create new column for each sample"""
def read_gene_set(gene_file: str) -> list[str]:
    """Reads gene set from file and returns list of genes"""
    if gene_file.endswith(".txt"):
        with open(gene_file, "r") as f:
            lines = f.readlines()
    elif gene_file.endswith(".csv"):
        with open(gene_file, "r") as f:
            # skip first line, header
            lines = f.readlines()[1:]
            # keep only first column
            lines = [x.split(",")[0] for x in lines]
    else:
        lines = gene_file.split(",")
        if len(lines) == 1:
            lines = [gene_file]
        logger.info(f"Gene set specified as comma separated list {', '.join(lines)}")
    genes = [x.strip().upper() for x in lines]
    return genes


def write_header(
    vcf: VCFParser,
    file: TextIO,
    include_info: bool,
    include_GT: bool,
    include_GQ: bool,
    include_DP: bool,
    include_AD: bool,
    include_VEP: bool,
    vep_columns: str,
    selected_samples: str,
) -> list:
    """
    Finds VCF header and writes to CSV

    Params:
        - remove_info: Removes non formatted INFO field
        - include_GT: Creates new column with GT information
        only extracted from sample fields (0/1,1/1...)
        - include_VEP: Creates VEP headers
    Returns:
        - columns to use in parsing from default VCF
    """
    # # get base csv header 
    base_header = vcf.header.copy() #type: list
    gt_columns = []
    expanded_header = base_header.copy()
    vcf_samples = vcf.metadata.individuals
    logger.debug(f"VCF samples: {', '.join(vcf_samples)}")

    if selected_samples:
        samples = get_samples(vcf_samples, selected_samples)
        logger.debug(f"{','.join(samples)}")
        logger.debug(f"Individuals selected: {','.join(samples)}")
        not_samples = [x for x in vcf_samples if x not in samples]
        logger.debug(f"Individuals not selected: {','.join(not_samples)}")
        for x in not_samples:
            expanded_header.remove(x)
    else:
        logger.info("No samples were specified, including all samples")
        samples = vcf_samples

    # 1. Unconditionally remove "INFO" (with a safety check to prevent ValueErrors)
    if "INFO" in expanded_header:
        expanded_header.remove("INFO")
        base_header.remove("INFO")
        logger.info("INFO field removed from header and will be replaced with individual INFO columns")


    # 2. Handle VEP logic first so we only compute it once
    if include_VEP:
        vep_columns = get_vep_header(vcf, vep_columns)
    else:
        vep_columns = []

    # 3. Handle INFO logic
    info_columns = []
    if not include_info:
        logger.info("INFO field removed")
    else:
        # Iterate directly over the dict keys, no need for list() or .keys()
        info_columns = [k for k in vcf.metadata.info_dict if k != "CSQ"]
        
        # Deduplicate if VEP is also included
        if include_VEP and vep_columns:
            logger.info("Duplicate columns between INFO and VEP will be included only as VEP columns")
            vep_cols_set = set(vep_columns) # Convert to set for much faster lookups
            dup_info_cols = [col for col in info_columns if col in vep_cols_set]
            logger.warning(f"Duplicate columns between INFO and VEP: {', '.join(dup_info_cols)} removed from INFO columns")
            info_columns = [col for col in info_columns if col not in vep_cols_set]
        # Add info columns to the header
        expanded_header.extend(info_columns)

    # 4. Finally, append VEP columns to the header
    if include_VEP:
        expanded_header.extend(vep_columns)

    # Dynamically add genotype-related columns based on selected flags
    genotype_fields = []
    if include_GT:
        genotype_fields.append("GT")
    if include_GQ:
        genotype_fields.append("GQ")
    if include_DP:
        genotype_fields.append("DP")
    if include_AD:
        genotype_fields.append("AD")

    if genotype_fields:
        logger.info(
            f"Including genotype columns for each sample: {', '.join(genotype_fields)}"
        )
        genotype_field_columns = [
            f"{field}_{sample}" for sample in samples for field in genotype_fields
        ]
        expanded_header = expanded_header + genotype_field_columns

    else:
        logger.info(
            "No genotype or genotype quality columns added, variants without genotypes for your samples will not be excluded"
        )

    logger.debug(f"Written header: {expanded_header}")
    file.write("\t".join(expanded_header) + "\n")  # write header
    # for value extraction we need to know which columns are base (from VCF) and which are added (INFO, VEP, genotype fields)
    expanded_wo_base = [x for x in expanded_header if x not in base_header]
    return base_header, info_columns, vep_columns, expanded_wo_base, samples


def write_values(
    vcf: VCFParser,
    file: TextIO,
    base_header: list[str],
    vep_cols: list[str],
    info_cols: list[str],
    include_gt: bool,
    include_gq: bool,
    include_dp: bool,
    include_ad: bool,
    selected_individuals: list[str],
    filter_by: str,
    selected_genes: list[str] = None,
):
    """
    Writes variant data to the output file.

    Params:
        - vcf: VCFParser object.
        - file: Output file object.
        - header: Base header columns.
        - vep_cols: VEP columns to include.
        - info_cols: INFO columns to include.
        - include_gt: Whether to include genotype columns.
        - include_gq: Whether to include genotype quality columns.
        - selected_individuals: List of individuals to include.
        - selected_genes: List of genes to filter by.
        - filter_by: The VEP column to filter by (default: "SYMBOL").
    """
    progress = 0
    skipped_variants = 0
    if selected_genes and not vep_cols:
        logger.error(
            f"Gene selection requires the {filter_by} VEP column. Please specify --add_vep. And make sure the VEP field is present in the VCF."
        )
        raise ValueError(
            f"Gene selection requires the {filter_by} VEP column. Please specify --add_vep. And make sure the VEP field is present in the VCF."
        )

    for progress, variant in enumerate(vcf, start=progress):
        logger.debug(f"Processing variant #{progress}: {variant['variant_id']}")
        try:
            base_variant = _get_base_values(variant, base_header)
            logger.debug(f"Base variant: {base_variant}")

            if vep_cols:
                skipped_variants += _process_vep_columns(
                    variant,
                    file,
                    base_variant,
                    vep_cols,
                    info_cols,
                    include_gt,
                    include_gq,
                    include_dp,
                    include_ad,
                    selected_individuals,
                    selected_genes,
                    filter_by
                )
            elif info_cols:
                info_values = _get_info_values(variant, info_cols)
                row_values = base_variant + info_values
                if include_gt or include_gq or include_dp or include_ad:
                    skipped_variants += _process_and_write_genotype_data(
                        variant,
                        file,
                        row_values,
                        include_gt,
                        include_gq,
                        include_dp,
                        include_ad,
                        selected_individuals,
                    )
                else:
                    file.write("\t".join(str(x) for x in row_values) + "\n")
            elif include_gt or include_gq or include_dp or include_ad:
                skipped_variants += _process_and_write_genotype_data(
                    variant,
                    file,
                    base_variant,
                    include_gt,
                    include_gq,
                    include_dp,
                    include_ad,
                    selected_individuals,
                )
            else:
                file.write("\t".join(str(x) for x in base_variant) + "\n")

            if progress % 1000000 == 0 and progress > 0:
                logger.info(f"{progress} variants processed...")

        except Exception as e:
            logger.error(f"Unexpected error at variant #{progress}: {e}")
            skipped_variants += 1
            continue
    else:
        logger.info(f"Processed {progress} variants.")

    logger.info(f"Variants skipped = {skipped_variants}")


def _get_base_values(variant: dict, header: list[str]) -> list[str]:
    """Extracts base variant VALUES as a list of values (not a string).
    Param:
    - variant: dictionary containing variant information
    - header: list of base header columns to extract (e.g., CHROM, POS, ID, REF, ALT, QUAL, FILTER)
    Returns:
    - List of values corresponding to the base header columns, in the same order as the header
    """
    logger.debug(f"[_get_base_values] The base header is: {header}")
    logger.debug(f"Extracting base values for variant {variant['variant_id']} using header {header}")
    values = [variant.get(h, ".") for h in header]
    logger.debug(f"Extracted base values: {dict(zip(header, values))}")
    return values


def _process_vep_columns(
    variant: dict,
    file: TextIO,
    base_variant: list[str],
    vep_cols: list[str],
    info_cols: list[str],
    include_gt: bool,
    include_gq: bool,
    include_dp: bool,
    include_ad: bool,
    selected_individuals: list[str],
    selected_genes: list[str],
    filter_by: str = "SYMBOL",
) -> int:
    """Processes VEP columns and writes to file."""
    skipped_variants = 0
    for transcript in variant.get("vep_info", {}).values():
        for entry in transcript:
            # Filter by `filter_by` column (gene) if values are provided (e.g., SYMBOL and ACTC1 or Feature and ENST00000361445)
            if selected_genes and (entry.get(filter_by) not in selected_genes):
                logger.debug(
                    f"[{variant['variant_id']}] Skipping variant because {filter_by} {entry.get(filter_by)} is not in the selected {filter_by.lower()} set"
                )
                skipped_variants += 1
                continue

            vep_values = [entry.get(col, ".") for col in vep_cols]
            logger.debug(f"[{variant['variant_id']}] VEP values for {dict((col, entry.get(col, '.')) for col in vep_cols)}")
            gnomAD_cols = [col for col in vep_cols if "gnomADv4" in col]
            logger.debug(f"[{variant['variant_id']}] gnomAD values: {dict((col, entry.get(col, '.')) for col in gnomAD_cols)}")
            if "gnomADv4_AF_grpmax_joint" in entry:
                # logger.debug(
                #     f"[{variant['variant_id']}] gnomADv4_AF_grpmax_joint value: {entry['gnomADv4_AF_grpmax_joint']}"
                # )
                pass
            # logger.debug(f"[{variant['variant_id']}] VEP values: {vep_values} that will be added")
            info_values = _get_info_values(variant, info_cols)
            # base_variant is now a list, so concatenate all values in order
            row_values = base_variant + info_values + vep_values

            if include_gt or include_gq or include_dp or include_ad:
                # Pass as a list, not a string
                skipped_variants += _process_and_write_genotype_data(
                    variant,
                    file,
                    row_values,
                    include_gt,
                    include_gq,
                    include_ad,
                    include_dp,
                    selected_individuals,
                )
            else:
                file.write("\t".join(str(x) for x in row_values) + "\n")
    return skipped_variants


def _get_info_values(variant: dict, info_cols: list[str]) -> list[str]:
    """Extracts INFO column values."""
    values = []
    logger.debug(f" Variant dictionary INFO: {variant.get('info_dict', {})}")
    for col in info_cols:
        value = variant.get("info_dict", {}).get(col, ".")
        if isinstance(value, list):
            value = "&".join(value)
        # logger.debug(f"INFO column {col} value: {value}")
        val_str = str(value) if value is not None else "."
        val_str = val_str.replace(r"\t", "&")
        if val_str.lower() == "nan" or val_str == ".":
            val_str = "."
        values.append(val_str)
    return values

def _process_and_write_genotype_data(
    variant: dict,
    file: TextIO,
    variant_line: list,
    include_gt: bool,
    include_gq: bool,
    include_dp: bool,
    include_ad: bool,
    selected_individuals: list[str],
) -> int:
    """Processes and writes genotype and genotype quality data."""
    skip_variant = True
    # variant_line is a list of values up to this point
    # For each selected individual, append genotype fields in the correct order
    gt_values = []
    for individual in selected_individuals:
        genotype = variant.get("genotypes", {}).get(individual)
        if genotype:
            if include_gt:
                gt_values.append(genotype.genotype)
            if include_gq:
                gt_values.append(str(genotype.genotype_quality))
            if include_dp:
                gt_values.append(str(genotype.depth_of_coverage))
            if include_ad:
                gt_values.append(str(genotype.allele_depth))
            if include_gt and "1" in genotype.genotype:
                skip_variant = False
        else:
            # If no genotype, fill with dots for all requested fields
            if include_gt:
                gt_values.append(".")
            if include_gq:
                gt_values.append(".")
            if include_dp:
                gt_values.append(".")
            if include_ad:
                gt_values.append(".")
    # If no GT requested, don't filter by called variants
    if not include_gt:
        skip_variant = False
    if skip_variant:
        logger.debug(
            f"[NOT CALLED FOR SAMPLES] - Skipping variant {variant['variant_id']}"
        )
    else:
        # Write the full row: base + info + vep + genotype columns
        file.write("\t".join(str(x) for x in variant_line + gt_values) + "\n")
    return 1 if skip_variant else 0


def file_handling(args):
    if os.path.exists(args.output) and not args.overwrite:
        logger.error(
            f"File {args.output} already exists. Please remove it or use a different name."
        )
        raise FileExistsError(
            f"File {args.output} already exists. Please remove it or use a different name."
        )
    # open file
    if args.output.endswith(".gz"):
        f = gzip.open(args.output, mode="wt", compresslevel=9)
    else:
        f = open(args.output, "w")
    logger.info(f"Creating file {f.name}")
    return f


def argparsing() -> argparse.ArgumentParser:
    """
    Parse args either from cli or passed as arguments"""

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", type=str, help="path to data")
    parser.add_argument("--output", type=str, help="file path and name")
    parser.add_argument(
        "--vep_columns",
        type=str,
        help="VEP columns to select",
    )
    parser.add_argument(
        "--sample_columns",
        type=str,
        help="path to txt containing sample names to select or\
comma separated list of sample names to select. Use 'all' to select all samples.\
WARNING: variants not called or homozygous reference for all individuals will be excluded",
    )
    parser.add_argument(
        "--add_info",
        action="store_true",
        help="Include all INFO field unparsed",
    )
    parser.add_argument(
        "--add_gt",
        action="store_true",
        help="Include one genotype column per sample (only GT)",
    )
    parser.add_argument(
        "--add_vep",
        action="store_true",
        help="Include parsed CSQ fields (specify which fields with --vep_columns param)",
    )
    parser.add_argument(
        "--logging_level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Log level to show [DEBUG,INFO,WARNING,ERROR]",
        default="INFO",
    )
    parser.add_argument(
        "--gene_set",
        type=str,
        help="Path to file containing gene names to select",
    )

    parser.add_argument(
        "--add_gq",
        action="store_true",
        help="Include one genotype quality column per sample (only GT)",
    )

    parser.add_argument(
        "--add_dp",
        action="store_true",
        help="Include depth column per sample (only GT)",
    )

    parser.add_argument(
        "--add_ad",
        action="store_true",
        help="Include allele depth column per sample (only GT)",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing file",
    )
    parser.add_argument(
        "--filter_by",
        type=str,
        help="VEP column to filter by when using --add_vep and gene selection (default: SYMBOL)",
        default="SYMBOL",
    )

    return parser


def main(args: list[str] | None):
    args = argparsing().parse_args(args)

    file_logging(path="_log", filename=args.output, save_in_log=True)
    logger.setLevel(LOGGING_MAP[args.logging_level])

    my_parser = VCFParser(infile=args.input, split_variants=True)
    # check if output file exists and raise error if it does
    f = file_handling(args)

    # read gene set file
    if args.gene_set:
        genes = read_gene_set(args.gene_set)
        logger.info(f"Gene-Set selected:\n{', '.join(genes)}")
    else:
        genes = None
        logger.info("No Gene-Set selected. All genes will be included.")

    # write header
    base_header, info_cols, vep_cols, expanded_header, individuals = write_header(
        vcf=my_parser,
        file=f,
        include_info=args.add_info,
        include_GT=args.add_gt,
        include_GQ=args.add_gq,
        include_AD=args.add_ad,
        include_DP=args.add_dp,
        include_VEP=args.add_vep,
        vep_columns=args.vep_columns,
        selected_samples=args.sample_columns,
    )
    logger.debug(f"The column ORDER IS: {base_header}+{info_cols}")
    logger.debug(f"[main] The base header is: {base_header}")
    # write values
    write_values(
        vcf=my_parser,
        file=f,
        base_header=base_header,
        vep_cols=vep_cols,
        info_cols=info_cols,
        include_gt=args.add_gt,
        include_gq=args.add_gq,
        include_dp=args.add_dp,
        include_ad=args.add_ad,
        selected_individuals=individuals,
        selected_genes=genes,
        filter_by=args.filter_by
    )
    # close file
    logger.info(f"Closing file {f.name}")
    f.close()


if __name__ == "__main__":
    main(None)
