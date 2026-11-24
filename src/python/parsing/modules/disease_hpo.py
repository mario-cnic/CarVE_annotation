"""
Gene-priority tagging and filtering from a user gene file (--genes_file).
"""
import pandas as pd
from . import logger


def get_dynamic_priority(row, context: str):
    logger.debug(f"Evaluating gene priority for row with SYMBOL {row['SYMBOL']} and context {context}")
    if pd.isna(row.get('EVIDENCE')) or str(row.get('EVIDENCE')).strip() in ('', '.', 'nan', 'NaN', 'None'):
        return 'others'
    source = str(row['SOURCE'])
    evidence = str(row['EVIDENCE'])
    
    if context == 'Peripheral_and_conduction':
        if 'Peripheral_and_conduction' in source:
            return 'prioritary'
        return evidence.lower()

    elif context == 'Cardiomyopathies':
        if 'Cardiomyopathies' in source:
            return evidence.lower()
        if 'Peripheral_and_conduction' in source:
            return 'candidate'

    elif context == 'Aortic_and_vascular':
        if 'Aortic_and_vascular' in source:
            return evidence.lower()
        if 'Cardiomyopathies' in source:
            return 'candidate'
        return evidence.lower()
    
    elif context == 'all':
        if 'prioritary' in evidence.lower():
            return 'prioritary'
        return evidence.lower()

    logger.warning(f"Unrecognized context '{context}' for SYMBOL {row['SYMBOL']}; using evidence fallback")
    return evidence.lower()


def prioritize_genes(
    data: pd.DataFrame, gene_file: str, context: str = "all"
) -> pd.DataFrame:
    genes_df = pd.read_csv(gene_file)
    genes_df['GENE'] = genes_df['GENE'].str.strip()

    validate_gene_categories(
        context,
        genes_df['SOURCE'].str.split(";", expand=True).stack().unique()
    )
    
    data = data.merge(genes_df[['GENE', 'SOURCE', 'EVIDENCE']], 
                     left_on="SYMBOL", right_on="GENE", how="left", 
                     suffixes=('', '_genes'))
    
    suffix_cols = {'SOURCE_genes': 'SOURCE', 'EVIDENCE_genes': 'EVIDENCE'}
    for old_col, new_col in suffix_cols.items():
        if old_col in data.columns:
            data[new_col] = data[old_col]
            data = data.drop(columns=[old_col])
    
    data["gene_priority"] = data.apply(lambda row: get_dynamic_priority(row, context), axis=1)
    data["gene_priority"] = data["gene_priority"].fillna("others")
    
    gene_categories = ["prioritary", "secondary", "candidate", "others"]
    data["gene_priority"] = pd.Categorical(
        data["gene_priority"],
        categories=gene_categories,
        ordered=True
    )
    
    return data.drop(columns=["GENE"])


def validate_gene_categories(context, source_categories):
    if context != "all":
        categories = context.split(",")
        choices = source_categories.tolist()
        if any(x not in choices for x in categories):
            logger.error(
                f"Invalid gene category {categories}, should be one of {choices}"
            )
            raise ValueError(
                f"Invalid gene category {categories}, should be one of {choices}"
            )


def filter_by_gene_priority(data: pd.DataFrame, priorities: list[str]):
    if not priorities:
        logger.warning("No priorities provided for filtering, skipping")
        return data
    data = data[data["gene_priority"].isin(priorities)]
    logger.info(f"Dataframe filtered by gene priority {priorities}")
    return data
