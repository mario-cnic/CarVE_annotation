import os
import logging
import re
import numpy as np
import pandas as pd
from typing import Literal

logger = logging.getLogger("clinical_prioritization")

LOGGING_MAP = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
}

GENE_PRIO_CHOICES = ["prioritary", "secondary", "candidate", "others", "all"]
GENOTYPE_QUAL_COLUMNS = "^GQ_*"


def file_logging(path="_log", filename=None):
    if not os.path.exists(path):
        os.makedirs(path)
    log_name = filename if filename else "filtering.log"
    log_path = os.path.join(path, os.path.basename(log_name) + ".log")
    fh = logging.FileHandler(log_path, mode="w")
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)



def _safe_numeric_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    return pd.Series(0.0, index=df.index)


def _safe_string_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col in df.columns:
        return df[col].astype(str)
    return pd.Series("", index=df.index)


def fix_dup_columns(data: pd.DataFrame):
    """
    If there is two columns, one with .1 and the other without merge them
    into a single one, getting the column with non empty value. do not combine them
    """
    logger.info("Fixing duplicated columns")
    cols = data.columns
    dup_cols = []
    for col in cols:
        if col.endswith(".1") or col.endswith(".2"):
            logger.debug(f"Found duplicated column {col}")
            data[col] = data[col].replace("", np.nan)
            data[col] = data[col].replace(".", np.nan)
            dup_cols.append(col)
            col_no_dot = col[:-2]
            if col_no_dot in cols:
                try:
                    data[col_no_dot] = data[col_no_dot].fillna(data[col])
                    logger.debug(f"Merging columns {col} into {col_no_dot}")
                except ValueError:
                    logger.warning(
                        f"Error merging columns {col} and {col_no_dot}, skipping"
                    )
                    continue
                data.drop(columns=[col], inplace=True)
                logger.debug(f"Columns {col} and {col_no_dot} merged into {col_no_dot}")
    logger.info(f"Duplicated columns {', '.join(dup_cols)} fixed")
    cols = data.columns.to_list()
    dup_cols = [col for col in cols if cols.count(col) > 1]
    if dup_cols:
        logger.warning(
            f"Duplicated columns {', '.join(dup_cols)} found, please check the data"
        )
        data = data.loc[:, ~data.columns.duplicated()]
    return data


def parse_genotype(data: pd.DataFrame):
    """Transform all GT_ columns to HET, HOMREF, HOMALT
    HET: 0/1, 1/0, 0|1, 1|0
    HOMREF: 0/0, 0|0
    HOMALT: 1/1, 1|1
    MISSING: ./., .|. or NaN
    """
    GENOTYPE_MAPPING = {
        r"^(\./\.|\.)$": "MISSING",
        r"^(0/1|1/0|0\|1|1\|0)$": "HET",
        r"^(0/0|0\|0)$": "HOMREF",
        r"^(1/1|1\|1)$": "HOMALT",
    }
    logger.info("Fixing genotype columns")
    gt_cols = [col for col in data.columns if col.startswith("GT_")]
    if not gt_cols:
        logger.warning("No GT_ columns found in data, skipping")
        return data
    logger.debug(f"First 5 rows of GT_ columns:\n{data[gt_cols].head()}")
    logger.debug(f"Unique values for GT_ columns: {data[gt_cols].nunique()}")
    data[gt_cols] = data[gt_cols].replace(
        GENOTYPE_MAPPING,
        regex=True,
    )
    logger.info(f"Genotype columns {''.join(gt_cols)} fixed")
    logger.debug(f"First 5 rows of GT_ columns:\n{data[gt_cols].head()}")
    return data


def build_locus(data):
    """
    Create new column Locus with CHROM,POS,REF,ALT
    """
    if "Locus" in data.columns:
        logger.info("Column `Locus` already exists, skipping creation")
        return data

    data["Locus"] = (
        data["CHROM"].astype(str)
        + ":"
        + data["POS"].astype(str)
        + "-"
        + data["REF"]
        + "-"
        + data["ALT"]
    )
    logger.info("New column `Locus` created using CHROM:POS-REF-ALT")
    return data


def build_cDNA_MANE(data: pd.DataFrame):
    """
    Create cDNA_NAME concatenating HGVSc + MANE_SELECT columns
    """
    colname = "CDNA_NAME"
    if colname in data.columns:
        logger.info(f"Column `{colname}` already exists, skipping creation")
        return data
    data['HGVSc'] = data['HGVSc'].astype(str)
    cdna = data["HGVSc"].str.extract(r"\w*(:c\..*)")
    data.insert(4, colname, data["MANE_SELECT"] + cdna.iloc[:, 0])
    logger.info(
        f"New column `{colname}` created using MANE_SELECT+HGVSc"
    )
    return data


def build_qc_status(data: pd.DataFrame) -> pd.DataFrame:
    """
    Evaluates pure VCF standard and variant caller quality metrics to assign QC_STATUS and QC_DETAIL.
    """
    logger.info("Building pure VCF caller QC_STATUS and QC_DETAIL columns...")

    def evaluate_qc(row):
        qual = row.get("QUAL", np.nan)
        v_filter = str(row.get("FILTER", "")).strip()
        qd = row.get("QD", np.nan)
        fs = row.get("FS", np.nan)
        mq = row.get("MQ", np.nan)
        
        qual_valid = pd.notna(qual) and str(qual) not in (".", "nan", "None", "")
        filter_valid = v_filter not in ("", ".", "nan", "None")
        qd_valid = pd.notna(qd) and str(qd) not in (".", "nan", "None", "")
        fs_valid = pd.notna(fs) and str(fs) not in (".", "nan", "None", "")
        mq_valid = pd.notna(mq) and str(mq) not in (".", "nan", "None", "")
        
        has_caller_qc = qual_valid or filter_valid or qd_valid or fs_valid or mq_valid
        
        if not has_caller_qc:
            return pd.Series(["UNSCORED", "No raw caller quality metrics present in VCF"])
        
        reasons = []
        is_low_qual = False
        
        if filter_valid and v_filter.upper() not in ("PASS", "."):
            is_low_qual = True
            reasons.append(f"FILTER={v_filter}")
            
        if qual_valid:
            try:
                q_val = float(qual)
                if q_val < 30.0:
                    is_low_qual = True
                    reasons.append(f"QUAL={q_val:.1f}<30")
            except ValueError:
                pass
                
        if qd_valid:
            try:
                qd_val = float(qd)
                if qd_val < 2.0:
                    is_low_qual = True
                    reasons.append(f"QD={qd_val:.1f}<2.0")
            except ValueError:
                pass

        if fs_valid:
            try:
                fs_val = float(fs)
                if fs_val > 60.0:
                    is_low_qual = True
                    reasons.append(f"FS={fs_val:.1f}>60")
            except ValueError:
                pass

        if mq_valid:
            try:
                mq_val = float(mq)
                if mq_val < 40.0:
                    is_low_qual = True
                    reasons.append(f"MQ={mq_val:.1f}<40")
            except ValueError:
                pass

        if is_low_qual:
            return pd.Series(["LOW_QUAL", "; ".join(reasons)])
            
        return pd.Series(["PASS", "Passed all available caller quality metrics"])

    qc_df = data.apply(evaluate_qc, axis=1)
    data["QC_STATUS"] = qc_df[0]
    data["QC_DETAIL"] = qc_df[1]
    
    logger.info(f"QC_STATUS counts:\n{data['QC_STATUS'].value_counts().to_string()}")
    return data


def filter_genotype_quality(
    data: pd.DataFrame, 
    columns: list[str] | str = GENOTYPE_QUAL_COLUMNS, 
    genotype_quality_score: float = 30,
    append_lowqual: bool = True
):
    if columns is None:
        logger.warning("No Genotype Quality columns provided for filtering")
        return data

    if isinstance(columns, str):
        columns = data.filter(regex=columns).columns.tolist()
        logger.debug(f"Genotype columns filtered by regex {columns}")
    if not columns:
        logger.warning(f"No columns found for pattern {columns}, skipping quality filtering")
        return data

    gq_mask = data[columns].ge(genotype_quality_score).any(axis=1)
    data = data[gq_mask].copy()
    
    logger.info(
        f"Dataframe filtered by quality score > {genotype_quality_score} for columns {columns} "
        f"Variants after quality > [{genotype_quality_score}] filtering: {data.shape[0]}"
    )

    if append_lowqual:
        data = _transform_low_gq(data, columns, gq_threshold=genotype_quality_score)

    return data


def _transform_low_gq(data, columns, gq_threshold: float = 10):
    for col in columns:
        if col.startswith("GQ_"):
            gt_col = col.replace("GQ_", "GT_")
            low_qual_mask = data[col] < gq_threshold
            data.loc[low_qual_mask, gt_col] = data[gt_col] + "_LOWQUAL"
            logger.info(
                f"Added `_LOWQUAL` to {gt_col} where {col} < {gq_threshold}: {data[gt_col].unique()}"
            )
    return data


def filter_variant_quality(
    data: pd.DataFrame, 
    qual_column: str = 'QUAL',
    variant_quality_score: float = 20
):
    if qual_column not in data.columns:
        logger.warning(f"Column {qual_column} not found in data, skipping variant quality filtering")
        return data
    qual_mask = data[qual_column].ge(variant_quality_score)
    data = data[qual_mask].copy()
    logger.info(
        f"Dataframe filtered by variant quality score > {variant_quality_score} for column {qual_column} "
        f"Variants after quality > [{variant_quality_score}] filtering: {data.shape[0]}"
    )
    return data


def select_columns(
    data: pd.DataFrame,
    column_file: str = None,
    formato: bool = False
):
    if column_file:
        cols = pd.read_csv(column_file, header=None)
        cols = cols[0].tolist()
        cols = [c.replace(",", "") for c in cols]
        cols += [c + ".1" for c in cols if c + ".1" in data.columns]
        cols += [c + ".2" for c in cols if c + ".2" in data.columns]
        logger.info(f"Selecting columns from file {column_file}")
        logger.debug(f"Columns to keep: {cols}")
        if formato:
            sample_cols = [col.replace('GT_','') for col in data.columns if col.startswith("GT_")]
            if not sample_cols:
                logger.error("No sample columns found in data")
                raise ValueError("No sample columns found in data")
            cols += ["FORMAT"]
            cols += sample_cols
        missing_cols = [col for col in cols if col not in data.columns]
        if missing_cols:
            logger.warning(
                f"The following columns are missing in the data: {missing_cols}"
            )
        genotype_columns = ("GT_", "GQ_","AD_", "DP_")
        cols += [col for col in data.columns if col.startswith(genotype_columns)]
        cols = list(dict.fromkeys(cols))
        data = data[[col for col in cols if col in data.columns]]
    else:
        logger.info("No column file provided, keeping all columns")
        cols = data.columns.tolist()
    logger.info(f"Dataframe shape after selecting columns: {data.shape}")
    return data


def select_samples(data: pd.DataFrame, samples: list[str]):
    if samples is None:
        logger.warning("No samples provided for filtering")
        return data
    gt_columns = data.filter(regex="^GT_").columns
    sample_gt_columns = gt_columns.intersection(samples)
    if sample_gt_columns.empty:
        logger.warning(f"No GT_ columns matched the provided samples: {samples} - available GT_ columns: {gt_columns.tolist()}")
        logger.warning("Skipping sample filtering")
        return data
    sample_gt_data = data.loc[:, sample_gt_columns]
    data = data[sample_gt_data.notna().ne("").any(axis=1)]
    sample_gt_data = data.loc[:, sample_gt_columns]
    data = data[~sample_gt_data.eq("0/0").all(axis=1)]
    data = data[~sample_gt_data.eq(".").all(axis=1)]
    data = data[~sample_gt_data.isna().all(axis=1)]
    data = data[~sample_gt_data.eq("").all(axis=1)]
    data = data[~sample_gt_data.eq("MISSING").all(axis=1)]
    data = data[~sample_gt_data.eq("HOMREF").all(axis=1)]
    logger.info(
        f"Dataframe shape after removing rows with all GT_ columns 0/0, empty, NaN or '.': {data.shape}"
    )
    return data


def merge_sample_cols(df, pattern="GT_", remove=True):
    sample_columns = [col for col in df.columns if col.startswith(pattern)]

    def process_sample_values(row, pattern, gt_columns):
        valid = {}
        for col in gt_columns:
            if pattern in ["GQ","AD","DP"]:
                sample_name = col.replace(pattern, "")
                if "merged_GT" in row and sample_name not in row["merged_GT"]:
                    continue
            if pd.notna(row[col]):
                if (
                    pattern == 'GT'
                    and row[col] != "HOMREF"
                    and row[col] != "MISSING"
                ) or (pattern != 'GT' and row[col] != 0):
                    sample_name = col.replace(pattern, "")
                    valid[sample_name] = row[col]
        merged = ", ".join(f"{k}:{v}" for k, v in valid.items())
        count = len(valid)
        return pd.Series([merged, count])

    stripped_pattern = pattern.strip("_")
    logger.info(
        f"Merging {pattern} columns into {f'merged_{stripped_pattern}'} and {f'count_{stripped_pattern}'}"
    )

    df[[f"merged_{stripped_pattern}", f"count_{stripped_pattern}"]] = df.apply(
        lambda row: process_sample_values(row, stripped_pattern, sample_columns), axis=1
    )

    if remove:
        df = df.drop(columns=sample_columns)
        logger.debug(f"Removed columns {sample_columns}")
    return df


def remove_columns(data: pd.DataFrame, columns: str):
    if columns is None:
        logger.warning("No columns provided for removing")
        return data
    elif isinstance(columns, str):
        cols = data.filter(regex=columns).columns.tolist()
        logger.info(f"Columns removed by regex {columns}")
        data.drop(columns=cols, axis=1, inplace=True, errors="ignore")
    elif isinstance(columns, list):
        data = data.drop(columns=columns, errors="ignore")
    else:
        raise ValueError("columns must be a string or a list")
    return data


def sort_by(
    data: pd.DataFrame,
    by: dict[str, bool],
    na_pos: Literal["first", "last"] = "first",
):
    valid_by = {}
    for col, ascending in by.items():
        if col not in data.columns:
            logger.warning(f"Column {col} not found in dataframe, skipping sort by this column")
        else:
            valid_by[col] = ascending
    
    if not valid_by:
        logger.warning("No valid columns found for sorting, returning unsorted dataframe")
        return data
    
    sorted_df = data.sort_values(
        by=list(valid_by.keys()), ascending=list(valid_by.values()), na_position=na_pos
    )
    logger.info(f"Sorted dataframe by {list(valid_by.keys())}")
    return sorted_df


def order_columns(data: pd.DataFrame, columns: list[str]):
    gt_cols = [x for x in data.columns if x.startswith("GT_")]
    if gt_cols:
        columns += gt_cols
    gq_cols = [x for x in data.columns if x.startswith("GQ_")]
    if gq_cols:
        columns += gq_cols
    ad_cols = [x for x in data.columns if x.startswith("AD_")]
    if ad_cols:
        columns += ad_cols
    dp_cols = [x for x in data.columns if x.startswith("DP_")]
    if dp_cols:
        columns += dp_cols
    missing_cols = [col for col in columns if col not in data.columns]
    if missing_cols:
        logger.warning(f"The following columns are missing in the data: {missing_cols}")
        columns = [col for col in columns if col in data.columns]
    logger.debug(f"Ordering columns by {columns}")
    data = data[columns + [col for col in data.columns if col not in columns]]
    logger.info(f"Dataframe ordered by columns {columns}")
    return data


def save_data(data: pd.DataFrame, out_file: str, log: str | None = None):
    logger.info(f"Saving data with shape {data.shape} to {out_file}")

    if out_file.endswith(".xlsx"):
        logger.info(f"Saving data to {out_file} as xlsx")
        logger.warning(
            "Excel (.xlsx) files have a hard limit of 1,048,576 rows and can be slow/memory-intensive to write. "
            f"Your data has {data.shape[0]} rows. For larger datasets, consider using '.parquet' or '.feather' formats."
        )
        if data.shape[0] >= 1048576:
            logger.error(
                f"DataFrame has {data.shape[0]} rows, which exceeds Excel's maximum limit of 1,048,576 rows (including header)."
            )
            raise ValueError(
                f"Cannot save to .xlsx: Data has {data.shape[0]} rows, which exceeds the Excel limit of 1,048,576 rows."
            )
        data.to_excel(
            out_file,
            index=False,
            engine="openpyxl",
            header=True,
        )
    elif out_file.endswith(".tsv"):
        logger.info(f"Saving data to {out_file} as tsv")
        data.to_csv(
            out_file,
            sep="\t",
            index=False,
        )
    elif out_file.endswith(".csv"):
        logger.info(f"Output file {out_file} saved as csv")
        data.to_csv(
            out_file,
            sep=",",
            index=False,
        )
    elif out_file.endswith(".parquet") or out_file.endswith(".pq"):
        logger.info(f"Saving data to {out_file} as parquet")
        try:
            data.to_parquet(out_file, index=False)
        except ImportError as e:
            logger.error("Failed to save Parquet file. Please ensure 'pyarrow' or 'fastparquet' is installed in your conda environment.")
            raise e
    elif out_file.endswith(".feather") or out_file.endswith(".arrow"):
        logger.info(f"Saving data to {out_file} as feather")
        try:
            data.reset_index(drop=True).to_feather(out_file)
        except ImportError as e:
            logger.error("Failed to save Feather file. Please ensure 'pyarrow' is installed in your conda environment.")
            raise e
    else:
        logger.error(
            f"Output file {out_file} has an unsupported format, should be .xlsx, .tsv, .csv, .parquet, or .feather"
        )
        raise ValueError(
            f"Output file {out_file} has an unsupported format, should be .xlsx, .tsv, .csv, .parquet, or .feather"
        )
