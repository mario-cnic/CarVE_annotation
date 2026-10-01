"""
Missense Pathogenicity Ensembles, AlphaMissense, REVEL, and dbNSFP Predictors.
"""
import numpy as np
import pandas as pd
from . import logger

DBNSFP_COLUMNS = [
    "MetaRNN_score",
    "MetaLR_score",
    "REVEL_score",
    "INFO",
    "APPRIS",
    "Ensembl_proteinid",
    "GENCODE_basic",
    "HGVSc_ANNOVAR",
    "HGVSc_VEP",
    "HGVSc_snpEff",
    "MPC_score",
    "MVP_score",
    "SIFT4G_pred",
    "SIFT4G_score",
    "TSL",
    "VEP_canonical",
    "VEST4_score",
    "gMVP_score",
    "Aloft_Confidence",
    "Aloft_Fraction_transcripts_affected",
    "Aloft_pred",
    "Aloft_prob_Dominant",
    "Aloft_prob_Recessive",
    "Aloft_prob_Tolerant",
    "FATHMM_pred",
    "FATHMM_score",
    "GENCODE_basic",
    "Interpro_domain",
    "MPC_score",
    "MVP_score",
    "MutationAssessor_score",
    "MutationTaster_score",
    "PROVEAN_pred",
    "PROVEAN_score",
    "Polyphen2_HDIV_pred",
    "Polyphen2_HDIV_score",
    "Polyphen2_HVAR_pred",
    "Polyphen2_HVAR_score",
    "REVEL_score",
    "SIFT_pred",
    "SIFT_score",
    "TSL",
    "Uniprot_acc",
    "Uniprot_entry",
    "aapos",
    "genename"
]

MISSENSE_COLUMNS = [
    "MetaRNN_score",
    "MetaLR_score",
    "REVEL_score"
]


def parse_missense(data, missense_cols: list[str]):
    """Transform missense columns to float, capturing MAX value if ambiguous."""
    for col in list(missense_cols):
        if col not in data.columns:
            logger.warning(f"Column {col} not found in data, won't be parsed")
            missense_cols.remove(col)

    def _parse_missense_value(x):
        if isinstance(x, float):
            return x
        if isinstance(x, str):
            numeric_values = []
            for value in x.split("&"):
                if value == ".":
                    continue
                try:
                    numeric_values.append(float(value))
                except ValueError:
                    continue
            return max(numeric_values) if numeric_values else np.nan
        return np.nan

    try:
        if hasattr(pd.DataFrame, "map"):
            data[missense_cols] = data[missense_cols].map(_parse_missense_value)
        else:
            data[missense_cols] = data[missense_cols].applymap(_parse_missense_value)
        logger.info(f"Transformed columns {missense_cols} to float (captured MAX value if ambiguous)")
    except KeyError:
        logger.warning("Missense pred columns couldn't be found in file, continuing...")

    return data


def parse_dbnsfp(data: pd.DataFrame, transcript: str, dbnsfp_cols: str):
    """Extract corresponding value from dbnsfp columns for a specified transcript/gene."""
    if not transcript:
        logger.warning("Transcript not provided, using the Feature column")
        transcript = data['Feature'].values if 'Feature' in data.columns else None
        if isinstance(transcript, np.ndarray) and len(transcript) == 1:
            transcript = transcript[0]
            logger.info(f"Using transcript {transcript} from Feature column for dbnsfp parsing")
        else:
            logger.warning("No transcript information found in Feature column, skipping dbnsfp parsing")
            return data
        if not transcript:
            logger.warning("Transcript information is empty, skipping dbnsfp parsing")
            return data

    protein = ""
    if transcript:
        protein_vals = data[data['Feature'] == transcript]['ENSP'].values
        if len(protein_vals) > 0:
            protein = protein_vals[0]
            logger.info(f"Transcript {transcript} corresponds to protein {protein}")
        elif len(protein_vals) == 0:
            logger.warning(f"Transcript {transcript} not found in Ensembl_transcriptid column, skipping dbnsfp parsing")
            return data
        elif len(protein_vals) > 1:
            logger.warning(f"Transcript {transcript} corresponds to multiple proteins {protein_vals}, skipping dbnsfp parsing")
            return data
    
    dbnsfp_cols_df = pd.read_csv(dbnsfp_cols, sep='\t')

    def extract_dbnsfp_value(row, col: str, target_id: str, id_column: str):
        col_val = row.get(col)
        id_val = row.get(id_column)
        if not target_id:
            if id_column == 'Ensembl_transcriptid':
                target_id = row.get('Feature') if 'Feature' in row else None
            elif id_column == 'Ensembl_proteinid':
                target_id = row.get('ENSP') if 'ENSP' in row else None
        if not isinstance(col_val, str) or not isinstance(id_val, str):
            return col_val if not pd.isna(col_val) else np.nan
            
        id_list = id_val.split("&")
        val_list = col_val.split("&")
        
        if target_id in id_list:
            pos = id_list.index(target_id)
            if pos < len(val_list) and val_list[pos] != ".":
                return val_list[pos]
                
        if all(v == "." for v in val_list):
            return np.nan
        return col_val

    for _, col_info in dbnsfp_cols_df.iterrows():
        col = col_info.get("column_name")
        depends_on = str(col_info.get("dependency_type", "transcript")).lower()
        
        if pd.isna(col) or col not in data.columns:
            logger.warning(f"Column {col} not found in data, won't be parsed")
            continue

        target_id = (protein if "proteinid" in depends_on else transcript) if transcript else None
        id_column = "Ensembl_proteinid" if "proteinid" in depends_on else "Ensembl_transcriptid"
        
        if id_column not in data.columns:
            logger.warning(f"{id_column} not in data, skipping {col}")
            continue

        try:
            data[col] = data.apply(
                lambda row: extract_dbnsfp_value(row, col, target_id, id_column),
                axis=1
            )
            logger.debug(f"Extracted dbnsfp column {col} values using {id_column}")
        except Exception as e:
            logger.warning(f"Error parsing column {col}: {e}")

    logger.info("Finished parsing dbnsfp columns")
    return data
