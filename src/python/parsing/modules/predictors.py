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


DBNSFP_MATCHED = "matched_transcript"
DBNSFP_NO_ENTRY = "no_entry_for_row_transcript"
DBNSFP_NO_ROW_TRANSCRIPT = "no_row_transcript"
DBNSFP_NO_PREDICTION = "no_prediction"
DBNSFP_MISALIGNED = "misaligned_list"


def load_aligned_columns(path: str) -> list[str]:
    """Column names (one per line, '#' comments allowed) listed in resources/dbnsfp_transcript_aligned_columns.txt."""
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]


def parse_dbnsfp_by_row_transcript(data: pd.DataFrame, aligned_columns: list[str] = None) -> pd.DataFrame:
    """Give each row the dbNSFP value of its OWN transcript.

    VEP's dbNSFP plugin attaches to every transcript row of a variant the full '&'-joined lists of
    per-transcript values, aligned position by position with `Ensembl_transcriptid`. Here, for the
    transcript-aligned columns (resources/dbnsfp_transcript_aligned_columns.txt, derived from data by
    src/tools/derive_dbnsfp_alignment.py), a row keeps the list element at the position of its own
    `Feature`; if its transcript is not in dbNSFP's list the aligned columns are left empty
    (`dbNSFP_match` = no_entry_for_row_transcript). Columns that are single-valued, or multi-valued but
    not transcript-aligned, are untouched. Adds `dbNSFP_match`, `dbNSFP_transcripts_all` (the original
    transcript list) and `<col>_anytranscript_max` for the numeric MISSENSE_COLUMNS that are aligned.
    """
    tcol = "Ensembl_transcriptid"
    if tcol not in data.columns:
        logger.info("No dbNSFP transcript column, skipping per-transcript dbNSFP selection")
        return data
    if aligned_columns is None:
        raise ValueError("parse_dbnsfp_by_row_transcript requires the aligned-column list "
                         "(--dbnsfp-aligned-columns); refusing to leave per-transcript lists on every row")
    cols = [c for c in aligned_columns if c in data.columns]
    missing = data[tcol].isna() | data[tcol].astype(str).isin([".", ""])
    data["dbNSFP_match"] = DBNSFP_NO_PREDICTION
    data["dbNSFP_transcripts_all"] = data[tcol].where(~missing)
    work = data.index[~missing]
    if len(work) == 0:
        return data
    strip = lambda x: str(x).split(".")[0]
    sub = {c: data.loc[work, c].tolist() for c in cols}
    feats = data.loc[work, "Feature"].tolist() if "Feature" in data.columns else [None] * len(work)
    tlists = sub[tcol] if tcol in sub else data.loc[work, tcol].tolist()
    new = {c: [np.nan] * len(work) for c in cols}
    anymax = {c: [np.nan] * len(work) for c in cols if c in MISSENSE_COLUMNS}
    status = []
    for i, (feat, tl) in enumerate(zip(feats, tlists)):
        if feat is None or (isinstance(feat, float) and np.isnan(feat)):
            status.append(DBNSFP_NO_ROW_TRANSCRIPT)
            continue
        ids = [strip(x) for x in str(tl).split("&")]
        n = len(ids)
        pos = ids.index(strip(feat)) if strip(feat) in ids else None
        st = DBNSFP_MATCHED if pos is not None else DBNSFP_NO_ENTRY
        for c in cols:
            v = sub[c][i]
            if not isinstance(v, str):
                continue
            items = v.split("&")
            if c in anymax:
                nums = []
                for x in items:
                    try:
                        nums.append(float(x))
                    except ValueError:
                        pass
                anymax[c][i] = max(nums) if nums else np.nan
            if len(items) == 1:
                new[c][i] = v if v != "." else np.nan        # single value: not a per-transcript list
            elif len(items) == n:
                if pos is not None and items[pos] != ".":
                    new[c][i] = items[pos]
            else:
                st = DBNSFP_MISALIGNED if st == DBNSFP_MATCHED else st
        status.append(st)
    for c in cols:
        data[c] = data[c].astype(object)
        data.loc[work, c] = new[c]
    for c, vals in anymax.items():
        data[f"{c}_anytranscript_max"] = np.nan
        data.loc[work, f"{c}_anytranscript_max"] = vals
    data.loc[work, "dbNSFP_match"] = status
    logger.info(f"dbNSFP per-transcript selection: {pd.Series(status).value_counts().to_dict()}")
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
