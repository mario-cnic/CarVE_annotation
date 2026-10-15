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
    src/tools/derive_dbnsfp_dependency.py), a row keeps the list element at the position of its own
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
    collapsed = {c: 0 for c in cols}
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
            if pos is None:
                continue                      # the row's transcript is not in dbNSFP's list: nothing here is its own
            if len(items) == n:
                if items[pos] != ".":
                    new[c][i] = items[pos]
            elif len(items) == 1:
                # one value for a multi-transcript list: variant-level (or collapsed upstream); kept, and counted
                new[c][i] = v if v != "." else np.nan
                collapsed[c] += 1
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
    single = {c: n for c, n in collapsed.items() if n}
    if single:
        logger.info(f"dbNSFP aligned columns holding ONE value for a multi-transcript list (kept as-is): {single}")
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


NMD_ESCAPE = "escape_predicted"
NMD_PREDICTED = "nmd_predicted"
NMD_NOT_EVALUATED = "not_evaluated"
NMD_NOT_ANNOTATED = "not_annotated"
NMD_NOT_APPLICABLE = "not_applicable"


def parse_nmd(data: pd.DataFrame) -> pd.DataFrame:
    """Derive NMD_STATUS per row from the VEP NMD plugin flag.

    The plugin only flags escape (last exon, last 50 bp of the penultimate exon, first 100
    coding bases, intronless transcript) for stop_gained/frameshift/canonical splice rows,
    judged at the variant position, and stays empty when it cannot place the variant in the
    CDS. Empty therefore means NMD predicted only for stop_gained/frameshift rows with a
    complete CDS position; canonical splice rows without a flag are not evaluated.
    """
    conseq = data["Consequence"].fillna("").astype(str) if "Consequence" in data.columns else pd.Series("", index=data.index)
    truncating = conseq.str.contains(r"stop_gained|frameshift_variant")
    in_scope = truncating | conseq.str.contains(r"splice_donor_variant|splice_acceptor_variant")
    status = pd.Series(NMD_NOT_APPLICABLE, index=data.index, dtype=object)
    if "NMD" not in data.columns:
        status.loc[in_scope] = NMD_NOT_ANNOTATED
        data["NMD_STATUS"] = status
        logger.warning("NMD column absent: NMD_STATUS set to not_annotated for truncating/splice rows")
        return data
    escape = data["NMD"].fillna("").astype(str).str.contains("NMD_escaping_variant")
    cds = data["CDS_position"].fillna("").astype(str).str.strip() if "CDS_position" in data.columns else pd.Series("", index=data.index)
    cds_placed = ~cds.isin(["", "-", ".", "nan"]) & ~cds.str.contains("?", regex=False)
    status.loc[in_scope] = NMD_NOT_EVALUATED
    status.loc[truncating & cds_placed] = NMD_PREDICTED
    status.loc[in_scope & escape] = NMD_ESCAPE
    data["NMD_STATUS"] = status
    logger.info(f"New column `NMD_STATUS` created with value counts: {dict(status.value_counts())}")
    return data


REMM_SCORED = "scored"
REMM_NOT_COVERED = "not_covered"
REMM_NOT_APPLICABLE = "not_applicable"
REMM_NOT_ANNOTATED = "not_annotated"
REMM_ERROR = "error"


def parse_remm(data: pd.DataFrame) -> pd.DataFrame:
    """Derive ReMM_score (float) and ReMM_status per row from the VEP --custom ReMM column.

    ReMM is one score per genomic position, not per allele, so only SNVs are scored and every
    ALT at the same position gets the same value. Non-SNVs are not_applicable; an SNV with no
    value in the track is not_covered (never 0). Informational only: no effect on NEW_IMPACT
    or PRIORITY_TIER.
    """
    if "ReMM" not in data.columns or "VARIANT_CLASS" not in data.columns:
        data["ReMM_score"] = np.nan
        data["ReMM_status"] = REMM_NOT_ANNOTATED
        logger.warning("ReMM or VARIANT_CLASS column absent: ReMM_status set to not_annotated")
        return data
    raw = data["ReMM"].astype(object).where(data["ReMM"].notna(), "").astype(str).str.strip()
    empty = raw.isin(["", ".", "-", "nan"])
    score = pd.to_numeric(raw.mask(empty), errors="coerce")
    unparsed = ~empty & score.isna()
    snv = data["VARIANT_CLASS"].fillna("").astype(str).eq("SNV")
    status = pd.Series(REMM_NOT_APPLICABLE, index=data.index, dtype=object)
    status.loc[snv & empty] = REMM_NOT_COVERED
    status.loc[snv & score.notna()] = REMM_SCORED
    status.loc[snv & unparsed] = REMM_ERROR
    # VEP gives a 1-bp deletion the deleted base's score; blank the raw value on non-SNV rows.
    data["ReMM"] = data["ReMM"].where(snv)
    data["ReMM_score"] = score.where(snv)
    data["ReMM_status"] = status
    if unparsed.any():
        logger.warning(f"ReMM: {int(unparsed.sum())} rows with an unparsable value, e.g. {raw[unparsed].iloc[0]!r}")
    logger.info(f"New column `ReMM_status` created with value counts: {dict(status.value_counts())}")
    return data
