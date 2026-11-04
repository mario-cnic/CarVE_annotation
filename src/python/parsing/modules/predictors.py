"""
Missense Pathogenicity Ensembles, AlphaMissense, REVEL, and dbNSFP Predictors.
"""
import os
import re

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
    
    # Reached only on a one-row frame (PRED-2). The CLI default path is not shipped and
    # TSV_TO_TABLE does not pass one, so a missing file skips this legacy step instead of failing.
    if not dbnsfp_cols or not os.path.isfile(dbnsfp_cols):
        logger.warning(f"dbNSFP column file {dbnsfp_cols!r} not found, skipping dbnsfp parsing")
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


AVI_SCORED = "scored"
AVI_NOT_COVERED = "not_covered"
AVI_NOT_APPLICABLE = "not_applicable"
AVI_NOT_ANNOTATED = "not_annotated"
AVI_ERROR = "error"


def parse_avi(data: pd.DataFrame) -> pd.DataFrame:
    """Convert the VEP AVI plugin columns AVI_PHRED and AVI_RAW to floats and derive AVI_status per row.

    AlphaGenome AVI scores exist for SNVs only, one per allele. The plugin leaves both columns empty
    for a non-SNV (not_applicable) and for an SNV absent from the file (not_covered, never 0).
    Informational only: no effect on NEW_IMPACT or PRIORITY_TIER.
    """
    if "AVI_PHRED" not in data.columns or "VARIANT_CLASS" not in data.columns:
        data["AVI_PHRED"] = np.nan
        data["AVI_RAW"] = np.nan
        data["AVI_status"] = AVI_NOT_ANNOTATED
        logger.warning("AVI_PHRED or VARIANT_CLASS column absent: AVI_status set to not_annotated")
        return data
    snv = data["VARIANT_CLASS"].fillna("").astype(str).eq("SNV")
    values = {}
    unparsed = pd.Series(False, index=data.index)
    for col in ["AVI_PHRED", "AVI_RAW"]:
        raw = data.get(col, pd.Series(np.nan, index=data.index))
        raw = raw.astype(object).where(raw.notna(), "").astype(str).str.strip()
        empty = raw.isin(["", ".", "-", "nan"])
        values[col] = pd.to_numeric(raw.mask(empty), errors="coerce")
        bad = ~empty & values[col].isna()
        if bad.any():
            logger.warning(f"{col}: {int(bad.sum())} rows with an unparsable value, e.g. {raw[bad].iloc[0]!r}")
        unparsed |= bad
    status = pd.Series(AVI_NOT_APPLICABLE, index=data.index, dtype=object)
    status.loc[snv & values["AVI_PHRED"].isna()] = AVI_NOT_COVERED
    status.loc[snv & values["AVI_PHRED"].notna()] = AVI_SCORED
    status.loc[snv & unparsed] = AVI_ERROR
    ok = snv & ~unparsed
    data["AVI_PHRED"] = values["AVI_PHRED"].where(ok)
    data["AVI_RAW"] = values["AVI_RAW"].where(ok)
    data["AVI_status"] = status
    logger.info(f"New column `AVI_status` created with value counts: {dict(status.value_counts())}")
    return data


UTR5_HIGH_IMPACT = "high_impact"
UTR5_OTHER = "other"
UTR5_UNCLASSIFIED = "unclassified"
UTR5_NO_EFFECT = "no_uorf_effect"
UTR5_NOT_EVALUATED = "not_evaluated"
UTR5_NOT_APPLICABLE = "not_applicable"
UTR5_NOT_ANNOTATED = "not_annotated"

_UTR5_OORF = {"inframe_oorf", "outofframe_oorf"}
_UTR5_KOZAK_HIGH = {"strong", "moderate"}
_UTR5_KOZAK_KNOWN = _UTR5_KOZAK_HIGH | {"weak"}
_UTR5_NEVER_HIGH = ("premature_start_codon_loss", "uorf_stop_codon_gain")


def _utr5_entries(annotation) -> list[dict]:
    """Split a 5UTR_annotation string into one key=value dict per affected uORF ('&' or ',' between uORFs)."""
    entries = []
    for chunk in re.split(r"[&,]", str(annotation)):
        fields = {}
        for item in chunk.split(":"):
            key, sep, value = item.partition("=")
            if sep:
                fields[key.strip()] = value.strip()
        if fields:
            entries.append(fields)
    return entries


def _utr5_start_or_evidence(entry: dict) -> bool | None:
    """True if the uORF start has a strong/moderate Kozak or translation evidence, False if neither, None if unknown."""
    kozak = entry.get("KozakStrength", "").lower()
    evidence = entry.get("Evidence", "").lower()
    if kozak in _UTR5_KOZAK_HIGH or evidence == "true":
        return True
    if kozak == "weak" and evidence == "false":
        return False
    return None


def _utr5_entry_class(consequence: str, entry: dict) -> str:
    """High-impact uORF-perturbing variant as defined by Whiffin et al. 2020 (Nat Commun 11:2523, Methods)."""
    if "premature_start_codon_gain" in consequence:
        if entry.get("type", "").lower() not in _UTR5_OORF:
            return UTR5_OTHER if entry.get("type") else UTR5_UNCLASSIFIED
        kozak = entry.get("KozakStrength", "").lower()
        if kozak in _UTR5_KOZAK_HIGH:
            return UTR5_HIGH_IMPACT
        return UTR5_OTHER if kozak in _UTR5_KOZAK_KNOWN else UTR5_UNCLASSIFIED
    if "uorf_stop_codon_loss" in consequence:
        alt_stop = entry.get("AltStop", "").lower()
        if alt_stop == "true":
            return UTR5_OTHER
        if alt_stop != "false":
            return UTR5_UNCLASSIFIED
        strong = _utr5_start_or_evidence(entry)
        return {True: UTR5_HIGH_IMPACT, False: UTR5_OTHER, None: UTR5_UNCLASSIFIED}[strong]
    if "uorf_frameshift" in consequence:
        ref_type = entry.get("ref_type", "").lower()
        alt_type = entry.get("alt_type", "").lower()
        if not ref_type or not alt_type:
            return UTR5_UNCLASSIFIED
        if ref_type in _UTR5_OORF or alt_type not in _UTR5_OORF:
            return UTR5_OTHER
        strong = _utr5_start_or_evidence(entry)
        return {True: UTR5_HIGH_IMPACT, False: UTR5_OTHER, None: UTR5_UNCLASSIFIED}[strong]
    return UTR5_UNCLASSIFIED


def _utr5_row_class(consequence, annotation) -> str:
    terms = [t.strip().lower() for t in re.split(r"[&,]", str(consequence)) if t.strip()]
    classes = []
    for term in terms:
        if any(k in term for k in _UTR5_NEVER_HIGH):
            classes.append(UTR5_OTHER)
            continue
        entries = _utr5_entries(annotation)
        classes += [_utr5_entry_class(term, e) for e in entries] or [UTR5_UNCLASSIFIED]
    for value in (UTR5_HIGH_IMPACT, UTR5_UNCLASSIFIED, UTR5_OTHER):
        if value in classes:
            return value
    return UTR5_UNCLASSIFIED


def parse_utrannotator(data: pd.DataFrame) -> pd.DataFrame:
    """Derive UTRAnnotator_status per row from the VEP UTRAnnotator columns.

    high_impact follows Whiffin et al. 2020: a uAUG gain that forms an in-frame or out-of-frame
    overlapping ORF (oORF) from a strong/moderate Kozak start; a uORF stop loss with no other stop
    in the UTR, from a strong/moderate Kozak start or a uORF with translation evidence; and, by the
    same rule, a frameshift that turns a uORF into an oORF. uAUG losses, uSTOP gains and other
    changes are `other`; `unclassified` when a deciding field is missing. 5'UTR rows the plugin
    evaluated without a uORF change are no_uorf_effect; not_evaluated when the plugin left no output.
    """
    if "5UTR_consequence" not in data.columns:
        data["UTRAnnotator_status"] = UTR5_NOT_ANNOTATED
        logger.warning("5UTR_consequence column absent: UTRAnnotator_status set to not_annotated")
        return data
    conseq = data["5UTR_consequence"].astype(object).where(data["5UTR_consequence"].notna(), "").astype(str).str.strip()
    has_effect = ~conseq.isin(["", "-", ".", "nan"])
    vep_conseq = data["Consequence"].fillna("").astype(str) if "Consequence" in data.columns else pd.Series("", index=data.index)
    utr5 = vep_conseq.str.contains("5_prime_UTR_variant", regex=False)
    evaluated = data["Existing_uORFs"].notna() if "Existing_uORFs" in data.columns else pd.Series(False, index=data.index)
    status = pd.Series(UTR5_NOT_APPLICABLE, index=data.index, dtype=object)
    status.loc[utr5] = UTR5_NOT_EVALUATED
    status.loc[utr5 & evaluated] = UTR5_NO_EFFECT
    if has_effect.any():
        annotation = data.loc[has_effect, "5UTR_annotation"] if "5UTR_annotation" in data.columns else pd.Series("", index=data.index[has_effect])
        status.loc[has_effect] = [
            _utr5_row_class(c, "" if pd.isna(a) else a) for c, a in zip(conseq[has_effect], annotation)
        ]
    data["UTRAnnotator_status"] = status
    logger.info(f"New column `UTRAnnotator_status` created with value counts: {dict(status.value_counts())}")
    return data


CONSTRAINT_MATCHED = "matched"
CONSTRAINT_NOT_IN_TABLE = "not_in_table"
CONSTRAINT_NOT_APPLICABLE = "not_applicable"
CONSTRAINT_NOT_ANNOTATED = "not_annotated"
# Output column -> column of the gene constraint table (build_gene_constraint_table.py)
CONSTRAINT_VALUES = {"pLI_gene_value": "pLI", "LOEUF_gene_value": "LOEUF", "mis_z_gene_value": "mis_z"}
CONSTRAINT_TEXT = {"gene_constraint_transcript": "transcript", "gene_constraint_flags": "constraint_flags"}


def load_gene_constraint(path: str) -> pd.DataFrame:
    """Gene constraint table indexed by Ensembl gene ID; refuses a table with a repeated gene."""
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in ["gene_id", *CONSTRAINT_VALUES.values(), *CONSTRAINT_TEXT.values()] if c not in table.columns]
    if missing:
        raise ValueError(f"Gene constraint table {path} lacks columns {missing}")
    if table.gene_id.duplicated().any():
        raise ValueError(f"Gene constraint table {path} repeats gene_id {table.gene_id[table.gene_id.duplicated()].iloc[0]}")
    for col in CONSTRAINT_VALUES.values():
        table[col] = pd.to_numeric(table[col].replace("", np.nan), errors="raise").astype("float64")
    return table.set_index("gene_id")


def annotate_gene_constraint(data: pd.DataFrame, table: pd.DataFrame | None) -> pd.DataFrame:
    """Attach gene-level pLI, LOEUF and missense Z to each row by its VEP `Gene` (Ensembl gene ID).

    Values are gnomAD's for the gene's MANE Select transcript (canonical when the gene has none); a
    gene without an estimate keeps an empty value, never 0. gene_constraint_status: matched,
    not_in_table (ENSG absent from the table), not_applicable (row without a gene) or not_annotated.
    """
    if "pLI_gene_value" in data.columns:
        logger.warning("Input carries pLI_gene_value (VEP pLI plugin); replaced by the gene constraint table")
    gene = data["Gene"].astype(object).where(data["Gene"].notna(), "").astype(str).str.strip() if "Gene" in data.columns else pd.Series("", index=data.index)
    has_gene = gene.str.startswith("ENSG")
    status = pd.Series(CONSTRAINT_NOT_APPLICABLE, index=data.index, dtype=object)
    if table is None:
        status.loc[has_gene] = CONSTRAINT_NOT_ANNOTATED
        for col in CONSTRAINT_VALUES:
            data[col] = np.nan
        for col in CONSTRAINT_TEXT:
            data[col] = None
        data["gene_constraint_status"] = status
        logger.warning("No gene constraint table: gene_constraint_status set to not_annotated for rows with a gene")
        return data
    matched = has_gene & gene.isin(table.index)
    status.loc[has_gene] = CONSTRAINT_NOT_IN_TABLE
    status.loc[matched] = CONSTRAINT_MATCHED
    looked_up = table.reindex(gene.where(matched))
    for out, col in CONSTRAINT_VALUES.items():
        data[out] = looked_up[col].to_numpy(dtype="float64")
    for out, col in CONSTRAINT_TEXT.items():
        data[out] = pd.Series(looked_up[col].to_numpy(dtype=object), index=data.index).where(matched, None)
    data["gene_constraint_status"] = status
    logger.info(f"New column `gene_constraint_status` created with value counts: {dict(status.value_counts())}")
    return data


def utr5_high_impact_mask(data: pd.DataFrame) -> pd.Series:
    """Rows whose UTRAnnotator_status is high_impact; derives the status first when the column is absent."""
    if "UTRAnnotator_status" not in data.columns:
        data = parse_utrannotator(data.copy())
    return data["UTRAnnotator_status"].eq(UTR5_HIGH_IMPACT)
