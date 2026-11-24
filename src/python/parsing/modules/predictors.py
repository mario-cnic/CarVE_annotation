"""
Missense Pathogenicity Ensembles, AlphaMissense, REVEL, and dbNSFP Predictors.
"""
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
    """Transform missense columns to float, capturing MAX value if ambiguous.
    `missense_cols` is not modified (callers pass the module-level MISSENSE_COLUMNS; MISC-43)."""
    for col in missense_cols:
        if col not in data.columns:
            logger.warning(f"Column {col} not found in data, won't be parsed")
    missense_cols = [c for c in missense_cols if c in data.columns]

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


NMD_ESCAPE = "escape_predicted"
NMD_PREDICTED = "nmd_predicted"
NMD_NOT_EVALUATED = "not_evaluated"
NMD_NOT_ANNOTATED = "not_annotated"
NMD_NOT_APPLICABLE = "not_applicable"

# A premature stop whose first nucleotide lies in the last NMD_PENULTIMATE_WINDOW_NT nucleotides of the
# penultimate exon escapes NMD (ClinGen PVS1 tree, Abou Tayoun 2018: "3'-most 50 nucleotides of the penultimate exon").
NMD_PENULTIMATE_WINDOW_NT = 50
# Stops within the first NMD_START_PROXIMAL_NT coding nucleotides are flagged start_proximal (possible
# re-initiation downstream); reported only, NMD_STATUS stays nmd_predicted.
NMD_START_PROXIMAL_NT = 100
NMD_RULE_COLUMNS = ["NMD_RULE", "NMD_PTC_BASIS"]

_AA = r"[A-Z][a-z]{2}"
_HGVSP_FS = re.compile(rf"^p\.\(?{_AA}(\d+){_AA}?fsTer(\d+|\?)")
_HGVSP_NONSENSE = re.compile(rf"^p\.\(?{_AA}(\d+)Ter\)?$")
_HGVSP_DELINS = re.compile(rf"^p\.\(?{_AA}(\d+)(?:_{_AA}\d+)?delins((?:{_AA})*?)Ter")
_HGVSP_INS = re.compile(rf"^p\.\(?{_AA}(\d+)_{_AA}\d+ins((?:{_AA})*?)Ter")
_HGVSC_RANGE = r"(\d+)(?:_(\d+))?"
_HGVSC_EDIT = re.compile(rf"^c\.{_HGVSC_RANGE}(delins|del|dup|ins)([ACGTN]*)$")
_HGVSC_SUB = re.compile(r"^c\.\d+[ACGTN]>[ACGTN]$")
_FIRST_INT = re.compile(r"^\s*(\d+)")


def load_transcript_exons(path: str) -> dict:
    """{ENST: (version, exon_cdna_ends tuple, cds_start_cdna, start_incomplete, end_incomplete)}
    from build_transcript_exon_table.py."""
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in ["transcript_id", "transcript_version", "exon_cdna_ends", "cds_start_cdna",
                           "start_incomplete", "end_incomplete"] if c not in table.columns]
    if missing:
        raise ValueError(f"Transcript exon table {path} lacks columns {missing}")
    if table.transcript_id.duplicated().any():
        raise ValueError(f"Transcript exon table {path} repeats {table.transcript_id[table.transcript_id.duplicated()].iloc[0]}")
    return {
        r.transcript_id: (r.transcript_version, tuple(int(x) for x in r.exon_cdna_ends.split(",")),
                          int(r.cds_start_cdna), r.start_incomplete == "1", r.end_incomplete == "1")
        for r in table.itertuples(index=False)
    }


def _first_int(value) -> int | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    m = _FIRST_INT.match(str(value))
    return int(m.group(1)) if m and "?" not in str(value).split("-")[0] else None


def _hgvs_part(value) -> str:
    text = "" if value is None or (isinstance(value, float) and np.isnan(value)) else str(value)
    return text.split(":", 1)[1] if ":" in text else text


def _stop_codon_aa(hgvsp: str) -> tuple[int | None, str]:
    """(position of the new stop in the variant protein, kind); kind is 'stop', 'no_stop' or 'unknown'."""
    m = _HGVSP_FS.match(hgvsp)
    if m:
        return (None, "no_stop") if m.group(2) == "?" else (int(m.group(1)) + int(m.group(2)) - 1, "stop")
    m = _HGVSP_NONSENSE.match(hgvsp)
    if m:
        return int(m.group(1)), "stop"
    m = _HGVSP_DELINS.match(hgvsp)
    if m:
        return int(m.group(1)) + len(m.group(2)) // 3, "stop"
    m = _HGVSP_INS.match(hgvsp)
    if m:
        return int(m.group(1)) + 1 + len(m.group(2)) // 3, "stop"
    return None, "unknown"


def _length_change(hgvsc: str) -> int | None:
    """Coding length change (alt - ref) of an exonic HGVSc edit; None for intronic/UTR or unparsed edits."""
    if _HGVSC_SUB.match(hgvsc):
        return 0
    m = _HGVSC_EDIT.match(hgvsc)
    if not m:
        return None
    start, end, edit, seq = int(m.group(1)), int(m.group(2) or m.group(1)), m.group(3), m.group(4)
    span = end - start + 1
    if edit in ("ins", "delins") and not seq:
        return None
    return {"del": -span, "dup": span, "ins": len(seq), "delins": len(seq) - span}[edit]


def _nmd_rule(row, exons: dict) -> tuple[str, str, str]:
    """(NMD_STATUS, NMD_RULE, NMD_PTC_BASIS) for one stop_gained/frameshift row."""
    model = exons.get(str(row.Feature).split(".")[0])
    if model is None:
        return NMD_NOT_EVALUATED, "transcript_not_in_model", ""
    version, ends, cds_start, start_incomplete, end_incomplete = model
    hgvsc_id = str(row.HGVSc).split(":", 1)[0] if ":" in str(row.HGVSc) else ""
    if not hgvsc_id:
        return NMD_NOT_EVALUATED, "hgvs_missing", ""
    if hgvsc_id.partition(".")[2] != version:
        return NMD_NOT_EVALUATED, "transcript_version_mismatch", ""
    if end_incomplete:
        return NMD_NOT_EVALUATED, "incomplete_3prime", ""
    cdna, cds = _first_int(row.cDNA_position), _first_int(row.CDS_position)
    if cdna is None or cds is None:
        return NMD_NOT_EVALUATED, "position_unknown", ""
    if cdna - cds != cds_start - 1:
        return NMD_NOT_EVALUATED, "model_mismatch", ""
    penultimate_end = ends[-2] if len(ends) > 1 else None
    stop_aa, kind = _stop_codon_aa(_hgvs_part(row.HGVSp))
    if kind == "no_stop":
        if penultimate_end is None or cdna > penultimate_end:
            return NMD_ESCAPE, "last_exon" if penultimate_end else "intronless", "variant"
        return NMD_NOT_EVALUATED, "no_stop_in_new_frame", ""
    delta = _length_change(_hgvs_part(row.HGVSc))
    if kind == "unknown" or delta is None:
        return NMD_NOT_EVALUATED, "ptc_position_unknown", ""
    # first nucleotide of the new stop codon, moved from variant-transcript to reference coordinates
    ptc_cds = 3 * (stop_aa - 1) + 1 - delta
    ptc_cdna = cds_start - 1 + ptc_cds
    if penultimate_end is None:
        return NMD_ESCAPE, "intronless", "ptc"
    if ptc_cdna > penultimate_end:
        return NMD_ESCAPE, "last_exon", "ptc"
    if ptc_cdna > penultimate_end - NMD_PENULTIMATE_WINDOW_NT:
        return NMD_ESCAPE, "penultimate_last_50nt", "ptc"
    if ptc_cds <= NMD_START_PROXIMAL_NT and not start_incomplete:
        return NMD_PREDICTED, "start_proximal", "ptc"
    return NMD_PREDICTED, "none", "ptc"


def parse_nmd(data: pd.DataFrame, exons: dict | None = None) -> pd.DataFrame:
    """Derive NMD_STATUS, NMD_RULE and NMD_PTC_BASIS for stop_gained/frameshift/canonical splice rows.

    The rules are judged at the premature stop codon: for a frameshift its position comes from HGVSp
    (fsTer<n>) and the row transcript's exon structure (`exons`, load_transcript_exons). NMD_STATUS keeps
    its values: escape_predicted (NMD_RULE last_exon, penultimate_last_50nt or intronless), nmd_predicted
    (NMD_RULE none or start_proximal), not_evaluated (NMD_RULE gives the reason; canonical splice rows are
    splice_outcome_unknown), not_annotated (no exon table), not_applicable. NMD_PTC_BASIS: ptc (judged at
    the new stop) or variant (fsTer? in the last exon). The VEP NMD plugin column is left as it is.
    """
    conseq = data["Consequence"].fillna("").astype(str) if "Consequence" in data.columns else pd.Series("", index=data.index)
    truncating = conseq.str.contains(r"stop_gained|frameshift_variant")
    in_scope = truncating | conseq.str.contains(r"splice_donor_variant|splice_acceptor_variant")
    status = pd.Series(NMD_NOT_APPLICABLE, index=data.index, dtype=object)
    rule = pd.Series("", index=data.index, dtype=object)
    basis = pd.Series("", index=data.index, dtype=object)
    if exons is None:
        status.loc[in_scope] = NMD_NOT_ANNOTATED
        rule.loc[in_scope] = "no_transcript_model"
        logger.warning("No transcript exon table: NMD_STATUS set to not_annotated for truncating/splice rows")
    else:
        status.loc[in_scope] = NMD_NOT_EVALUATED
        rule.loc[in_scope] = "splice_outcome_unknown"
        if truncating.any():
            cols = ["Feature", "HGVSc", "HGVSp", "cDNA_position", "CDS_position"]
            sub = data.loc[truncating].reindex(columns=cols)
            res = [_nmd_rule(r, exons) for r in sub.itertuples(index=False)]
            status.loc[truncating] = [r[0] for r in res]
            rule.loc[truncating] = [r[1] for r in res]
            basis.loc[truncating] = [r[2] for r in res]
    data["NMD_STATUS"] = status
    data["NMD_RULE"] = rule
    data["NMD_PTC_BASIS"] = basis
    logger.info(f"NMD_STATUS counts: {dict(status.value_counts())}; NMD_RULE counts: {dict(rule[in_scope].value_counts())}")
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
    """Worst class over the row's uORF effects. Entries carrying `consequence=` (UTRAnnotatorPGP) are
    graded only under their own term; entries without it are graded under every term."""
    terms = [t.strip().lower() for t in re.split(r"[&,]", str(consequence)) if t.strip()]
    entries = _utr5_entries(annotation)
    linked = any("consequence" in e for e in entries)
    classes = []
    for term in terms:
        if any(k in term for k in _UTR5_NEVER_HIGH):
            classes.append(UTR5_OTHER)
            continue
        own = [e for e in entries if e.get("consequence", "").lower() == term] if linked else entries
        classes += [_utr5_entry_class(term, e) for e in own] or [UTR5_UNCLASSIFIED]
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
