"""
ACMG/ClinGen 2015 Criteria Tagging, Upgraded NEW_IMPACT, and ClinVar Significance Extraction.
"""
import numpy as np
import pandas as pd
from . import logger
from .io_qc import _safe_numeric_series, _safe_string_series
from .splicing import build_spliceMAX, _extract_spliceai_score, NEW_SPIP_COLUMNS


def _extract_clinvar_status(df: pd.DataFrame):
    """Extract boolean masks for ClinVar Pathogenic and Benign variants."""
    is_clinvar_path = pd.Series(False, index=df.index)
    is_clinvar_benign = pd.Series(False, index=df.index)
    for cln_col in ["CLINVAR_DISPLAY", "clinvar_clnsig", "CLNSIG", "ClinicalSignificance", "ClinVar", "clinvar"]:
        if cln_col in df.columns:
            cln_str = df[cln_col].astype(str).str.lower()
            is_clinvar_path = cln_str.str.contains("pathogenic") & (~cln_str.str.contains("conflicting|benign"))
            is_clinvar_benign = cln_str.str.contains("benign") & (~cln_str.str.contains("conflicting|pathogenic"))
            break
    return is_clinvar_path, is_clinvar_benign


def build_newImpact(data: pd.DataFrame):
    """
    Create upgraded NEW_IMPACT column with the most severe functional impact
    by integrating VEP base impact, deep splicing predictors, missense consensus, UTR disruption, and ClinVar.
    """
    data = build_spliceMAX(data)

    if "IMPACT" not in data.columns:
        data["IMPACT"] = "MODIFIER"
    data["NEW_IMPACT"] = data["IMPACT"].astype(str).copy()
    logger.debug(f"Initial NEW_IMPACT values from VEP: {data['NEW_IMPACT'].unique()}")

    # 1. SPLICING IMPACT
    spliceai_score = _extract_spliceai_score(data)
    spliceai_mod_mask = (data["NEW_IMPACT"] != "HIGH") & (spliceai_score > 0.20)
    spliceai_high_mask = spliceai_score > 0.50

    spip_exists = {
        NEW_SPIP_COLUMNS["prediction"],
        NEW_SPIP_COLUMNS["interpretation"]
    }.issubset(data.columns)
    
    if spip_exists:
        spip_pred = _safe_numeric_series(data, NEW_SPIP_COLUMNS["prediction"])
        spip_interp = _safe_string_series(data, NEW_SPIP_COLUMNS["interpretation"])
        spip_mod_mask = (data["NEW_IMPACT"] != "HIGH") & (spip_pred > 0.20) & (spip_interp != "NTR")
        spip_high_mask = (spip_pred > 0.50) & (spip_interp != "NTR")
    else:
        spip_mod_mask = pd.Series(False, index=data.index)
        spip_high_mask = pd.Series(False, index=data.index)

    pangolin_score = _safe_numeric_series(data, "Pangolin_max_score")
    pangolin_mod_mask = (data["NEW_IMPACT"] != "HIGH") & (pangolin_score >= 0.20)
    pangolin_high_mask = pangolin_score >= 0.50

    sv_status = _safe_string_series(data, "SpliceVault_status")
    sv_aberrant = sv_status == "aberrant_event_detected"
    sv_mod_mask = (data["NEW_IMPACT"] != "HIGH") & sv_aberrant
    sv_high_mask = sv_aberrant & (spliceai_score >= 0.20)

    bp_status = _safe_string_series(data, "Branchpoint_status")
    labranchor_score = _safe_numeric_series(data, "LaBranchoR_score")
    bp_high_mask = (bp_status == "disrupted") | (labranchor_score >= 0.50)

    if "splicevardb_classification" in data.columns:
        splicevardb_high_mask = data["splicevardb_classification"].notna() & (data["splicevardb_classification"].astype(str).str.strip() != "") & (data["splicevardb_classification"].astype(str) != "nan")
    else:
        splicevardb_high_mask = pd.Series(False, index=data.index)

    splicing_mod = spliceai_mod_mask | spip_mod_mask | pangolin_mod_mask | sv_mod_mask
    splicing_high = spliceai_high_mask | spip_high_mask | pangolin_high_mask | sv_high_mask | bp_high_mask | splicevardb_high_mask

    data.loc[splicing_mod, "NEW_IMPACT"] = "MODERATE"
    data.loc[splicing_high, "NEW_IMPACT"] = "HIGH"

    # 2. PATHOGENIC MISSENSE CONSENSUS
    am_score = _safe_numeric_series(data, "am_pathogenicity")
    am_class = _safe_string_series(data, "am_class").str.lower()
    revel_score = _safe_numeric_series(data, "REVEL_score")

    missense_mod_mask = (data["NEW_IMPACT"].isin(["LOW", "MODIFIER"])) & ((am_score >= 0.564) | (am_class == "pathogenic") | (revel_score >= 0.50))
    data.loc[missense_mod_mask, "NEW_IMPACT"] = "MODERATE"

    missense_high_mask = ((am_score >= 0.564) | (am_class == "pathogenic")) & (revel_score >= 0.75)
    data.loc[missense_high_mask, "NEW_IMPACT"] = "HIGH"

    # 3. UTR ANNOTATIONS
    utr_annot_cols = [
        "lost_stop_codon", "lost_start_codon", "mrl_gainedOrLost",
        "utr_num_kozak_gainedOrLost", "utr_num_uAUG_gainedOrLost",
        "num_polyA_signal_gainedOrLost", "tss_kozak_score_gainedOrLost"
    ]
    existing_utr_cols = [col for col in utr_annot_cols if col in data.columns]
    if existing_utr_cols:
        utr_values_of_interest = ["gained", "lost", "true"]
        utrann_high_mask = (
            data[existing_utr_cols]
            .astype(str)
            .apply(lambda col: col.str.lower())
            .isin(utr_values_of_interest)
            .any(axis=1)
        )
        data.loc[utrann_high_mask, "NEW_IMPACT"] = "HIGH"

    if "5UTR_consequence" in data.columns:
        utrant_high_mask = data["5UTR_consequence"].notna() & (data["5UTR_consequence"].astype(str).str.strip() != "") & (data["5UTR_consequence"].astype(str) != "nan")
        data.loc[utrant_high_mask, "NEW_IMPACT"] = "HIGH"

    # 4. CLINVAR PATHOGENIC EVIDENCE
    is_clinvar_path, _ = _extract_clinvar_status(data)
    data.loc[is_clinvar_path, "NEW_IMPACT"] = "HIGH"

    data["NEW_IMPACT"] = pd.Categorical(
        data["NEW_IMPACT"],
        categories=["HIGH", "MODERATE", "LOW", "MODIFIER"],
        ordered=True,
    )
    logger.info(f"New column `NEW_IMPACT` created with value counts: {dict(data['NEW_IMPACT'].value_counts())}")
    return data


def build_acmg_criteria(data: pd.DataFrame) -> pd.DataFrame:
    """Assigns explicit ACMG/ClinGen 2015 evidence code tags to variants."""
    logger.info("Computing explicit ACMG/ClinGen 2015 criteria codes...")
    conseq = _safe_string_series(data, "Consequence").str.lower()
    symbol_series = _safe_string_series(data, "SYMBOL").str.upper()
    af_series = _safe_numeric_series(data, "gnomADv4_AF_grpmax_joint") if "gnomADv4_AF_grpmax_joint" in data.columns else _safe_numeric_series(data, "MAX_AF")
    pli_val = _safe_numeric_series(data, "pLI_gene_value")

    is_lof = conseq.str.contains("stop_gained|frameshift|splice_donor|splice_acceptor")
    is_pvs1 = is_lof & (pli_val >= 0.9) & (symbol_series != "MYH7")

    is_clinvar_path, _ = _extract_clinvar_status(data)
    is_ps1 = is_clinvar_path & conseq.str.contains("missense")

    is_pm2 = af_series < 0.0001
    is_pm4 = conseq.str.contains("inframe_insertion|inframe_deletion|stop_lost|start_lost")

    splice_max = _extract_spliceai_score(data)
    spip_score = _safe_numeric_series(data, "SPiP_prediction") if "SPiP_prediction" in data.columns else _safe_numeric_series(data, "SPiP")
    am_score = _safe_numeric_series(data, "am_pathogenicity")
    revel_score = _safe_numeric_series(data, "REVEL_score")

    is_pp3 = (splice_max >= 0.50) | (spip_score >= 0.50) | (am_score >= 0.564) | (revel_score >= 0.75)
    is_bp4 = (splice_max < 0.10) & (spip_score < 0.10) & (am_score > 0) & (am_score < 0.34) & (revel_score < 0.25)
    is_bp7 = conseq.str.contains("synonymous") & (splice_max < 0.10)

    is_ba1 = af_series > 0.05
    is_bs1 = (af_series > 0.01) & (af_series <= 0.05)

    acmg_lists = []
    for idx in range(len(data)):
        codes = []
        if is_ba1.iloc[idx]: codes.append("BA1")
        elif is_bs1.iloc[idx]: codes.append("BS1")

        if is_pvs1.iloc[idx]: codes.append("PVS1")
        if is_ps1.iloc[idx]: codes.append("PS1")
        if is_pm2.iloc[idx]: codes.append("PM2_Supporting")
        if is_pm4.iloc[idx]: codes.append("PM4")
        if is_pp3.iloc[idx]: codes.append("PP3_Strong")
        if is_bp4.iloc[idx]: codes.append("BP4")
        if is_bp7.iloc[idx]: codes.append("BP7")

        acmg_lists.append(", ".join(codes) if codes else "Unclassified")

    data["ACMG_CRITERIA"] = acmg_lists
    logger.info("ACMG criteria assigned successfully.")
    return data
