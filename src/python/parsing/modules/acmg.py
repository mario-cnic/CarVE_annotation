"""
Upgraded NEW_IMPACT and ClinVar significance extraction.
"""
import pandas as pd
from . import logger
from .io_qc import _safe_numeric_series, _safe_string_series
from .clinvar import annotate_clinvar, clinvar_masks
from .splicing import build_spliceMAX, _extract_spliceai_score, NEW_SPIP_COLUMNS
from .predictors import utr5_high_impact_mask




def _extract_clinvar_status(df: pd.DataFrame):
    """Extract source-selected ClinVar masks and expose selection provenance."""
    classification, _, _, _ = annotate_clinvar(df)
    return clinvar_masks(classification)


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

    if "splicevardb_classification" in data.columns:
        splicevardb_high_mask = data["splicevardb_classification"].notna() & (data["splicevardb_classification"].astype(str).str.strip() != "") & (data["splicevardb_classification"].astype(str) != "nan")
    else:
        splicevardb_high_mask = pd.Series(False, index=data.index)

    splicing_mod = spliceai_mod_mask | spip_mod_mask | pangolin_mod_mask | sv_mod_mask
    splicing_high = spliceai_high_mask | spip_high_mask | pangolin_high_mask | sv_high_mask | splicevardb_high_mask

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

    # 3. 5'UTR uORF: a high-impact change raises LOW/MODIFIER to MODERATE
    utr5_mod_mask = data["NEW_IMPACT"].isin(["LOW", "MODIFIER"]) & utr5_high_impact_mask(data)
    data.loc[utr5_mod_mask, "NEW_IMPACT"] = "MODERATE"

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
