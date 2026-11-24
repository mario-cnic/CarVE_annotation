"""
Interpretation inputs (gene curation, frequency status, unified SpliceAI), column ordering and frequency filters.
"""
import numpy as np
import pandas as pd
from . import logger
from .io_qc import build_qc_status
from .splicing import _extract_spliceai_details
from .population import annotate_population_af, resolve_population_af, rarity_af, ABSENT, GNOMAD_GRPMAX, GNOMAD_JOINT


def build_interpretation_inputs(data: pd.DataFrame, freq_col: str = GNOMAD_GRPMAX) -> pd.DataFrame:
    """Adds the annotation Module 3 scores from: population-frequency status and the unified SpliceAI
    score. Variant scores, tiers, ACMG tags, gene curation and inheritance are computed by Module 3."""
    if "QC_STATUS" not in data.columns:
        data = build_qc_status(data)

    annotate_population_af(data, freq_col)

    splice_max, splice_method = _extract_spliceai_details(data)
    data["SPLICE_MAX_UNIFIED"] = splice_max
    data["SPLICEAI_METHOD"] = splice_method
    splice_errors = data.get("SpliceAI_status", pd.Series("", index=data.index)).eq("error")
    splice_errors |= data.get("spliceai_custom_match", pd.Series("", index=data.index)).eq("predictor_error")
    data["SPLICEAI_STATUS"] = np.where(splice_max.notna(), "scored", np.where(splice_errors, "error", "not_covered"))

    return order_clinical_columns(data)


def order_clinical_columns(data: pd.DataFrame) -> pd.DataFrame:
    """Re-orders DataFrame columns so primary clinical triage tags are at the front."""
    front_cols = [
        "SYMBOL", "Locus", "HGVSc", "HGVSp", "Consequence", "NEW_IMPACT", "IMPACT",
        "GENE_TIER", "CLINGEN_VALIDITY", "GAMMA_WEIGHT", "PATHOGENIC_MECHANISM", "PVS1_ELIGIBILITY",
        "COMPOUND_HET_STATUS", "COMPOUND_HET_PAIR", "CRITICAL_DOMAINS_HOTSPOTS"
    ]
    
    evidence_cols = [
        "CLINVAR_SELECTED_CLNSIG", "CLINVAR_SOURCE", "CLINVAR_SELECTED_REVIEW", "ClinVar_Stars",
        "ClinVarVCF_CLNSIG", "CLINVAR_DISPLAY", "clinvar_clnsig", "CLNSIG",
        "gnomADv4_AF_grpmax_joint", "MAX_AF", "af", "gnomad_af", "AF",
        "pLI_gene_value", "pLI", "NMD_STATUS", "NMD_RULE", "NMD_PTC_BASIS", "NMD",
        "spliceai_custom_MAX", "spliceAI_MAX", "SPLICE_MAX_UNIFIED", "SPLICEAI_METHOD", "SPLICEAI_STATUS",
        "SPiP_prediction", "SPiP", "Pangolin_max_score", "Pangolin",
        "SpliceVault_status", "LaBranchoR_status",
        "am_pathogenicity", "AlphaMissense", "REVEL_score", "REVEL",
        "UTRAnnotator_status", "5UTR_consequence",
        "QC_STATUS", "QC_DETAIL"
    ]

    existing_cols = list(data.columns)
    
    ordered = []
    for c in front_cols:
        if c in existing_cols and c not in ordered:
            ordered.append(c)
            
    for c in existing_cols:
        if c not in ordered and c not in evidence_cols and not c.startswith(("GT_", "GQ_", "AD_", "DP_")):
            ordered.append(c)

    for c in evidence_cols:
        if c in existing_cols and c not in ordered:
            ordered.append(c)

    for c in existing_cols:
        if c not in ordered:
            ordered.append(c)

    return data[ordered]


def audit_annotation_availability(data: pd.DataFrame) -> dict[str, dict]:
    """Audits which bioinformatic tools, annotators, and predictors are present vs missing in the dataset."""
    tool_map = {
        "gnomAD Allele Frequency": [GNOMAD_GRPMAX, GNOMAD_JOINT],
        "SpliceAI Custom Window (-D 4999 / ~10kb)": ["spliceai_custom_MAX", "SpliceAI_custom_MAX", "spliceai_custom_DS_AG", "SpliceAI_custom_DS_AG", "SpliceAI_D10000", "SpliceAI_10k"],
        "SpliceAI VEP Standard (~50bp)": ["spliceAI_MAX", "SpliceAI_pred_DS_AG", "SpliceAI_DS_AG", "SpliceAI"],
        "SPiP Deep Splicing": ["SPiP_prediction", "SPiP"],
        "Pangolin Non-Coding": ["Pangolin_max_score", "Pangolin"],
        "AlphaMissense": ["am_pathogenicity", "AlphaMissense"],
        "REVEL Missense Ensemble": ["REVEL_score", "REVEL"],
        "SpliceVault RNA-Seq Cohort": ["SpliceVault_status", "SpliceVault"],
        "LaBranchoR branchpoint (positional, informational)": ["LaBranchoR_score"],
        "SpliceAI2 precomputed (informational)": ["SpliceAI2_score"],
        "ReMM Non-Coding (informational)": ["ReMM_score"],
        "AlphaGenome AVI (informational)": ["AVI_PHRED"],
        "ClinVar Pathogenicity": ["ClinVarVCF_CLNSIG", "clinvar_clnsig", "CLIN_SIG", "CLNSIG", "CLINVAR_DISPLAY", "ClinVar"],
        "pLI Gene LoF Intolerance": ["pLI_gene_value", "pLI"],
        "5' UTR uORF (UTRAnnotator)": ["5UTR_consequence", "UTRAnnotator_status"]
    }

    audit_results = {}
    cols_lower = {c.lower(): c for c in data.columns}

    for tool_name, candidates in tool_map.items():
        if tool_name == "gnomAD Allele Frequency":
            af, _, af_status = resolve_population_af(data)
            count = int((af.notna() | af_status.eq(ABSENT)).sum())
            audit_results[tool_name] = {
                "status": "Available" if count else "Unknown (no frequency evidence)",
                "column": "custom v4 grpmax / joint",
                "coverage": f"{count:,}/{len(data):,} ({count/max(1, len(data))*100:.1f}%)",
                "is_present": any(c in data.columns for c in candidates),
            }
            continue
        matched_col = None
        for cand in candidates:
            if cand in data.columns:
                matched_col = cand
                break
            elif cand.lower() in cols_lower:
                matched_col = cols_lower[cand.lower()]
                break

        if matched_col is not None:
            non_null = data[matched_col].notna().sum()
            audit_results[tool_name] = {
                "status": "Available",
                "column": matched_col,
                "coverage": f"{non_null:,}/{len(data):,} ({non_null/max(1, len(data))*100:.1f}%)",
                "is_present": True
            }
        else:
            audit_results[tool_name] = {
                "status": "Not Present (Baseline Fallback 0.0 Used)",
                "column": "None",
                "coverage": "0%",
                "is_present": False
            }

    return audit_results


def filter_freq(
    data: pd.DataFrame,
    freq,
    freq_column: str = "gnomADv4_AF_grpmax_joint",
):
    expected = [GNOMAD_GRPMAX, GNOMAD_JOINT] if freq_column == GNOMAD_GRPMAX else [freq_column]
    if not any(column in data.columns for column in expected):
        raise ValueError(f"Population frequency column(s) missing: {expected}")
    result = data.copy()
    af = rarity_af(annotate_population_af(result, freq_column), result["POPULATION_AF_STATUS"])
    unknown = af.isna()
    new_df = result.loc[unknown | (af < freq)].copy()
    absent = int(result["POPULATION_AF_STATUS"].eq(ABSENT).sum())
    logger.info(f"Frequency < {freq}: retained {len(new_df)} rows, including {absent} absent from gnomAD v4 and {int(unknown.sum())} with unknown AF")
    return new_df


def classification_dominant(data: pd.DataFrame, freq_col: str, freq: float):
    thresholds = {
        ("prioritary", "HIGH"): 0.01,
        ("prioritary", "MODERATE"): 1e-3,
        ("prioritary", "MODIFIER"): 1e-3,
        ("prioritary", "LOW"): 5e-4,
        ("secondary", "HIGH"): 1e-4,
        ("secondary", "MODERATE"): 1e-4,
        ("secondary", "LOW"): 5e-4,
        ("candidate", "HIGH"): 1e-4,
        ("candidate", "MODERATE"): 1e-4,
    }

    data = data.copy()
    annotate_population_af(data, freq_col)

    def evaluate_variant(row):
        prio = row["gene_priority"]
        impact = row["NEW_IMPACT"]
        freq_val = 0.0 if row["POPULATION_AF_STATUS"] == ABSENT else row["POPULATION_AF"]
        if pd.isna(freq_val):
            return "Unknown frequency"
        limit = thresholds.get((prio, impact))
        if limit and freq_val < limit:
            return "Candidate"
        return "Excluded"

    data["Variant_rank"] = data.apply(evaluate_variant, axis=1)
    logger.info("New column `Variant_rank` created based on gene priority, impact and frequency")
    return data
