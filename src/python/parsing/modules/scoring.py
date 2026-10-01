"""
Priority Score Calculation (0.0 - 100.0), Discrete Priority Tiers (Tier 1-4), Rationale Generation, and Column Ordering.
"""
import numpy as np
import pandas as pd
from . import logger
from .io_qc import _safe_numeric_series, _safe_string_series, build_qc_status
from .splicing import _extract_spliceai_details, NEW_SPIP_COLUMNS
from .disease_hpo import apply_disease_phenotype_curation, HPO_GENE_PANELS
from .acmg import build_acmg_criteria, _extract_clinvar_status


def build_priority_tier(
    data: pd.DataFrame,
    freq_col: str = "gnomADv4_AF_grpmax_joint",
    hpo_terms: list[str] | None = None,
    disease_phenotype: str = "cardiomyopathies_all"
) -> pd.DataFrame:
    """
    Build standardized ACMG/AMP and ClinGen-aligned variant priority tiers (PRIORITY_TIER),
    explicit ACMG evidence tags (ACMG_CRITERIA), disease curation evidence levels (DISEASE_GENE_EVIDENCE),
    and continuous prioritization score (VARIANT_PRIORITY_SCORE, 0.0 - 100.0).
    """
    logger.info(f"Calculating multi-evidence PRIORITY_TIER and VARIANT_PRIORITY_SCORE for disease phenotype '{disease_phenotype}'...")
    
    if "QC_STATUS" not in data.columns:
        data = build_qc_status(data)

    data = apply_disease_phenotype_curation(data, disease_phenotype)

    if freq_col in data.columns:
        af_series = _safe_numeric_series(data, freq_col)
    else:
        freq_candidates = [c for c in data.columns if "af" in c.lower() or "freq" in c.lower()]
        if freq_candidates:
            af_series = _safe_numeric_series(data, freq_candidates[0])
        else:
            af_series = pd.Series(0.0, index=data.index)

    conseq = _safe_string_series(data, "Consequence").str.lower()
    impact = _safe_string_series(data, "NEW_IMPACT").str.upper()
    pvs1_elig = _safe_string_series(data, "PVS1_ELIGIBILITY") if "PVS1_ELIGIBILITY" in data.columns else pd.Series("Allowed", index=data.index)
    symbol_series = _safe_string_series(data, "SYMBOL").str.upper()
    is_pvs1_vetoed = pvs1_elig.str.contains("Vetoed", case=False, na=False) | (symbol_series == "MYH7")

    is_lof = conseq.str.contains("stop_gained|frameshift|splice_donor|splice_acceptor")
    is_active_lof = is_lof & (~is_pvs1_vetoed)
    
    splice_max, splice_method = _extract_spliceai_details(data)
    data["SPLICE_MAX_UNIFIED"] = splice_max
    data["SPLICEAI_METHOD"] = splice_method

    spip_col_name = NEW_SPIP_COLUMNS.get("prediction", "SPiP") if "NEW_SPIP_COLUMNS" in globals() else "SPiP"
    spip_score = _safe_numeric_series(data, spip_col_name) if spip_col_name in data.columns else _safe_numeric_series(data, "SPiP")
    pangolin_score = _safe_numeric_series(data, "Pangolin_max_score")
    revel_score = _safe_numeric_series(data, "REVEL_score")
    am_score = _safe_numeric_series(data, "am_pathogenicity")
    pli_val = _safe_numeric_series(data, "pLI_gene_value")
    
    sv_aberrant = _safe_string_series(data, "SpliceVault_status") == "aberrant_event_detected"
    bp_disrupt = _safe_string_series(data, "Branchpoint_status") == "disrupted"

    utr_annot_cols = [
        "lost_stop_codon", "lost_start_codon", "mrl_gainedOrLost",
        "utr_num_kozak_gainedOrLost", "utr_num_uAUG_gainedOrLost",
        "num_polyA_signal_gainedOrLost", "tss_kozak_score_gainedOrLost"
    ]
    existing_utr_cols = [col for col in utr_annot_cols if col in data.columns]
    utr_disrupt = pd.Series(False, index=data.index)
    if existing_utr_cols:
        utr_values_of_interest = ["gained", "lost", "true"]
        utr_disrupt = (
            data[existing_utr_cols]
            .astype(str)
            .apply(lambda col: col.str.lower())
            .isin(utr_values_of_interest)
            .any(axis=1)
        )

    if "5UTR_consequence" in data.columns:
        utr_5prime_mask = data["5UTR_consequence"].notna() & (data["5UTR_consequence"].astype(str).str.strip() != "") & (data["5UTR_consequence"].astype(str) != "nan")
        utr_disrupt = utr_disrupt | utr_5prime_mask

    is_clinvar_path, is_clinvar_benign = _extract_clinvar_status(data)
    clinvar_stars = _safe_numeric_series(data, "ClinVar_Stars") if "ClinVar_Stars" in data.columns else pd.Series(0.0, index=data.index)

    score = np.zeros(len(data), dtype=float)
    score += np.where(is_active_lof, 35.0, 0.0)
    
    splicing_composite = np.maximum(splice_max, np.maximum(spip_score, pangolin_score))
    score += splicing_composite * 30.0
    score += np.where(sv_aberrant, 5.0, 0.0)
    score += np.where(bp_disrupt, 5.0, 0.0)
    
    # Phase 1: Calibrated In Silico Predictors (REVEL / AlphaMissense) - Inhibited if LoF
    insilico_points = np.zeros(len(data), dtype=float)
    for i in range(len(data)):
        if is_active_lof.iloc[i]:
            continue
        r_val = revel_score.iloc[i]
        am_val = am_score.iloc[i]
        if not np.isnan(r_val) and r_val > 0.0:
            if r_val >= 0.932: insilico_points[i] = 25.0
            elif r_val >= 0.773: insilico_points[i] = 18.0
            elif r_val >= 0.644: insilico_points[i] = 10.0
            elif r_val < 0.183: insilico_points[i] = -15.0
            elif r_val <= 0.289: insilico_points[i] = -8.0
        elif not np.isnan(am_val) and am_val > 0.0:
            if am_val >= 0.750: insilico_points[i] = 18.0
            elif am_val >= 0.564: insilico_points[i] = 10.0
            elif am_val < 0.340: insilico_points[i] = -8.0

    score += insilico_points

    score += np.where(utr_disrupt, 20.0, 0.0)
    score += np.where(pli_val >= 0.9, 10.0, 0.0)
    
    # Phase 1: ClinVar VCEP Star Ponderation
    vcep_path_mask = is_clinvar_path & (clinvar_stars >= 3.0)
    concordant_path_mask = is_clinvar_path & (~vcep_path_mask)
    score += np.where(vcep_path_mask, 50.0, np.where(concordant_path_mask, 45.0, 0.0))
    score -= np.where(is_clinvar_benign, 40.0, 0.0)
    
    # Phase 1: Disease-Specific FAF PopMax & MCAF Frequency Gating
    ba1_thresh = _safe_numeric_series(data, "FAF_POPMAX_BA1") if "FAF_POPMAX_BA1" in data.columns else pd.Series(0.002, index=data.index)
    bs1_thresh = _safe_numeric_series(data, "FAF_POPMAX_BS1") if "FAF_POPMAX_BS1" in data.columns else pd.Series(0.0001, index=data.index)
    mcaf_thresh = bs1_thresh / 2.5

    is_ba1 = af_series >= ba1_thresh
    is_bs1 = (~is_ba1) & (af_series >= bs1_thresh)
    is_pm2 = (~is_ba1) & (~is_bs1) & (af_series < mcaf_thresh)

    af_score = np.where(is_ba1, -100.0, np.where(is_bs1, -25.0, np.where(is_pm2, 15.0, 0.0)))
    score += af_score

    # Quality Score Penalty (-25.0 pts for LOW_QUAL VCF caller metrics)
    is_low_qual = (data["QC_STATUS"] == "LOW_QUAL") if "QC_STATUS" in data.columns else pd.Series(False, index=data.index)
    score -= np.where(is_low_qual, 25.0, 0.0)

    is_compound_het = (
        data["COMPOUND_HET_STATUS"].isin(["Phased TRANS Pair", "Unphased Candidate Pair"])
        if "COMPOUND_HET_STATUS" in data.columns else pd.Series(False, index=data.index)
    )

    if "INHERITANCE_MODEL" in data.columns:
        inh_series = data["INHERITANCE_MODEL"].astype(str)
        score += np.where(inh_series == "De Novo", 20.0, 0.0)
        score += np.where(
            (~is_compound_het) & inh_series.isin([
                "Autosomal Recessive (Hom)", "Autosomal Dominant / Heterozygous", "X-Linked"
            ]), 15.0, 0.0
        )
        score += np.where((~is_compound_het) & (inh_series == "Homozygous Variant"), 10.0, 0.0)

    if "COMPOUND_HET_STATUS" in data.columns:
        c_status_series = data["COMPOUND_HET_STATUS"].astype(str)
        score += np.where(c_status_series == "Phased TRANS Pair", 15.0, 0.0)
        score += np.where(c_status_series == "Unphased Candidate Pair", 10.0, 0.0)

    disease_ev_score = _safe_numeric_series(data, "DISEASE_EVIDENCE_SCORE")
    score += disease_ev_score

    hpo_match_mask = pd.Series(False, index=data.index)
    if hpo_terms and len(hpo_terms) > 0:
        target_genes = set()
        for term in hpo_terms:
            term_clean = term.split(" - ")[0].strip().upper()
            if not term_clean.startswith("HP:"):
                term_clean = f"HP:{term_clean}"
            if term_clean in HPO_GENE_PANELS:
                target_genes.update(HPO_GENE_PANELS[term_clean])
        if "SYMBOL" in data.columns and target_genes:
            hpo_match_mask = data["SYMBOL"].isin(target_genes)

    data["HPO_MATCH"] = hpo_match_mask
    score += np.where(hpo_match_mask, 15.0, 0.0)

    # ClinGen Tier Score Caps (Tier_C <= 45.0, Tier_D <= 35.0)
    gene_tier_series = _safe_string_series(data, "GENE_TIER") if "GENE_TIER" in data.columns else pd.Series("Tier_A", index=data.index)
    is_tier_c = gene_tier_series == "Tier_C"
    is_tier_d = gene_tier_series == "Tier_D"

    score = np.where(is_tier_c, np.minimum(score, 45.0), score)
    score = np.where(is_tier_d, np.minimum(score, 35.0), score)
    score = np.clip(score, 0.0, 100.0).round(1)
    data["VARIANT_PRIORITY_SCORE"] = score

    data = build_acmg_criteria(data)

    has_functional_signal = (
        (splice_max >= 0.20) |
        (spip_score >= 0.20) |
        (pangolin_score >= 0.20) |
        (sv_aberrant) |
        (bp_disrupt) |
        (utr_disrupt) |
        (revel_score >= 0.25) |
        (am_score >= 0.34) |
        (is_active_lof) |
        (is_clinvar_path) |
        (impact.isin(["MODERATE", "HIGH"])) |
        (conseq.str.contains("missense|splice|stop|frameshift|inframe|utr|5_prime|3_prime"))
    )

    tiers = pd.Series("Tier 4 (Benign / Tolerated)", index=data.index)
    has_splice_effect = (splice_max >= 0.20) | (spip_score >= 0.20) | (pangolin_score >= 0.20)

    tier3_mask = (~is_ba1) & (has_functional_signal) & ((score >= 25.0) | (has_splice_effect))
    tiers[tier3_mask] = "Tier 3 (VUS / Moderate Potential)"

    # LOW_QUAL, Tier_C, Tier_D, BA1, and BS1 variants are capped out of Tier 1 and Tier 2
    is_capped_out = is_low_qual | is_tier_c | is_tier_d | is_ba1 | is_bs1
    tier2_mask = (~is_capped_out) & (has_functional_signal) & (
        (score >= 50.0) |
        (splice_max >= 0.50) |
        (spip_score >= 0.30) |
        (pangolin_score >= 0.30) |
        (bp_disrupt) |
        (utr_disrupt) |
        (am_score >= 0.750) |
        (revel_score >= 0.773) |
        (is_active_lof)
    )
    tiers[tier2_mask] = "Tier 2 (Likely Deleterious / Strong Candidate)"

    tier1_mask = (~is_capped_out) & (is_pm2) & (
        (score >= 75.0) |
        ((is_active_lof) & (pli_val >= 0.9)) |
        (is_clinvar_path) |
        ((splice_max >= 0.80) & (sv_aberrant | is_active_lof | (pli_val >= 0.9))) |
        (revel_score >= 0.932) |
        ((am_score >= 0.750) & (revel_score >= 0.773))
    )
    tiers[tier1_mask] = "Tier 1 (Critical Pathogenic Candidate)"

    data["PRIORITY_TIER"] = pd.Categorical(
        tiers,
        categories=[
            "Tier 1 (Critical Pathogenic Candidate)",
            "Tier 2 (Likely Deleterious / Strong Candidate)",
            "Tier 3 (VUS / Moderate Potential)",
            "Tier 4 (Benign / Tolerated)"
        ],
        ordered=True
    )

    explanations = []
    rationales = []

    for i in range(len(data)):
        parts = []
        if is_lof.iloc[i]: parts.append("LoF (+35.0)")
        s_comp = splicing_composite.iloc[i]
        if s_comp > 0.05:
            m_label = splice_method.iloc[i]
            if "Custom" in m_label:
                parts.append(f"Splicing (+{s_comp*30.0:.1f}, SpliceAI {splice_max.iloc[i]:.2f} via Custom -D 10000)")
            elif "VEP" in m_label:
                parts.append(f"Splicing (+{s_comp*30.0:.1f}, SpliceAI {splice_max.iloc[i]:.2f} via VEP Standard ~50bp)")
            else:
                parts.append(f"Splicing (+{s_comp*30.0:.1f})")
        if sv_aberrant.iloc[i]: parts.append("SpliceVault (+5.0)")
        if bp_disrupt.iloc[i]: parts.append("Branchpoint (+5.0)")

        # Phase 1 In Silico Explanation
        in_pts = insilico_points[i]
        if in_pts != 0.0:
            if in_pts > 0: parts.append(f"InSilico PP3 (+{in_pts:.1f})")
            else: parts.append(f"InSilico BP4 ({in_pts:.1f})")

        if utr_disrupt.iloc[i]: parts.append("UTR (+20.0)")
        if pli_val.iloc[i] >= 0.9: parts.append("LoF-Intolerant (+10.0)")
        
        if vcep_path_mask.iloc[i]: parts.append("ClinVar VCEP Pathogenic (+50.0)")
        elif concordant_path_mask.iloc[i]: parts.append("ClinVar Pathogenic (+45.0)")
        elif is_clinvar_benign.iloc[i]: parts.append("ClinVar Benign (-40.0)")

        af_v = af_series.iloc[i]
        if is_ba1.iloc[i]: parts.append(f"BA1 Exceeds FAF_PopMax Threshold (Terminal Tier 4)")
        elif is_bs1.iloc[i]: parts.append(f"BS1_Strong AF > FAF_BS1 (-25.0)")
        elif is_pm2.iloc[i]: parts.append(f"PM2_Supporting AF < MCAF (+15.0)")

        c_status = str(data["COMPOUND_HET_STATUS"].iloc[i]) if "COMPOUND_HET_STATUS" in data.columns else "No Pair"
        if c_status in ["Phased TRANS Pair", "Unphased Candidate Pair"]:
            if c_status == "Phased TRANS Pair":
                parts.append("Compound Het TRANS Phased (+15.0)")
            elif c_status == "Unphased Candidate Pair":
                parts.append("Compound Het Unphased Candidate (+10.0)")
        elif "INHERITANCE_MODEL" in data.columns:
            inh_val = str(data["INHERITANCE_MODEL"].iloc[i])
            if inh_val == "De Novo":
                parts.append("De Novo (+20.0)")
            elif inh_val in ["Autosomal Recessive (Hom)", "Autosomal Dominant / Heterozygous", "X-Linked"]:
                parts.append(f"{inh_val} (+15.0)")
            elif inh_val == "Homozygous Variant":
                parts.append("Homozygous Variant (+10.0)")

        if is_low_qual.iloc[i]:
            qc_det = str(data["QC_DETAIL"].iloc[i]) if "QC_DETAIL" in data.columns else "Low QUAL/GQ/FILTER"
            parts.append(f"Low Quality Caller Metric (-25.0, {qc_det})")

        exp_str = " | ".join(parts) if parts else "Neutral Rarity / Benign"
        explanations.append(exp_str)

        t_str = str(tiers.iloc[i])
        final_score = data["VARIANT_PRIORITY_SCORE"].iloc[i]
        if "Tier 1" in t_str:
            if vcep_path_mask.iloc[i]:
                rat = "ClinGen VCEP Pathogenic / Likely Pathogenic Expert Entry"
            elif is_clinvar_path.iloc[i]:
                rat = "ClinVar Pathogenic / Likely Pathogenic Entry"
            elif is_lof.iloc[i] and pli_val.iloc[i] >= 0.9:
                rat = f"High-Confidence Loss-of-Function in Constrained Gene (pLI = {pli_val.iloc[i]:.2f} >= 0.90)"
            elif splice_max.iloc[i] >= 0.80:
                rat = f"Ultra-High Splicing Alteration (SpliceAI = {splice_max.iloc[i]:.2f} >= 0.80 via {splice_method.iloc[i]})"
            elif revel_score.iloc[i] >= 0.932:
                rat = f"Strong Pathogenic Missense (REVEL = {revel_score.iloc[i]:.3f} >= 0.932, PP3_Strong)"
            else:
                rat = f"High-Priority Composite Score ({final_score:.1f}/100 >= 75.0) & Ultra-Rare Population Frequency (AF < MCAF)"
        elif "Tier 2" in t_str:
            if is_lof.iloc[i]:
                rat = f"Loss-of-Function Alteration in Rare Variant (AF = {af_v:.4f})"
            elif splice_max.iloc[i] >= 0.50:
                rat = f"High Splicing Delta Score (SpliceAI = {splice_max.iloc[i]:.2f} >= 0.50 via {splice_method.iloc[i]})"
            elif revel_score.iloc[i] >= 0.773:
                rat = f"Moderate Pathogenic Missense (REVEL = {revel_score.iloc[i]:.3f} >= 0.773, PP3_Moderate)"
            elif am_score.iloc[i] >= 0.750:
                rat = f"Moderate Pathogenic Missense (AlphaMissense = {am_score.iloc[i]:.3f} >= 0.750, PP3_Moderate)"
            else:
                rat = f"Strong Priority Composite Score ({final_score:.1f}/100 >= 50.0) & Rare Population Frequency"
        elif "Tier 3" in t_str:
            if has_splice_effect.iloc[i]:
                rat = f"Moderate Splicing Delta Alteration (SpliceAI = {splice_max.iloc[i]:.2f} >= 0.20 via {splice_method.iloc[i]})"
            else:
                rat = f"Moderate Priority Composite Score ({final_score:.1f}/100 >= 25.0) & Rare Population Frequency"
        else:
            if is_ba1.iloc[i]:
                rat = f"Stand-Alone Benign: AF ({af_v:.4f}) Exceeds FAF_PopMax Threshold ({ba1_thresh.iloc[i]:.4f})"
            elif is_bs1.iloc[i]:
                rat = f"Strong Benign: AF ({af_v:.4f}) Exceeds FAF_BS1 Threshold ({bs1_thresh.iloc[i]:.4f})"
            elif is_clinvar_benign.iloc[i]:
                rat = "ClinVar Benign / Likely Benign Entry"
            else:
                rat = f"Low Priority Composite Score ({final_score:.1f}/100 < 25.0) Without Functional Evidence"

        if is_low_qual.iloc[i]:
            qc_det = str(data["QC_DETAIL"].iloc[i]) if "QC_DETAIL" in data.columns else "Low QUAL/GQ/FILTER"
            rat += f" [⚠️ LOW_QUAL Caller Warning: {qc_det}]"

        rationales.append(rat)

    data["SCORE_EXPLANATION"] = explanations
    data["TIER_RATIONALE"] = rationales

    data = order_clinical_columns(data)
    logger.info(f"PRIORITY_TIER assigned with counts: {dict(data['PRIORITY_TIER'].value_counts())}")
    return data


def order_clinical_columns(data: pd.DataFrame) -> pd.DataFrame:
    """Re-orders DataFrame columns so primary clinical triage tags are at the front."""
    front_cols = [
        "SYMBOL", "Locus", "HGVSc", "HGVSp", "Consequence", "NEW_IMPACT", "IMPACT",
        "PRIORITY_TIER", "VARIANT_PRIORITY_SCORE", "TIER_RATIONALE", "SCORE_EXPLANATION",
        "ACMG_CRITERIA", "GENE_TIER", "CLINGEN_VALIDITY", "GAMMA_WEIGHT", "PATHOGENIC_MECHANISM", "PVS1_ELIGIBILITY",
        "INHERITANCE_MODEL", "COMPOUND_HET_STATUS", "COMPOUND_HET_PAIR",
        "DISEASE_GENE_EVIDENCE", "DISEASE_CURATION_SOURCE", "HIGH_RISK_DISEASE_GENE",
        "ACMG_SF_V3_2", "HPO_MATCH", "CRITICAL_DOMAINS_HOTSPOTS"
    ]
    
    evidence_cols = [
        "CLINVAR_DISPLAY", "clinvar_clnsig", "CLNSIG",
        "gnomADv4_AF_grpmax_joint", "MAX_AF", "af", "gnomad_af", "AF",
        "pLI_gene_value", "pLI",
        "spliceai_custom_MAX", "spliceAI_MAX", "SPLICE_MAX_UNIFIED", "SPLICEAI_METHOD",
        "SPiP_prediction", "SPiP", "Pangolin_max_score", "Pangolin",
        "SpliceVault_status", "Branchpoint_status",
        "am_pathogenicity", "AlphaMissense", "REVEL_score", "REVEL",
        "lost_stop_codon", "lost_start_codon", "mrl_gainedOrLost", "utr_num_kozak_gainedOrLost", "5UTR_consequence",
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
        "gnomAD Allele Frequency": ["gnomADv4_AF_grpmax_joint", "MAX_AF", "af", "gnomad_af", "AF"],
        "SpliceAI Custom Window (-D 10000 / 20kb)": ["spliceai_custom_MAX", "SpliceAI_custom_MAX", "spliceai_custom_DS_AG", "SpliceAI_custom_DS_AG", "SpliceAI_D10000", "SpliceAI_10k"],
        "SpliceAI VEP Standard (~50bp)": ["spliceAI_MAX", "SpliceAI_pred_DS_AG", "SpliceAI_DS_AG", "SpliceAI"],
        "SPiP Deep Splicing": ["SPiP_prediction", "SPiP"],
        "Pangolin Non-Coding": ["Pangolin_max_score", "Pangolin"],
        "AlphaMissense": ["am_pathogenicity", "AlphaMissense"],
        "REVEL Missense Ensemble": ["REVEL_score", "REVEL"],
        "SpliceVault RNA-Seq Cohort": ["SpliceVault_status", "SpliceVault"],
        "Branchpoint Disruption": ["Branchpoint_status", "Branchpoint"],
        "ClinVar Pathogenicity": ["clinvar_clnsig", "CLINVAR_DISPLAY", "CLNSIG", "ClinVar"],
        "pLI Gene LoF Intolerance": ["pLI_gene_value", "pLI"],
        "5'/3' UTR & Kozak Regulatory": ["lost_stop_codon", "5UTR_consequence", "utr_num_kozak_gainedOrLost", "mrl_gainedOrLost"]
    }

    audit_results = {}
    cols_lower = {c.lower(): c for c in data.columns}

    for tool_name, candidates in tool_map.items():
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
    if data[freq_column].isnull().sum() == data.shape[0]:
        logger.debug(f"Values for {freq_column}: {data[freq_column].unique()}")
        logger.error(
            f"All values in {freq_column} are NaN, please check the column"
        )
        raise ValueError(f"All values in {freq_column} are NaN, please check the column")
    try:
        data[freq_column] = data[freq_column].replace(".", np.nan)
    except KeyError as e:
        logger.error(f"Column {freq_column} not found in data")
        freq_similar = [
            col for col in data.columns if freq_column.lower() in col.lower()
        ]
        logger.warning(f"Maybe you meant {freq_similar}")
        raise e
    data[freq_column] = data[freq_column].fillna(0)
    data[freq_column] = data[freq_column].astype(float)
    logger.info(f"NA values filled with 0 in column {freq_column}")
    new_df = data[data[freq_column] < freq]
    logger.info(f"Dataframe shape after frequency [{freq}] filtering {new_df.shape[0]}")
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

    def evaluate_variant(row):
        prio = row["gene_priority"]
        impact = row["NEW_IMPACT"]
        freq_val = row[freq_col]
        limit = thresholds.get((prio, impact))
        if limit and freq_val < limit:
            return "Candidate"
        return "Excluded"

    data["Variant_rank"] = data.apply(evaluate_variant, axis=1)
    logger.info("New column `Variant_rank` created based on gene priority, impact and frequency")
    return data
