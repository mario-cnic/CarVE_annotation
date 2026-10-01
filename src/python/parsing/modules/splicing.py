"""
Deep Splicing Predictors, Intron Offset Calculation, and RNA-Seq Cohort Annotations.
"""
import re
import numpy as np
import pandas as pd
from . import logger
from .io_qc import _safe_numeric_series

SPLICING_COLUMNS = [
    "spliceai_custom_DS_AL",
    "spliceai_custom_DS_AG",
    "spliceai_custom_DS_DL",
    "spliceai_custom_DS_DG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_DS_DG",
]

NEW_SPIP_COLUMNS = {
    "interpretation": "SPiP_interpretation",
    "ci": "SPiP_ci",
    "prediction": "SPiP_prediction",
    "confidence": "SPiP_confidence",
    "score": "SPiP_score",
    "mechanism": "SPiP_mechanism",
    "mut_in_pb": "SPiP_mutInPBarea",
    "delta_esr": "SPiP_deltaESRscore",
    "proba_crypt": "SPiP_probaCryptMut",
    "status": "SPiP_status"
}

NEW_PANGOLIN_COLUMNS = {
    "max_score": "Pangolin_max_score",
    "heart_lv": "Pangolin_heart_lv_score",
    "heart_aa": "Pangolin_heart_aa_score",
    "status": "Pangolin_status"
}

NEW_BRANCHPOINT_COLUMNS = {
    "branchpointer_prob": "Branchpointer_prob",
    "branchpointer_energy": "Branchpointer_U2_energy",
    "branchpoint_disrupted": "Branchpoint_disrupted",
    "branchpoint_status": "Branchpoint_status",
    "labranchor_score": "LaBranchoR_score",
    "labranchor_acc_dist": "LaBranchoR_acc_dist",
    "labranchor_status": "LaBranchoR_status",
}

# Informational only (see walkthrough/20260910_predictor_inventory_and_audit.md, PRED-11):
# MaxEntScan_ref/_alt/_diff are VEP-computed raw pass-through values. This module derives a
# normalized decrease and a disruption flag from them, but neither currently feeds
# build_newImpact/build_priority_tier -- by explicit decision, MaxEntScan does not yet move
# NEW_IMPACT or PRIORITY_TIER, unlike Pangolin/SpliceAI/SPiP/SpliceVarDB.
NEW_MAXENTSCAN_COLUMNS = {
    "pct_decrease": "MaxEntScan_pct_decrease",
    "disrupted": "MaxEntScan_disrupted",
    "status": "MaxEntScan_status",
}
MAXENTSCAN_DISRUPTION_THRESHOLD = 0.15  # >=15% relative decrease from ref to alt


def parse_spip(data: pd.DataFrame, spip_col: str = "SPiP"):
    """
    For each row, extract SPiP annotations matching the target SYMBOL or transcript.
    """
    def extract_spip_fields(row):
        spip_raw = str(row.get(spip_col, "."))
        if pd.isna(spip_raw) or spip_raw in (".", "nan", "N/A", ""):
            return pd.Series([np.nan, np.nan, np.nan, np.nan, np.nan, "none", np.nan, np.nan, np.nan, "not_covered"])
        
        if "caused an error in SPiP execution" in spip_raw:
            return pd.Series(["Error", np.nan, np.nan, np.nan, np.nan, "error", np.nan, np.nan, np.nan, "error"])

        target_symbol = str(row.get("SYMBOL", ""))
        spip_entries = re.split(r'[,&]', spip_raw)
        matched_fields = None

        for entry in spip_entries:
            fields = entry.split("|")
            if len(fields) > 12 and (not target_symbol or fields[12] == target_symbol):
                matched_fields = fields
                break
        
        if matched_fields is None and len(spip_entries) > 0:
            matched_fields = spip_entries[0].split("|")

        if not matched_fields or len(matched_fields) < 5:
            return pd.Series([np.nan, np.nan, np.nan, np.nan, np.nan, "none", np.nan, np.nan, np.nan, "not_covered"])

        interpretation = matched_fields[2] if len(matched_fields) > 2 else np.nan
        ci = matched_fields[3] if len(matched_fields) > 3 else np.nan
        
        try:
            prediction = float(matched_fields[4]) if len(matched_fields) > 4 and matched_fields[4] != "." else np.nan
        except ValueError:
            prediction = np.nan

        score = prediction
        confidence = ci
        mut_in_pb = matched_fields[20] if len(matched_fields) > 20 else np.nan
        
        try:
            delta_esr = float(matched_fields[21]) if len(matched_fields) > 21 and matched_fields[21] != "." else np.nan
        except ValueError:
            delta_esr = np.nan

        try:
            proba_crypt = float(matched_fields[24]) if len(matched_fields) > 24 and matched_fields[24] != "." else np.nan
        except ValueError:
            proba_crypt = np.nan

        class_crypt = matched_fields[25] if len(matched_fields) > 25 else ""

        mechanism = "none"
        if not pd.isna(proba_crypt) and proba_crypt >= 0.20 or class_crypt in ("Cryptic", "High Risk", "Medium Risk"):
            mechanism = "cryptic_site"
        elif str(mut_in_pb).lower() in ("1", "true", "yes"):
            mechanism = "bp_disruption"
        elif not pd.isna(delta_esr) and abs(delta_esr) > 0.05:
            mechanism = "esr_alteration"
        elif not pd.isna(prediction) and prediction >= 0.20:
            mechanism = "complex_splicing"

        status = "scored" if not pd.isna(prediction) else "not_covered"
        return pd.Series([interpretation, ci, prediction, confidence, score, mechanism, mut_in_pb, delta_esr, proba_crypt, status])

    logger.info("Building expanded SPiP columns with mechanism classification")
    if spip_col not in data.columns:
        logger.warning(f"Column {spip_col} not found in data")
        return data

    try:
        spip_extracted = data.apply(extract_spip_fields, axis=1)
        spip_extracted.columns = list(NEW_SPIP_COLUMNS.values())
        for col_name in spip_extracted.columns:
            data[col_name] = spip_extracted[col_name]
    except Exception as e:
        logger.error(f"Error extracting SPiP fields: {e}")
        if data.shape[0] == 0:
            logger.error("Dataframe is empty, cannot extract SPiP fields")
            raise ValueError("Dataframe is empty, cannot extract SPiP fields")
        else:
            raise e

    logger.info("SPiP columns extracted, categorized, and status contract applied")
    return data


def decode_splicevault_prediction(pred_str: str, var_pos: int, strand: str = '+') -> str:
    """Parses colon-separated values into human-readable details."""
    if not pred_str or pd.isna(pred_str) or pred_str in ("N/A", "nan", "."):
        return "No prediction"
    
    events = []
    for event in re.split(r'[,&]', str(pred_str)):
        event = event.strip()
        if not event or event == ".":
            continue
        parts = event.split(':')
        if len(parts) < 3:
            events.append(event)
            continue
        
        rank = parts[0]
        ev_type = parts[1]
        impact = parts[2]
        pct = parts[3] if len(parts) > 3 else None
        frame = parts[4] if len(parts) > 4 else None
        
        type_map = {
            "CD": "Cryptic Donor",
            "CA": "Cryptic Acceptor",
            "ES": "Exon Skipping",
            "IR": "Intron Retention"
        }
        full_type = type_map.get(ev_type, ev_type)
        
        if ev_type == "ES":
            impact_desc = f"Exon {impact} skipping"
        else:
            try:
                offset = int(impact)
                predicted_site = var_pos - offset if strand == '-' else var_pos + offset
                direction = "downstream" if offset > 0 else "upstream"
                abs_offset = abs(offset)
                rel_desc = "at mutation site" if offset == 0 else f"{abs_offset} bp {direction} (pos:{predicted_site})"
                impact_desc = f"Cryptic Site at {rel_desc}"
            except ValueError:
                impact_desc = f"Offset {impact} bp"
        
        ev_desc = f"Rank {rank}: {full_type} ({impact_desc})"
        if pct:
            ev_desc += f", {pct}% support"
        if frame:
            ev_desc += f", {frame.replace('_', ' ')}"
        events.append(ev_desc)
        
    return "; ".join(events) if events else "No prediction"


def parse_splicevault(data: pd.DataFrame):
    """Decodes SpliceVault predictions and assigns status contract semantics."""
    if "SpliceVault_top_events" not in data.columns:
        logger.info("SpliceVault_top_events not found in data, skipping SpliceVault decoding")
        return data

    logger.info("Parsing SpliceVault predictions and building status contract...")
    
    def process_row(row):
        pred_str = row.get("SpliceVault_top_events")
        site_cnt = row.get("SpliceVault_site_sample_count")
        site_pos = row.get("SpliceVault_site_pos")
        
        has_pred = not pd.isna(pred_str) and str(pred_str) not in ("", ".", "nan", "N/A", "None")
        has_coverage = not pd.isna(site_cnt) and str(site_cnt) not in ("", ".", "nan", "N/A", "None") or \
                       not pd.isna(site_pos) and str(site_pos) not in ("", ".", "nan", "N/A", "None")
        
        if not has_pred:
            if has_coverage:
                return pd.Series(["No aberrant events detected in RNA-seq cohort", "no_events_found"])
            else:
                return pd.Series(["Not covered", "not_covered"])
        
        pos_val = row.get("POS")
        try:
            var_pos = int(pos_val)
        except (ValueError, TypeError):
            var_pos = 0
            
        strand_val = row.get("STRAND")
        if pd.isna(strand_val) or str(strand_val) in ("", ".", "N/A"):
            strand_val = row.get("strand", "+")
            
        strand_str = str(strand_val).strip()
        strand = '-' if strand_str in ("-1", "-", "negative") else '+'
        
        decoded = decode_splicevault_prediction(str(pred_str), var_pos, strand)
        if decoded and decoded != "No prediction":
            status = "aberrant_event_detected"
        elif has_coverage:
            decoded = "No aberrant events detected in RNA-seq cohort"
            status = "no_events_found"
        else:
            decoded = "Not covered"
            status = "not_covered"
            
        return pd.Series([decoded, status])

    try:
        res = data.apply(process_row, axis=1)
        res.columns = ["SpliceVault_Predictions_Decoded", "SpliceVault_status"]
        data["SpliceVault_Predictions_Decoded"] = res["SpliceVault_Predictions_Decoded"]
        data["SpliceVault_status"] = res["SpliceVault_status"]
        logger.info("SpliceVault predictions decoded and status contract applied successfully")
    except Exception as e:
        logger.error(f"Error parsing SpliceVault predictions: {e}")

    return data


def parse_pangolin(data: pd.DataFrame, pangolin_col: str = "Pangolin"):
    """Parses Pangolin splicing scores."""
    if pangolin_col not in data.columns and "Pangolin_max_score" not in data.columns:
        logger.info("Pangolin column not found in data, skipping Pangolin parsing")
        return data

    logger.info("Parsing Pangolin predictions and status contract...")

    def process_row(row):
        pang_raw = str(row.get(pangolin_col, "."))
        if pd.isna(pang_raw) or pang_raw in (".", "nan", "N/A", ""):
            max_s = row.get("Pangolin_max_score")
            if not pd.isna(max_s) and str(max_s) not in (".", "nan"):
                try:
                    val = float(max_s)
                    return pd.Series([val, float(row.get("Pangolin_heart_lv_score", val)), float(row.get("Pangolin_heart_aa_score", val)), "scored"])
                except ValueError:
                    pass
            return pd.Series([np.nan, np.nan, np.nan, "not_covered"])

        parts = pang_raw.split("|")
        scores = []
        for p in parts:
            if ":" in p:
                sub_parts = p.split(":")
                if len(sub_parts) == 2:
                    try:
                        val = abs(float(sub_parts[1]))
                        if 0.0 <= val <= 1.0:
                            scores.append(val)
                    except ValueError:
                        pass

        if not scores:
            return pd.Series([np.nan, np.nan, np.nan, "not_covered"])

        max_score = max(scores)
        return pd.Series([max_score, max_score, max_score, "scored"])

    try:
        res = data.apply(process_row, axis=1)
        res.columns = list(NEW_PANGOLIN_COLUMNS.values())
        for col_name in res.columns:
            data[col_name] = res[col_name]
        logger.info("Pangolin columns parsed and status contract applied successfully")
    except Exception as e:
        logger.error(f"Error parsing Pangolin: {e}")

    return data


def parse_branchpointer(data: pd.DataFrame):
    """Parses Branchpointer and LaBranchoR predicted scores."""
    check_cols = ["Branchpointer_prob", "LaBranchoR_score", "Branchpoint_disrupted"]
    if not any(c in data.columns for c in check_cols):
        logger.info("Branchpoint columns not found in data, skipping Branchpointer/LaBranchoR parsing")
        return data

    logger.info("Parsing Branchpoint predictions and status contract...")

    def process_row(row):
        bp_prob = row.get("Branchpointer_prob")
        u2_energy = row.get("Branchpointer_U2_energy")
        disrupted = row.get("Branchpoint_disrupted", "NO")
        labranchor = row.get("LaBranchoR_score")
        acc_dist = row.get("LaBranchoR_acc_dist")

        has_bp = not pd.isna(bp_prob) and str(bp_prob) not in (".", "nan", "")
        has_lab = not pd.isna(labranchor) and str(labranchor) not in (".", "nan", "")

        bp_prob_val = float(bp_prob) if has_bp else np.nan
        u2_energy_val = float(u2_energy) if (not pd.isna(u2_energy) and str(u2_energy) not in (".", "nan", "")) else np.nan
        disrupted_val = str(disrupted) if str(disrupted) in ("YES", "NO") else ("YES" if (has_bp and bp_prob_val >= 0.5) else "NO")
        bp_status = "scored" if (has_bp or disrupted_val == "YES") else "not_covered"

        lab_score_val = float(labranchor) if has_lab else np.nan
        acc_dist_val = int(acc_dist) if (not pd.isna(acc_dist) and str(acc_dist) not in (".", "nan", "")) else np.nan
        lab_status = "scored" if has_lab else "not_covered"

        return pd.Series([bp_prob_val, u2_energy_val, disrupted_val, bp_status, lab_score_val, acc_dist_val, lab_status])

    try:
        res = data.apply(process_row, axis=1)
        res.columns = list(NEW_BRANCHPOINT_COLUMNS.values())
        for col_name in res.columns:
            data[col_name] = res[col_name]
        logger.info("Branchpoint columns parsed and status contract applied successfully")
    except Exception as e:
        logger.error(f"Error parsing Branchpoint columns: {e}")

    return data


def parse_maxentscan(data: pd.DataFrame):
    """Parses VEP's MaxEntScan_ref/_alt scores into a normalized decrease and a disruption
    flag, plus a scored/not_covered status contract (a non-positive reference score is a
    real, meaningful outcome -- not a predictor error).

    Informational only: does not currently feed build_newImpact/build_priority_tier
    (see NEW_MAXENTSCAN_COLUMNS comment above). MAXENTSCAN_DISRUPTION_THRESHOLD (0.15)
    flags a weakened splice site only (ref > alt); a newly created/strengthened cryptic
    site (alt > ref) is not flagged by this metric.
    """
    check_cols = ["MaxEntScan_ref", "MaxEntScan_alt"]
    if not all(c in data.columns for c in check_cols):
        logger.info("MaxEntScan columns not found in data, skipping MaxEntScan parsing")
        return data

    logger.info("Parsing MaxEntScan scores and status contract...")

    def _to_float(val):
        if pd.isna(val) or str(val).strip() in (".", "nan", ""):
            return np.nan
        try:
            return float(val)
        except (TypeError, ValueError):
            return np.nan

    def process_row(row):
        ref_val = _to_float(row.get("MaxEntScan_ref"))
        alt_val = _to_float(row.get("MaxEntScan_alt"))

        if pd.isna(ref_val) or pd.isna(alt_val):
            return pd.Series([np.nan, "NO", "not_covered"])

        if ref_val <= 0:
            # A non-positive reference score is a normal, meaningful result (the reference
            # sequence isn't recognized as a splice-site consensus at all -- exactly the
            # situation for a variant creating a cryptic site), not a predictor failure.
            # The relative-decrease metric is undefined here, so this is not_covered, not error.
            return pd.Series([np.nan, "NO", "not_covered"])

        pct_decrease = (ref_val - alt_val) / ref_val
        disrupted = "YES" if pct_decrease >= MAXENTSCAN_DISRUPTION_THRESHOLD else "NO"
        return pd.Series([pct_decrease, disrupted, "scored"])

    try:
        res = data.apply(process_row, axis=1)
        res.columns = list(NEW_MAXENTSCAN_COLUMNS.values())
        for col_name in res.columns:
            data[col_name] = res[col_name]
        logger.info("MaxEntScan columns parsed and status contract applied successfully")
    except Exception as e:
        logger.error(f"Error parsing MaxEntScan columns: {e}")

    return data


def build_signed_intron_offset(data: pd.DataFrame):
    """Creates `intron_offset_signed` (Integer) and `splice_side` ('donor', 'acceptor', 'exonic')."""
    if "intron_offset_signed" in data.columns and "splice_side" in data.columns:
        logger.info("Columns `intron_offset_signed` and `splice_side` already exist")
        return data

    logger.info("Building signed intron offset and splice_side columns...")
    
    def process_row(row):
        if "INTRON_OFFSET_SIGNED" in row and not pd.isna(row["INTRON_OFFSET_SIGNED"]):
            try:
                val = int(row["INTRON_OFFSET_SIGNED"])
                side = str(row.get("SPLICE_SIDE", "donor" if val > 0 else ("acceptor" if val < 0 else "exonic")))
                return pd.Series([val, side])
            except ValueError:
                pass

        hgvsc = str(row.get("HGVSc", ""))
        if pd.isna(hgvsc) or hgvsc in ("", ".", "nan", "N/A"):
            hgvsc = str(row.get("CDNA_NAME", ""))
            
        match = re.search(r'c\.\d+([+-])(\d+)', hgvsc)
        if match:
            sign, val_str = match.groups()
            val = int(val_str)
            if sign == '+':
                return pd.Series([val, "donor"])
            elif sign == '-':
                return pd.Series([-val, "acceptor"])
                
        offset_val = row.get("HGVS_OFFSET")
        if not pd.isna(offset_val) and str(offset_val) not in ("", ".", "nan"):
            try:
                off_int = int(float(offset_val))
                if "-" in hgvsc:
                    return pd.Series([-abs(off_int), "acceptor"])
                elif "+" in hgvsc:
                    return pd.Series([abs(off_int), "donor"])
                return pd.Series([off_int, "exonic" if off_int == 0 else "donor"])
            except ValueError:
                pass
                
        return pd.Series([0, "exonic"])

    res = data.apply(process_row, axis=1)
    res.columns = ["intron_offset_signed", "splice_side"]
    data["intron_offset_signed"] = res["intron_offset_signed"]
    data["splice_side"] = res["splice_side"]
    return data


def parse_spliceai_custom(data: pd.DataFrame) -> pd.DataFrame:
    """Unpack custom SpliceAI (-D 10000 / 20kb window) prediction string into dedicated columns."""
    col = "SpliceAI"
    if col not in data.columns:
        logger.info(f"Column `{col}` not found in data, skipping custom SpliceAI unpacking")
        return data

    logger.info("Unpacking custom SpliceAI (-D 10000) predictions into dedicated `spliceai_custom_*` columns...")

    def extract_fields(val):
        if pd.isna(val) or val is None or str(val).strip() in ("", ".", "nan", "None", "-"):
            return [np.nan] * 10
        first_entry = str(val).split(",")[0].strip()
        parts = first_entry.split("|")
        if len(parts) >= 10:
            try:
                symbol = parts[1] if parts[1] != "." else np.nan
                ds_ag = float(parts[2]) if parts[2] != "." else np.nan
                ds_al = float(parts[3]) if parts[3] != "." else np.nan
                ds_dg = float(parts[4]) if parts[4] != "." else np.nan
                ds_dl = float(parts[5]) if parts[5] != "." else np.nan
                dp_ag = float(parts[6]) if parts[6] != "." else np.nan
                dp_al = float(parts[7]) if parts[7] != "." else np.nan
                dp_dg = float(parts[8]) if parts[8] != "." else np.nan
                dp_dl = float(parts[9]) if parts[9] != "." else np.nan
                scores = [x for x in [ds_ag, ds_al, ds_dg, ds_dl] if not np.isnan(x)]
                max_score = max(scores) if scores else np.nan
                return [symbol, ds_ag, ds_al, ds_dg, ds_dl, dp_ag, dp_al, dp_dg, dp_dl, max_score]
            except Exception:
                return [np.nan] * 10
        return [np.nan] * 10

    extracted = pd.DataFrame(
        data[col].apply(extract_fields).tolist(),
        index=data.index,
        columns=[
            "spliceai_custom_SYMBOL",
            "spliceai_custom_DS_AG",
            "spliceai_custom_DS_AL",
            "spliceai_custom_DS_DG",
            "spliceai_custom_DS_DL",
            "spliceai_custom_DP_AG",
            "spliceai_custom_DP_AL",
            "spliceai_custom_DP_DG",
            "spliceai_custom_DP_DL",
            "spliceai_custom_MAX",
        ],
    )

    for c in extracted.columns:
        data[c] = extracted[c]

    logger.info("Custom SpliceAI unpacking completed successfully")
    return data


def build_spliceMAX(
    splicemax: pd.DataFrame,
    splice_cols: list = [
        "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_DS_DL",
        "SpliceAI_pred_DS_DG",
        "SpliceAI_DS_AL",
        "SpliceAI_DS_AG",
        "SpliceAI_DS_DL",
        "SpliceAI_DS_DG",
    ],
):
    avail_cols = [x for x in splice_cols if x in splicemax.columns]
    if avail_cols:
        splicemax["spliceAI_MAX"] = splicemax[avail_cols].max(axis=1)
        splicemax["spliceAI_MAX"] = splicemax["spliceAI_MAX"].astype(float)
        splicemax["SpliceAI_status"] = np.where(splicemax["spliceAI_MAX"].isna(), "not_covered", "scored")
        splicemax["spliceAI_MAX"] = splicemax["spliceAI_MAX"].fillna(0)
        logger.info(f"New column `spliceAI_MAX` and `SpliceAI_status` created using columns {avail_cols}")
    else:
        logger.info("skipping spliceAI max calculation...")
        splicemax["SpliceAI_status"] = "not_covered"

    return splicemax


def _extract_spliceai_details(df: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    splice_custom = _safe_numeric_series(df, "spliceai_custom_MAX")
    splice_vep = _safe_numeric_series(df, "spliceAI_MAX")
    final_scores = np.where(splice_custom > 0, splice_custom, splice_vep)
    
    methods = []
    has_custom = "spliceai_custom_MAX" in df.columns
    has_vep = "spliceAI_MAX" in df.columns or "SpliceAI_pred_DS_AG" in df.columns or "SpliceAI_DS_AG" in df.columns

    for idx in range(len(df)):
        c_s = splice_custom.iloc[idx]
        v_s = splice_vep.iloc[idx]
        if c_s > 0:
            methods.append("Custom Window (-D 10000 / 20kb)")
        elif v_s > 0 or has_vep:
            methods.append("VEP Standard (~50bp)")
        elif has_custom:
            methods.append("Custom Window (-D 10000 / 20kb)")
        else:
            methods.append("Not Covered")

    return pd.Series(final_scores, index=df.index), pd.Series(methods, index=df.index)


def _extract_spliceai_score(df: pd.DataFrame) -> pd.Series:
    scores, _ = _extract_spliceai_details(df)
    return scores
