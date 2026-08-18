#!/usr/bin/env python3
"""
Medical & Geneticist-Grade 4-Tab Clinical Prioritization Report Generator

Generates a standalone, interactive HTML dashboard for a specific gene with 4 tabs:
  1. 📊 General / Overview (Executive summary, global triage KPIs, top candidates)
  2. 🧬 Splicing Tab (Dedicated splicing metrics, custom 20kb delta channels, SpliceVault, Branchpoint)
  3. 🎯 UTR Tab (5' and 3' UTR variants, uAUG/Kozak alterations, translation consequences)
  4. 🔬 Missense Tab (AlphaMissense, REVEL, CADD, consensus pathogenicity)

Features:
  - Dynamically generated <thead> ensuring exact column header alignment.
  - Multi-candidate ClinVar extraction (clinvar_clnsig, CLIN_SIG, CLNSIG).
  - Dedicated gnomADv4 AF grpmax with fallback and explicit warning banner.
  - SpliceAI 20kb context window & tool availability disclaimers.
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["ARROW_IO_THREADS"] = "1"

import sys
import glob
import argparse
import logging
import pandas as pd
import numpy as np
from datetime import datetime

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("clinical_report_generator")

def _safe_str(val, max_len=40):
    if pd.isna(val) or val is None or str(val).strip() == "" or str(val).lower() == "nan" or str(val).lower() == "none":
        return "-"
    s = str(val).strip()
    if len(s) > max_len:
        return s[:max_len-3] + "..."
    return s

def _safe_float(val, digits=3):
    try:
        f = float(val)
        if np.isnan(f):
            return "-"
        return f"{f:.{digits}f}"
    except (ValueError, TypeError):
        return "-"

def get_col_numeric(df, col, default=0.0):
    if col in df.columns:
        return pd.to_numeric(df[col], errors="coerce").fillna(default)
    return pd.Series(default, index=df.index)

def get_col_string(df, col, default="-"):
    if col in df.columns:
        return df[col].astype(str)
    return pd.Series(default, index=df.index)

def _extract_clinvar_display(df):
    cln_candidates = ["clinvar_clnsig", "CLIN_SIG", "CLNSIG"]
    for c in cln_candidates:
        if c in df.columns:
            s = df[c].astype(str).replace({"nan": "-", "None": "-", "": "-"})
            if (s != "-").sum() > 0:
                return s
    return pd.Series("-", index=df.index)

def _extract_gnomadv4_af_grpmax(df):
    """
    Extracts gnomADv4 AF grpmax with graceful fallback.
    Returns: (formatted_series, numeric_series, is_fallback_used)
    """
    has_gnomadv4 = ("gnomADv4_AF_grpmax_joint" in df.columns) and (df["gnomADv4_AF_grpmax_joint"].notna().sum() > 0)
    
    if has_gnomadv4:
        af_series = pd.to_numeric(df["gnomADv4_AF_grpmax_joint"], errors="coerce")
        pop_series = df.get("gnomADv4_grpmax_joint", pd.Series("", index=df.index)).fillna("").astype(str)
        is_fallback = False
    else:
        # Fallback to MAX_AF -> gnomADg_AF -> gnomADe_AF -> AF
        af_series = pd.Series(np.nan, index=df.index)
        for c in ["MAX_AF", "gnomADg_AF", "gnomADe_AF", "AF"]:
            if c in df.columns:
                s = pd.to_numeric(df[c], errors="coerce")
                af_series = af_series.fillna(s)
        pop_series = df.get("MAX_AF_POPS", pd.Series("", index=df.index)).fillna("").astype(str)
        is_fallback = True

    formatted = []
    for idx, val in af_series.items():
        if pd.isna(val) or val == 0:
            formatted.append("-")
        else:
            raw_pop = str(pop_series.get(idx, "")).replace("gnomADe_", "").replace("gnomADg_", "").strip()
            pop_tag = f" ({raw_pop})" if raw_pop and raw_pop not in ["nan", "None", ""] else ""
            if val < 0.001:
                formatted.append(f"{val:.2e}{pop_tag}")
            else:
                formatted.append(f"{val:.4f}{pop_tag}")

    return pd.Series(formatted, index=df.index), af_series, is_fallback

def render_table_html(df_subset, col_defs, table_id):
    """
    Dynamically generates the complete <table> HTML with perfectly aligned <thead> and <tbody>
    """
    headers = "".join(f"<th>{h}</th>" for h, _, _ in col_defs)
    rows = []
    for _, r in df_subset.iterrows():
        cells = []
        for header, key, fmt in col_defs:
            raw_val = r.get(key, None)
            if fmt == "locus":
                loc = _safe_str(raw_val, 24)
                parts = str(raw_val).split("-")[0].split(":")
                if len(parts) == 2:
                    link = f"https://gnomad.broadinstitute.org/variant/{parts[0]}-{parts[1]}-{r.get('REF','')}-{r.get('ALT','')}?dataset=gnomad_r4"
                else:
                    link = "#"
                cells.append(f'<td><strong><a href="{link}" target="_blank" style="color:#2563eb; text-decoration:none;">{loc}</a></strong></td>')
            elif fmt == "tier":
                tier_str = str(raw_val)
                badge = "badge-tier4"
                if "Tier 1" in tier_str: badge = "badge-tier1"
                elif "Tier 2" in tier_str: badge = "badge-tier2"
                elif "Tier 3" in tier_str: badge = "badge-tier3"
                cells.append(f'<td><span class="badge {badge}">{tier_str.split(" (")[0]}</span></td>')
            elif fmt == "impact":
                imp = _safe_str(raw_val, 15)
                cells.append(f'<td><span class="impact-{imp.lower()}">{imp}</span></td>')
            elif fmt == "code":
                cells.append(f'<td><code>{_safe_str(raw_val, 32)}</code></td>')
            elif fmt == "score1":
                cells.append(f'<td><strong>{_safe_float(raw_val, 1)}</strong></td>')
            elif fmt == "float2":
                cells.append(f'<td>{_safe_float(raw_val, 2)}</td>')
            elif fmt == "float6":
                cells.append(f'<td>{_safe_float(raw_val, 6)}</td>')
            elif fmt == "bold_float2":
                cells.append(f'<td><strong>{_safe_float(raw_val, 2)}</strong></td>')
            elif fmt == "clinvar":
                cln = _safe_str(raw_val, 35)
                if "pathogenic" in cln.lower():
                    cells.append(f'<td><strong style="color:#dc2626;">{cln}</strong></td>')
                elif "benign" in cln.lower():
                    cells.append(f'<td><span style="color:#16a34a;">{cln}</span></td>')
                else:
                    cells.append(f'<td>{cln}</td>')
            else:
                cells.append(f'<td>{_safe_str(raw_val, 30)}</td>')
        rows.append(f"<tr>{''.join(cells)}</tr>")
    
    tbody = "\n".join(rows)
    return f"""
    <table id="{table_id}">
        <thead>
            <tr>
                {headers}
            </tr>
        </thead>
        <tbody>
            {tbody}
        </tbody>
    </table>
    """

def generate_gene_report(pq_path, out_html_path):
    logger.info(f"Loading {pq_path} for 4-tab clinical report generation...")
    df = pd.read_parquet(pq_path)
    
    gene_name = os.path.basename(pq_path).split(".")[0]
    if "SYMBOL" in df.columns and df["SYMBOL"].dropna().nunique() > 0:
        gene_name = df["SYMBOL"].dropna().iloc[0]

    n_total = len(df)
    logger.info(f"Processing {gene_name} ({n_total:,} variants)...")

    # Gene Meta
    pli_val = df["pLI_gene_value"].dropna().iloc[0] if "pLI_gene_value" in df.columns and len(df["pLI_gene_value"].dropna()) > 0 else "N/A"
    gene_prio = df["gene_priority"].dropna().iloc[0] if "gene_priority" in df.columns and len(df["gene_priority"].dropna()) > 0 else "N/A"
    mane_tx = df["MANE_SELECT"].dropna().iloc[0] if "MANE_SELECT" in df.columns and len(df["MANE_SELECT"].dropna()) > 0 else (df["Feature"].dropna().iloc[0] if "Feature" in df.columns and len(df["Feature"].dropna()) > 0 else "-")

    # Unify ClinVar Display Column
    df["CLINVAR_DISPLAY"] = _extract_clinvar_display(df)

    # Unify gnomADv4 AF grpmax Display Column with Fallback Detection
    df["GNOMADV4_AF_GRPMAX_DISPLAY"], af_numeric, is_af_fallback = _extract_gnomadv4_af_grpmax(df)

    # --- Tool Availability & Disclaimers ---
    disclaimers = []
    
    # SpliceAI Provenance
    has_custom_spliceai = ("spliceai_custom_MAX" in df.columns) and (df["spliceai_custom_MAX"].notna().sum() > 0)
    if has_custom_spliceai:
        disclaimers.append(('provenance', '🧬 <strong>Custom SpliceAI Active</strong>: Evaluated with local <strong>±10,000 bp (20 kb) context window</strong> (-D 10000) for comprehensive deep intronic & non-canonical junction discovery.'))
    else:
        disclaimers.append(('info', 'ℹ️ <strong>VEP SpliceAI Active</strong>: Evaluated with Illumina precomputed lookup table (<strong>500 bp window</strong>). Custom 20kb inference was not computed for this batch.'))

    # gnomADv4 AF grpmax Provenance / Fallback Disclaimer
    if is_af_fallback:
        disclaimers.append(('warning', '⚠️ <strong>gnomADv4 Joint grpmax</strong>: Column not present or unindexed for this dataset; fallback applied to VEP global PopMax AF (MAX_AF).'))
    else:
        disclaimers.append(('provenance', '📊 <strong>gnomAD v4.1 Active</strong>: PopMax allele frequencies extracted from joint exomes + genomes callset (AF_grpmax_joint).'))

    # Tool checks
    tool_status_checks = [
        ("AlphaMissense", "am_pathogenicity" in df.columns and df["am_pathogenicity"].notna().sum() > 0),
        ("REVEL", "REVEL_score" in df.columns and df["REVEL_score"].notna().sum() > 0),
        ("SPiP", ("SPiP_prediction" in df.columns or "SPiP" in df.columns) and (df.get("SPiP_prediction", df.get("SPiP", pd.Series(dtype=float))).notna().sum() > 0)),
        ("SpliceVault", "SpliceVault_status" in df.columns and (get_col_string(df, "SpliceVault_status") == "aberrant_event_detected").sum() > 0),
        ("Pangolin", "Pangolin_max_score" in df.columns and df["Pangolin_max_score"].notna().sum() > 0),
        ("Branchpointer / LaBranchoR", "Branchpoint_status" in df.columns or "LaBranchoR_score" in df.columns),
        ("UTRAnnotator (5' UTR)", "5UTR_consequence" in df.columns and df["5UTR_consequence"].notna().sum() > 0),
        ("UTR.annotation (Kozak/uAUG/PolyA)", any(c in df.columns for c in ["lost_start_codon", "lost_stop_codon", "utr_num_kozak_gainedOrLost", "mrl_gainedOrLost"])),
    ]
    for tool_name, is_present in tool_status_checks:
        if not is_present:
            disclaimers.append(('warning', f'⚠️ <strong>{tool_name}</strong>: Column not present or unindexed for this gene dataset; prioritization applied graceful fallback.'))

    # Numeric Extractions
    splice_custom = get_col_numeric(df, "spliceai_custom_MAX", 0.0)
    splice_vep = get_col_numeric(df, "spliceAI_MAX", 0.0)
    splice_score = np.where(splice_custom > 0, splice_custom, splice_vep)
    df["SPLICE_MAX_UNIFIED"] = splice_score
    
    spip_score = get_col_numeric(df, "SPiP_prediction", get_col_numeric(df, "SPiP", 0.0))
    pangolin_score = get_col_numeric(df, "Pangolin_max_score", 0.0)
    revel_score = get_col_numeric(df, "REVEL_score", 0.0)
    am_score = get_col_numeric(df, "am_pathogenicity", 0.0)
    cadd_score = get_col_numeric(df, "CADD_PHRED", get_col_numeric(df, "CADD_phred", 0.0))
    intron_offset = get_col_numeric(df, "intron_offset_signed", 0.0)
    prio_score = get_col_numeric(df, "VARIANT_PRIORITY_SCORE", 0.0)

    conseq = get_col_string(df, "Consequence", "").str.lower()
    tier_col = get_col_string(df, "PRIORITY_TIER", "Tier 4 (Benign / Tolerated)")
    impact_col = get_col_string(df, "NEW_IMPACT", get_col_string(df, "IMPACT", "MODIFIER"))

    # Global Tier Counts
    n_tier1 = (tier_col.str.contains("Tier 1")).sum()
    n_tier2 = (tier_col.str.contains("Tier 2")).sum()
    n_tier3 = (tier_col.str.contains("Tier 3")).sum()
    n_tier4 = (tier_col.str.contains("Tier 4")).sum()

    import plotly.express as px
    import plotly.graph_objects as go

    # =========================================================================
    # TAB 1: GENERAL / OVERVIEW CHARTS
    # =========================================================================
    tier_df = pd.DataFrame({
        "Tier": ["Tier 1 (Critical Pathogenic)", "Tier 2 (Likely Deleterious)", "Tier 3 (VUS / Moderate)", "Tier 4 (Benign / Tolerated)"],
        "Count": [n_tier1, n_tier2, n_tier3, n_tier4]
    })
    fig_donut = px.pie(
        tier_df, values="Count", names="Tier",
        color="Tier",
        color_discrete_map={
            "Tier 1 (Critical Pathogenic)": "#dc2626",
            "Tier 2 (Likely Deleterious)": "#ea580c",
            "Tier 3 (VUS / Moderate)": "#d97706",
            "Tier 4 (Benign / Tolerated)": "#16a34a"
        },
        hole=0.45,
        title=f"Clinical Triage Stratification ({n_total:,} Variants)"
    )
    fig_donut.update_traces(textposition='inside', textinfo='percent+label')
    fig_donut.update_layout(showlegend=False, margin=dict(t=40, b=20, l=20, r=20), height=330)

    # Top Candidate Table (Tier 1 & 2)
    top_overall_df = df[tier_col.str.contains("Tier 1|Tier 2") | (prio_score >= 35.0)].sort_values("VARIANT_PRIORITY_SCORE", ascending=False).head(200)
    general_cols = [
        ("Locus (GRCh38)", "Locus", "locus"),
        ("HGVSc", "HGVSc", "code"),
        ("HGVSp", "HGVSp", "code"),
        ("Priority Tier", "PRIORITY_TIER", "tier"),
        ("Score", "VARIANT_PRIORITY_SCORE", "score1"),
        ("NEW_IMPACT", "NEW_IMPACT", "impact"),
        ("Consequence", "Consequence", "str"),
        ("gnomADv4 AF grpmax", "GNOMADV4_AF_GRPMAX_DISPLAY", "str"),
        ("SpliceAI Δ", "SPLICE_MAX_UNIFIED", "bold_float2"),
        ("AlphaMissense", "am_pathogenicity", "bold_float2"),
        ("REVEL", "REVEL_score", "float2"),
        ("ClinVar", "CLINVAR_DISPLAY", "clinvar"),
    ]
    general_table_html = render_table_html(top_overall_df, general_cols, "generalTable")

    # =========================================================================
    # TAB 2: SPLICING TAB
    # =========================================================================
    splice_mask = (
        (splice_score >= 0.20) |
        (spip_score >= 0.20) |
        (pangolin_score >= 0.20) |
        (get_col_string(df, "SpliceVault_status") == "aberrant_event_detected") |
        (get_col_string(df, "Branchpoint_status") == "disrupted") |
        (conseq.str.contains("splice"))
    )
    df_splice = df[splice_mask].copy().sort_values("VARIANT_PRIORITY_SCORE", ascending=False)
    n_splice_total = len(df_splice)
    n_splice_high = (splice_score >= 0.50).sum()
    n_splice_deep = ((splice_score >= 0.20) & (intron_offset.abs() > 500)).sum()
    n_splice_vault = (get_col_string(df, "SpliceVault_status") == "aberrant_event_detected").sum()
    n_branchpoint = (get_col_string(df, "Branchpoint_status") == "disrupted").sum()

    # Splicing Offset Scatter Plot
    if len(df_splice) > 0:
        df_splice_plot = df_splice.copy()
        df_splice_plot["Splice_Score"] = df_splice_plot["SPLICE_MAX_UNIFIED"].astype(float)
        df_splice_plot["Offset"] = get_col_numeric(df_splice_plot, "intron_offset_signed", 0.0)
        df_splice_plot["Tier_Label"] = get_col_string(df_splice_plot, "PRIORITY_TIER", "Tier 4")

        fig_splice_scatter = px.scatter(
            df_splice_plot,
            x="Offset", y="Splice_Score",
            color="Tier_Label",
            color_discrete_map={
                "Tier 1 (Critical Pathogenic Candidate)": "#dc2626",
                "Tier 2 (Likely Deleterious / Strong Candidate)": "#ea580c",
                "Tier 3 (VUS / Moderate Potential)": "#d97706",
                "Tier 4 (Benign / Tolerated)": "#94a3b8"
            },
            hover_data=["HGVSc", "Consequence"],
            labels={"Offset": "Distance to Splice Site (bp, signed)", "Splice_Score": "SpliceAI Max Δ Score"},
            title="SpliceAI Δ vs Intronic Distance"
        )
        fig_splice_scatter.add_hline(y=0.50, line_dash="dash", line_color="#dc2626", annotation_text="High SpliceAI ≥ 0.50")
        fig_splice_scatter.add_hline(y=0.20, line_dash="dot", line_color="#d97706", annotation_text="Moderate SpliceAI ≥ 0.20")
        fig_splice_scatter.update_layout(margin=dict(t=40, b=20, l=20, r=20), height=330)
    else:
        fig_splice_scatter = go.Figure()
        fig_splice_scatter.update_layout(title="No Splicing Variants Detected", height=330)

    splice_table_cols = [
        ("Locus (GRCh38)", "Locus", "locus"),
        ("HGVSc", "HGVSc", "code"),
        ("HGVSp", "HGVSp", "code"),
        ("Consequence", "Consequence", "str"),
        ("Offset (bp)", "intron_offset_signed", "str"),
        ("Priority Tier", "PRIORITY_TIER", "tier"),
        ("Score", "VARIANT_PRIORITY_SCORE", "score1"),
        ("gnomADv4 AF grpmax", "GNOMADV4_AF_GRPMAX_DISPLAY", "str"),
        ("SpliceAI Custom Δ", "SPLICE_MAX_UNIFIED", "bold_float2"),
        ("DS_AG", "spliceai_custom_DS_AG", "float2"),
        ("DS_AL", "spliceai_custom_DS_AL", "float2"),
        ("DS_DG", "spliceai_custom_DS_DG", "float2"),
        ("DS_DL", "spliceai_custom_DS_DL", "float2"),
        ("SPiP Pred", "SPiP_prediction", "float2"),
        ("SpliceVault", "SpliceVault_status", "str"),
        ("Branchpoint", "Branchpoint_status", "str"),
        ("ClinVar", "CLINVAR_DISPLAY", "clinvar"),
    ]
    splice_table_html = render_table_html(df_splice.head(250), splice_table_cols, "spliceTable")

    # =========================================================================
    # TAB 3: UTR TAB
    # =========================================================================
    utr_mask = (
        (conseq.str.contains("utr|5_prime|3_prime|prime_utr")) |
        (df.get("5UTR_consequence", pd.Series(dtype=str)).notna()) |
        (get_col_string(df, "utr_num_kozak_gainedOrLost").str.lower().isin(["gained", "lost", "true"])) |
        (get_col_string(df, "utr_num_uAUG_gainedOrLost").str.lower().isin(["gained", "lost", "true"]))
    )
    df_utr = df[utr_mask].copy().sort_values("VARIANT_PRIORITY_SCORE", ascending=False)
    n_utr_total = len(df_utr)
    n_5utr_conseq = (df.get("5UTR_consequence", pd.Series(dtype=str)).notna() & (get_col_string(df, "5UTR_consequence") != "nan") & (get_col_string(df, "5UTR_consequence") != "-") & (get_col_string(df, "5UTR_consequence") != "")).sum() if "5UTR_consequence" in df.columns else 0
    n_kozak = (get_col_string(df, "utr_num_kozak_gainedOrLost").str.lower().isin(["gained", "lost", "true"])).sum()
    n_uaug = (get_col_string(df, "utr_num_uAUG_gainedOrLost").str.lower().isin(["gained", "lost", "true"])).sum()
    n_start_stop = (conseq.str.contains("start_lost|stop_lost|lost_start|lost_stop")).sum()

    # UTR Consequence Chart
    if len(df_utr) > 0:
        utr_counts = df_utr["Consequence"].value_counts().head(8).reset_index()
        utr_counts.columns = ["Consequence", "Count"]
        fig_utr_bar = px.bar(
            utr_counts, x="Count", y="Consequence", orientation="h",
            title="UTR Variant Consequence Distribution",
            color="Count", color_continuous_scale="Teal"
        )
        fig_utr_bar.update_layout(yaxis={'categoryorder': 'total ascending'}, margin=dict(t=40, b=20, l=20, r=20), height=330)
    else:
        fig_utr_bar = go.Figure()
        fig_utr_bar.update_layout(title="No UTR Variants Detected", height=330)

    utr_table_cols = [
        ("Locus (GRCh38)", "Locus", "locus"),
        ("HGVSc", "HGVSc", "code"),
        ("Consequence", "Consequence", "str"),
        ("Priority Tier", "PRIORITY_TIER", "tier"),
        ("Score", "VARIANT_PRIORITY_SCORE", "score1"),
        ("NEW_IMPACT", "NEW_IMPACT", "impact"),
        ("5UTR Consequence", "5UTR_consequence", "str"),
        ("Kozak Alt", "utr_num_kozak_gainedOrLost", "str"),
        ("uAUG Alt", "utr_num_uAUG_gainedOrLost", "str"),
        ("MRL Alt", "mrl_gainedOrLost", "str"),
        ("PolyA Alt", "num_polyA_signal_gainedOrLost", "str"),
        ("gnomADv4 AF grpmax", "GNOMADV4_AF_GRPMAX_DISPLAY", "str"),
        ("ClinVar", "CLINVAR_DISPLAY", "clinvar"),
    ]
    utr_table_html = render_table_html(df_utr.head(250), utr_table_cols, "utrTable")

    # =========================================================================
    # TAB 4: MISSENSE TAB
    # =========================================================================
    missense_mask = conseq.str.contains("missense")
    df_missense = df[missense_mask].copy().sort_values("VARIANT_PRIORITY_SCORE", ascending=False)
    n_missense_total = len(df_missense)
    n_am_path = (am_score >= 0.564).sum()
    n_revel_path = (revel_score >= 0.75).sum()
    n_consensus_path = ((am_score >= 0.564) & (revel_score >= 0.75)).sum()
    n_cadd_high = (cadd_score >= 25.0).sum()

    # Missense Quadrant Scatter
    if len(df_missense) > 0:
        df_missense_plot = df_missense[(df_missense["am_pathogenicity"].notna()) | (df_missense["REVEL_score"].notna())].copy()
        df_missense_plot["AM_val"] = get_col_numeric(df_missense_plot, "am_pathogenicity", 0.0)
        df_missense_plot["REVEL_val"] = get_col_numeric(df_missense_plot, "REVEL_score", 0.0)
        df_missense_plot["Tier_Label"] = get_col_string(df_missense_plot, "PRIORITY_TIER", "Tier 4")

        fig_missense_scatter = px.scatter(
            df_missense_plot,
            x="REVEL_val", y="AM_val",
            color="Tier_Label",
            color_discrete_map={
                "Tier 1 (Critical Pathogenic Candidate)": "#dc2626",
                "Tier 2 (Likely Deleterious / Strong Candidate)": "#ea580c",
                "Tier 3 (VUS / Moderate Potential)": "#d97706",
                "Tier 4 (Benign / Tolerated)": "#94a3b8"
            },
            hover_data=["HGVSc", "HGVSp", "Consequence"],
            labels={"REVEL_val": "REVEL Score (0-1)", "AM_val": "AlphaMissense Score (0-1)"},
            title="Missense Pathogenicity (REVEL vs AlphaMissense)"
        )
        fig_missense_scatter.add_vline(x=0.75, line_dash="dash", line_color="#dc2626", annotation_text="REVEL ≥ 0.75")
        fig_missense_scatter.add_hline(y=0.564, line_dash="dash", line_color="#ea580c", annotation_text="AlphaMissense ≥ 0.564")
        fig_missense_scatter.update_layout(margin=dict(t=40, b=20, l=20, r=20), height=330)
    else:
        fig_missense_scatter = go.Figure()
        fig_missense_scatter.update_layout(title="No Missense Variants Detected", height=330)

    missense_table_cols = [
        ("Locus (GRCh38)", "Locus", "locus"),
        ("HGVSc", "HGVSc", "code"),
        ("HGVSp", "HGVSp", "code"),
        ("Priority Tier", "PRIORITY_TIER", "tier"),
        ("Score", "VARIANT_PRIORITY_SCORE", "score1"),
        ("AlphaMissense", "am_pathogenicity", "bold_float2"),
        ("AM Class", "am_class", "str"),
        ("REVEL", "REVEL_score", "bold_float2"),
        ("CADD Phred", "CADD_PHRED", "float2"),
        ("Polyphen2", "Polyphen2_HVAR_pred", "str"),
        ("SIFT", "SIFT_pred", "str"),
        ("gnomADv4 AF grpmax", "GNOMADV4_AF_GRPMAX_DISPLAY", "str"),
        ("ClinVar", "CLINVAR_DISPLAY", "clinvar"),
    ]
    missense_table_html = render_table_html(df_missense.head(250), missense_table_cols, "missenseTable")

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Render Disclaimers HTML
    disclaimer_items = []
    for dtype, msg in disclaimers:
        dclass = "disc-info" if dtype == "info" else ("disc-prov" if dtype == "provenance" else "disc-warn")
        disclaimer_items.append(f'<div class="disclaimer-banner {dclass}">{msg}</div>')
    disclaimers_html = "\n".join(disclaimer_items)

    # Assemble HTML
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clinical Prioritization Dashboard: {gene_name} | GRCh38</title>
    <script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
    <style>
        :root {{
            --bg-color: #f8fafc;
            --card-bg: #ffffff;
            --text-main: #1e293b;
            --text-muted: #64748b;
            --primary: #2563eb;
            --tier1: #dc2626;
            --tier2: #ea580c;
            --tier3: #d97706;
            --tier4: #16a34a;
            --border-color: #e2e8f0;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-main);
            margin: 0;
            padding: 24px;
        }}
        .header-banner {{
            background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 100%);
            color: white;
            padding: 24px 30px;
            border-radius: 12px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.08);
            margin-bottom: 20px;
        }}
        .header-banner h1 {{
            margin: 0 0 6px 0;
            font-size: 26px;
            font-weight: 700;
        }}
        .meta-tags {{
            display: flex;
            gap: 10px;
            flex-wrap: wrap;
            margin-top: 10px;
        }}
        .meta-tag {{
            background: rgba(255,255,255,0.15);
            padding: 4px 10px;
            border-radius: 16px;
            font-size: 12px;
            font-weight: 500;
        }}
        .disclaimers-container {{
            margin-bottom: 20px;
            display: flex;
            flex-direction: column;
            gap: 8px;
        }}
        .disclaimer-banner {{
            padding: 10px 16px;
            border-radius: 8px;
            font-size: 13px;
            line-height: 1.4;
        }}
        .disc-prov {{ background: #ecfdf5; border-left: 4px solid #059669; color: #065f46; }}
        .disc-info {{ background: #eff6ff; border-left: 4px solid #3b82f6; color: #1e40af; }}
        .disc-warn {{ background: #fffbeb; border-left: 4px solid #f59e0b; color: #92400e; }}
        
        /* Navigation Tabs */
        .tab-nav {{
            display: flex;
            gap: 8px;
            border-bottom: 2px solid var(--border-color);
            margin-bottom: 24px;
            background: var(--card-bg);
            padding: 8px 12px 0 12px;
            border-radius: 10px 10px 0 0;
        }}
        .tab-btn {{
            padding: 10px 20px;
            border: none;
            background: transparent;
            font-size: 14px;
            font-weight: 600;
            color: var(--text-muted);
            cursor: pointer;
            border-bottom: 3px solid transparent;
            transition: all 0.2s ease;
            border-radius: 6px 6px 0 0;
        }}
        .tab-btn:hover {{
            color: var(--primary);
            background: #f1f5f9;
        }}
        .tab-btn.active {{
            color: var(--primary);
            border-bottom: 3px solid var(--primary);
            background: #eff6ff;
        }}
        .tab-content {{
            display: none;
        }}
        .tab-content.active {{
            display: block;
        }}

        /* KPI Cards */
        .kpi-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }}
        .kpi-card {{
            background: var(--card-bg);
            padding: 16px 20px;
            border-radius: 10px;
            border: 1px solid var(--border-color);
            box-shadow: 0 2px 6px rgba(0,0,0,0.03);
            text-align: center;
        }}
        .kpi-card.t1 {{ border-top: 4px solid var(--tier1); }}
        .kpi-card.t2 {{ border-top: 4px solid var(--tier2); }}
        .kpi-card.t3 {{ border-top: 4px solid var(--tier3); }}
        .kpi-card.t4 {{ border-top: 4px solid var(--tier4); }}
        .kpi-card.blue {{ border-top: 4px solid var(--primary); }}
        .kpi-title {{
            font-size: 12px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--text-muted);
            margin-bottom: 4px;
        }}
        .kpi-value {{
            font-size: 28px;
            font-weight: 700;
        }}
        .kpi-value.t1 {{ color: var(--tier1); }}
        .kpi-value.t2 {{ color: var(--tier2); }}
        .kpi-value.t3 {{ color: var(--tier3); }}
        .kpi-value.t4 {{ color: var(--tier4); }}
        .kpi-value.blue {{ color: var(--primary); }}

        /* Charts */
        .charts-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
            gap: 20px;
            margin-bottom: 24px;
        }}
        .chart-card {{
            background: var(--card-bg);
            padding: 14px;
            border-radius: 10px;
            border: 1px solid var(--border-color);
            box-shadow: 0 2px 6px rgba(0,0,0,0.03);
        }}

        /* Table */
        .table-section {{
            background: var(--card-bg);
            padding: 20px;
            border-radius: 10px;
            border: 1px solid var(--border-color);
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
            margin-bottom: 30px;
        }}
        .table-header-flex {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 14px;
            flex-wrap: wrap;
            gap: 10px;
        }}
        .search-box {{
            padding: 8px 14px;
            border-radius: 20px;
            border: 1px solid #cbd5e1;
            font-size: 13px;
            width: 260px;
        }}
        .table-container {{
            overflow-x: auto;
            max-height: 550px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 12.5px;
            text-align: left;
        }}
        th {{
            background-color: #f8fafc;
            color: #475569;
            font-weight: 600;
            padding: 9px 12px;
            border-bottom: 2px solid var(--border-color);
            position: sticky;
            top: 0;
            z-index: 10;
        }}
        td {{
            padding: 9px 12px;
            border-bottom: 1px solid #f1f5f9;
            white-space: nowrap;
        }}
        tr:hover {{
            background-color: #f8fafc;
        }}
        .badge {{
            display: inline-block;
            padding: 2px 7px;
            border-radius: 10px;
            font-size: 11px;
            font-weight: 700;
            color: white;
        }}
        .badge-tier1 {{ background-color: var(--tier1); }}
        .badge-tier2 {{ background-color: var(--tier2); }}
        .badge-tier3 {{ background-color: var(--tier3); }}
        .badge-tier4 {{ background-color: var(--tier4); }}
        .impact-high {{ color: #dc2626; font-weight: 700; }}
        .impact-moderate {{ color: #ea580c; font-weight: 600; }}
        .impact-low {{ color: #16a34a; }}
        .impact-modifier {{ color: #64748b; }}
        .footer {{
            text-align: center;
            font-size: 12px;
            color: var(--text-muted);
            margin-top: 30px;
        }}
    </style>
</head>
<body>
    <div class="header-banner">
        <h1>🧬 Gene Clinical Prioritization Dashboard: {gene_name}</h1>
        <div>Transcript: <strong>{mane_tx}</strong> | Genome Build: <strong>GRCh38</strong> | Generated: <strong>{timestamp}</strong></div>
        <div class="meta-tags">
            <span class="meta-tag">Gene Constraint (pLI): <strong>{pli_val}</strong></span>
            <span class="meta-tag">ClinGen Priority: <strong>{gene_prio}</strong></span>
            <span class="meta-tag">Total Cohort Variants: <strong>{n_total:,}</strong></span>
            <span class="meta-tag">Splice Alterations: <strong>{n_splice_total:,}</strong></span>
            <span class="meta-tag">UTR Alterations: <strong>{n_utr_total:,}</strong></span>
            <span class="meta-tag">Missense Variants: <strong>{n_missense_total:,}</strong></span>
        </div>
    </div>

    <!-- Disclaimers & Provenance -->
    <div class="disclaimers-container">
        {disclaimers_html}
    </div>

    <!-- Tab Navigation -->
    <div class="tab-nav">
        <button class="tab-btn active" onclick="switchTab('tab-general')">📊 General / Overview</button>
        <button class="tab-btn" onclick="switchTab('tab-splicing')">🧬 Splicing Alterations ({n_splice_total:,})</button>
        <button class="tab-btn" onclick="switchTab('tab-utr')">🎯 5'/3' UTR & Translation ({n_utr_total:,})</button>
        <button class="tab-btn" onclick="switchTab('tab-missense')">🔬 Missense Pathogenicity ({n_missense_total:,})</button>
    </div>

    <!-- ==================== TAB 1: GENERAL ==================== -->
    <div id="tab-general" class="tab-content active">
        <div class="kpi-grid">
            <div class="kpi-card t1">
                <div class="kpi-title">Tier 1: Critical Pathogenic</div>
                <div class="kpi-value t1">{n_tier1:,}</div>
                <small>Immediate diagnostic review</small>
            </div>
            <div class="kpi-card t2">
                <div class="kpi-title">Tier 2: Likely Deleterious</div>
                <div class="kpi-value t2">{n_tier2:,}</div>
                <small>Strong candidate for validation</small>
            </div>
            <div class="kpi-card t3">
                <div class="kpi-title">Tier 3: VUS / Moderate</div>
                <div class="kpi-value t3">{n_tier3:,}</div>
                <small>Secondary research candidate</small>
            </div>
            <div class="kpi-card t4">
                <div class="kpi-title">Tier 4: Benign / Tolerated</div>
                <div class="kpi-value t4">{n_tier4:,}</div>
                <small>Common / non-damaging</small>
            </div>
        </div>

        <div class="charts-grid">
            <div class="chart-card">
                {fig_donut.to_html(full_html=False, include_plotlyjs=False)}
            </div>
            <div class="chart-card">
                {fig_missense_scatter.to_html(full_html=False, include_plotlyjs=False)}
            </div>
            <div class="chart-card">
                {fig_splice_scatter.to_html(full_html=False, include_plotlyjs=False)}
            </div>
        </div>

        <div class="table-section">
            <div class="table-header-flex">
                <div>
                    <h3 style="margin:0 0 2px 0;">Top Prioritized Candidate Variants ({len(top_overall_df)} Candidates)</h3>
                    <small style="color:var(--text-muted);">Includes Tier 1, Tier 2, and high-scoring Tier 3 candidates ranked by multi-evidence score</small>
                </div>
                <input type="text" class="search-box" placeholder="🔍 Search table..." onkeyup="filterTable(this, 'generalTable')">
            </div>
            <div class="table-container">
                {general_table_html}
            </div>
        </div>
    </div>

    <!-- ==================== TAB 2: SPLICING ==================== -->
    <div id="tab-splicing" class="tab-content">
        <div class="kpi-grid">
            <div class="kpi-card blue">
                <div class="kpi-title">Total Splice Candidates</div>
                <div class="kpi-value blue">{n_splice_total:,}</div>
                <small>SpliceAI / SPiP / Pangolin ≥ 0.20</small>
            </div>
            <div class="kpi-card t1">
                <div class="kpi-title">High-Impact SpliceAI (≥ 0.50)</div>
                <div class="kpi-value t1">{n_splice_high:,}</div>
                <small>High confidence splice disruption</small>
            </div>
            <div class="kpi-card t2">
                <div class="kpi-title">Deep Intronic Events (>500bp)</div>
                <div class="kpi-value t2">{n_splice_deep:,}</div>
                <small>Custom 20kb discovery</small>
            </div>
            <div class="kpi-card blue">
                <div class="kpi-title">SpliceVault RNA Events</div>
                <div class="kpi-value blue">{n_splice_vault:,}</div>
                <small>Empirical RNA mis-splicing</small>
            </div>
            <div class="kpi-card blue">
                <div class="kpi-title">Branchpoint Disruptions</div>
                <div class="kpi-value blue">{n_branchpoint:,}</div>
                <small>Branchpointer / LaBranchoR</small>
            </div>
        </div>

        <div class="charts-grid">
            <div class="chart-card" style="grid-column: 1 / -1;">
                {fig_splice_scatter.to_html(full_html=False, include_plotlyjs=False)}
            </div>
        </div>

        <div class="table-section">
            <div class="table-header-flex">
                <div>
                    <h3 style="margin:0 0 2px 0;">Splicing Altering Variants ({len(df_splice)} Candidates)</h3>
                    <small style="color:var(--text-muted);">Displaying only variants with predicted splicing effects or junction disruptions</small>
                </div>
                <input type="text" class="search-box" placeholder="🔍 Search splicing table..." onkeyup="filterTable(this, 'spliceTable')">
            </div>
            <div class="table-container">
                {splice_table_html}
            </div>
        </div>
    </div>

    <!-- ==================== TAB 3: UTR ==================== -->
    <div id="tab-utr" class="tab-content">
        <div class="kpi-grid">
            <div class="kpi-card blue">
                <div class="kpi-title">Total UTR Variants</div>
                <div class="kpi-value blue">{n_utr_total:,}</div>
                <small>5' & 3' UTR sequences</small>
            </div>
            <div class="kpi-card t1">
                <div class="kpi-title">5' UTR Translation Consequences</div>
                <div class="kpi-value t1">{n_5utr_conseq:,}</div>
                <small>5UTR.annotator hits</small>
            </div>
            <div class="kpi-card t2">
                <div class="kpi-title">Kozak / uAUG Alterations</div>
                <div class="kpi-value t2">{n_kozak + n_uaug:,}</div>
                <small>Upstream ORF / Kozak shifts</small>
            </div>
            <div class="kpi-card t1">
                <div class="kpi-title">Start / Stop Codon Lost</div>
                <div class="kpi-value t1">{n_start_stop:,}</div>
                <small>Initiation / termination loss</small>
            </div>
        </div>

        <div class="charts-grid">
            <div class="chart-card" style="grid-column: 1 / -1;">
                {fig_utr_bar.to_html(full_html=False, include_plotlyjs=False)}
            </div>
        </div>

        <div class="table-section">
            <div class="table-header-flex">
                <div>
                    <h3 style="margin:0 0 2px 0;">5' & 3' UTR Variants ({len(df_utr)} Candidates)</h3>
                    <small style="color:var(--text-muted);">Displaying only non-coding regulatory and translation-modifying UTR variants</small>
                </div>
                <input type="text" class="search-box" placeholder="🔍 Search UTR table..." onkeyup="filterTable(this, 'utrTable')">
            </div>
            <div class="table-container">
                {utr_table_html}
            </div>
        </div>
    </div>

    <!-- ==================== TAB 4: MISSENSE ==================== -->
    <div id="tab-missense" class="tab-content">
        <div class="kpi-grid">
            <div class="kpi-card blue">
                <div class="kpi-title">Total Missense Variants</div>
                <div class="kpi-value blue">{n_missense_total:,}</div>
                <small>Non-synonymous coding</small>
            </div>
            <div class="kpi-card t1">
                <div class="kpi-title">Consensus Pathogenic (AM+REVEL)</div>
                <div class="kpi-value t1">{n_consensus_path:,}</div>
                <small>AM ≥ 0.564 & REVEL ≥ 0.75</small>
            </div>
            <div class="kpi-card t2">
                <div class="kpi-title">AlphaMissense Pathogenic</div>
                <div class="kpi-value t2">{n_am_path:,}</div>
                <small>DeepMind AM ≥ 0.564</small>
            </div>
            <div class="kpi-card t2">
                <div class="kpi-title">REVEL Pathogenic</div>
                <div class="kpi-value t2">{n_revel_path:,}</div>
                <small>Ensemble REVEL ≥ 0.75</small>
            </div>
            <div class="kpi-card blue">
                <div class="kpi-title">CADD Phred ≥ 25</div>
                <div class="kpi-value blue">{n_cadd_high:,}</div>
                <small>Top 0.3% conserved</small>
            </div>
        </div>

        <div class="charts-grid">
            <div class="chart-card" style="grid-column: 1 / -1;">
                {fig_missense_scatter.to_html(full_html=False, include_plotlyjs=False)}
            </div>
        </div>

        <div class="table-section">
            <div class="table-header-flex">
                <div>
                    <h3 style="margin:0 0 2px 0;">Missense Variants ({len(df_missense)} Candidates)</h3>
                    <small style="color:var(--text-muted);">Displaying only non-synonymous amino acid substitutions</small>
                </div>
                <input type="text" class="search-box" placeholder="🔍 Search missense table..." onkeyup="filterTable(this, 'missenseTable')">
            </div>
            <div class="table-container">
                {missense_table_html}
            </div>
        </div>
    </div>

    <script>
        function switchTab(tabId) {{
            document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
            document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));
            
            document.getElementById(tabId).classList.add('active');
            event.currentTarget.classList.add('active');
            
            window.dispatchEvent(new Event('resize'));
        }}

        function filterTable(input, tableId) {{
            let filter = input.value.toLowerCase();
            let table = document.getElementById(tableId);
            let tr = table.getElementsByTagName("tr");

            for (let i = 1; i < tr.length; i++) {{
                let text = tr[i].textContent.toLowerCase();
                if (text.indexOf(filter) > -1) {{
                    tr[i].style.display = "";
                }} else {{
                    tr[i].style.display = "none";
                }}
            }}
        }}
    </script>

    <div class="footer">
        Generated by <strong>Cardiovascular Multi-Omics Annotation Pipeline</strong> • Version 2026.08 • GRCh38
    </div>
</body>
</html>
"""

    out_dir = os.path.dirname(out_html_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(out_html_path, "w") as f:
        f.write(html_content)

    logger.info(f"Successfully generated 4-tab clinical dashboard at: {out_html_path}")
    return True

def main():
    parser = argparse.ArgumentParser(description="Generate 4-Tab Medical Clinical Prioritization Dashboard per Gene")
    parser.add_argument("--input", help="Path to single gene parquet file (.parsed.clean.pq)")
    parser.add_argument("--output", help="Path for output HTML dashboard")
    parser.add_argument("--run-dir", default="RUNS/run_20260813_1028", help="Run directory containing results/ folder to process in batch")
    parser.add_argument("--all", action="store_true", help="Generate reports for all completed gene parquets in run-dir")
    args = parser.parse_args()

    if args.input:
        out_path = args.output if args.output else os.path.join(os.path.dirname(args.input), "..", "reports", f"{os.path.basename(args.input).split('.')[0]}_clinical_prioritization_report.html")
        generate_gene_report(args.input, out_path)
    elif args.all or args.run_dir:
        res_dir = os.path.join(args.run_dir, "results")
        rep_dir = os.path.join(args.run_dir, "reports")
        os.makedirs(rep_dir, exist_ok=True)
        
        pq_files = sorted(glob.glob(os.path.join(res_dir, "*.parsed.clean.pq")))
        logger.info(f"Found {len(pq_files)} gene parquets in {res_dir} to process...")
        
        for pq in pq_files:
            gname = os.path.basename(pq).split(".")[0]
            out_html = os.path.join(rep_dir, f"{gname}_clinical_prioritization_report.html")
            try:
                generate_gene_report(pq, out_html)
            except Exception as e:
                logger.error(f"Error generating clinical report for {gname}: {e}")

if __name__ == "__main__":
    main()
