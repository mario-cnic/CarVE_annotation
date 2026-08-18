#!/usr/bin/env python3
"""
Medical & Geneticist-Grade Clinical Prioritization Report Generator

Generates a self-contained, interactive HTML dashboard for a specific gene
documenting prioritized pathogenic candidates across Tier 1, Tier 2, and Tier 3,
multi-omics evidence plots, and an interactive searchable clinical table.
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
import json
import pandas as pd
import numpy as np
from datetime import datetime

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("clinical_report_generator")

def _safe_str(val, max_len=40):
    if pd.isna(val) or val is None or str(val).strip() == "" or str(val).lower() == "nan":
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

def generate_gene_report(pq_path, out_html_path):
    logger.info(f"Loading {pq_path} for clinical report generation...")
    df = pd.read_parquet(pq_path)
    
    gene_name = os.path.basename(pq_path).split(".")[0]
    if "SYMBOL" in df.columns and df["SYMBOL"].dropna().nunique() > 0:
        gene_name = df["SYMBOL"].dropna().iloc[0]

    n_total = len(df)
    logger.info(f"Processing {gene_name} ({n_total:,} variants)...")

    # Extract Key Metrics
    pli_val = df["pLI_gene_value"].dropna().iloc[0] if "pLI_gene_value" in df.columns and len(df["pLI_gene_value"].dropna()) > 0 else "N/A"
    gene_prio = df["gene_priority"].dropna().iloc[0] if "gene_priority" in df.columns and len(df["gene_priority"].dropna()) > 0 else "N/A"
    mane_tx = df["MANE_SELECT"].dropna().iloc[0] if "MANE_SELECT" in df.columns and len(df["MANE_SELECT"].dropna()) > 0 else (df["Feature"].dropna().iloc[0] if "Feature" in df.columns and len(df["Feature"].dropna()) > 0 else "-")

    # Tier Counts
    tier_col = df["PRIORITY_TIER"].astype(str) if "PRIORITY_TIER" in df.columns else pd.Series("Tier 4 (Benign / Tolerated)", index=df.index)
    n_tier1 = (tier_col.str.contains("Tier 1")).sum()
    n_tier2 = (tier_col.str.contains("Tier 2")).sum()
    n_tier3 = (tier_col.str.contains("Tier 3")).sum()
    n_tier4 = (tier_col.str.contains("Tier 4")).sum()

    # High Splicing & High Missense counts
    splice_score = pd.to_numeric(df.get("spliceai_custom_MAX", df.get("spliceAI_MAX", 0)), errors="coerce").fillna(0)
    revel_score = pd.to_numeric(df.get("REVEL_score", 0), errors="coerce").fillna(0)
    am_score = pd.to_numeric(df.get("am_pathogenicity", 0), errors="coerce").fillna(0)
    sv_events = (df.get("SpliceVault_status", "") == "aberrant_event_detected").sum() if "SpliceVault_status" in df.columns else 0
    bp_disrupt = (df.get("Branchpoint_status", "") == "disrupted").sum() if "Branchpoint_status" in df.columns else 0

    n_splice_path = (splice_score >= 0.50).sum()
    n_missense_path = ((am_score >= 0.564) & (revel_score >= 0.75)).sum()

    # Interactive Plots generation via Plotly
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # 1. Tier Donut Chart
    tier_df = pd.DataFrame({
        "Tier": ["Tier 1 (Critical Pathogenic)", "Tier 2 (Likely Deleterious)", "Tier 3 (VUS / Moderate)", "Tier 4 (Benign / Tolerated)"],
        "Count": [n_tier1, n_tier2, n_tier3, n_tier4],
        "Color": ["#e74c3c", "#e67e22", "#f1c40f", "#2ecc71"]
    })
    fig_donut = px.pie(
        tier_df,
        values="Count",
        names="Tier",
        color="Tier",
        color_discrete_map={
            "Tier 1 (Critical Pathogenic)": "#e74c3c",
            "Tier 2 (Likely Deleterious)": "#e67e22",
            "Tier 3 (VUS / Moderate)": "#f1c40f",
            "Tier 4 (Benign / Tolerated)": "#2ecc71"
        },
        hole=0.45,
        title=f"Clinical Triage Stratification ({n_total:,} variants)"
    )
    fig_donut.update_traces(textposition='inside', textinfo='percent+label')
    fig_donut.update_layout(showlegend=False, margin=dict(t=40, b=20, l=20, r=20), height=340)

    # 2. Pathogenicity Scatter (AlphaMissense vs REVEL, color by Tier)
    sub_missense = df[(revel_score > 0) | (am_score > 0)].copy()
    if len(sub_missense) > 0:
        sub_missense["AM_val"] = pd.to_numeric(sub_missense.get("am_pathogenicity", 0), errors="coerce").fillna(0)
        sub_missense["REVEL_val"] = pd.to_numeric(sub_missense.get("REVEL_score", 0), errors="coerce").fillna(0)
        sub_missense["Tier_Label"] = sub_missense["PRIORITY_TIER"].astype(str)
        sub_missense["HGVSc_clean"] = sub_missense["HGVSc"].astype(str).str.split(":").str[-1]
        sub_missense["HGVSp_clean"] = sub_missense["HGVSp"].astype(str).str.split(":").str[-1]

        fig_missense = px.scatter(
            sub_missense,
            x="REVEL_val",
            y="AM_val",
            color="Tier_Label",
            color_discrete_map={
                "Tier 1 (Critical Pathogenic Candidate)": "#e74c3c",
                "Tier 2 (Likely Deleterious / Strong Candidate)": "#e67e22",
                "Tier 3 (VUS / Moderate Potential)": "#f1c40f",
                "Tier 4 (Benign / Tolerated)": "#95a5a6"
            },
            hover_data=["HGVSc_clean", "HGVSp_clean", "Consequence"],
            labels={"REVEL_val": "REVEL Score (0-1)", "AM_val": "AlphaMissense Pathogenicity (0-1)"},
            title="Missense Pathogenicity Correlation (REVEL vs AlphaMissense)"
        )
        fig_missense.add_vline(x=0.75, line_dash="dash", line_color="#e74c3c", annotation_text="REVEL ≥ 0.75")
        fig_missense.add_hline(y=0.564, line_dash="dash", line_color="#e67e22", annotation_text="AlphaMissense ≥ 0.564")
        fig_missense.update_layout(margin=dict(t=40, b=20, l=20, r=20), height=340)
    else:
        fig_missense = go.Figure()
        fig_missense.update_layout(title="No Missense Predictors Available", height=340)

    # 3. Splicing Landscape (Signed Intron Offset vs SpliceAI Custom Δ)
    sub_splice = df[splice_score >= 0.10].copy()
    if len(sub_splice) > 0:
        sub_splice["Splice_Score"] = pd.to_numeric(sub_splice.get("spliceai_custom_MAX", sub_splice.get("spliceAI_MAX", 0)), errors="coerce").fillna(0)
        sub_splice["Offset"] = pd.to_numeric(sub_splice.get("intron_offset_signed", 0), errors="coerce").fillna(0)
        sub_splice["Tier_Label"] = sub_splice["PRIORITY_TIER"].astype(str)
        sub_splice["HGVSc_clean"] = sub_splice["HGVSc"].astype(str).str.split(":").str[-1]

        fig_splice = px.scatter(
            sub_splice,
            x="Offset",
            y="Splice_Score",
            color="Tier_Label",
            color_discrete_map={
                "Tier 1 (Critical Pathogenic Candidate)": "#e74c3c",
                "Tier 2 (Likely Deleterious / Strong Candidate)": "#e67e22",
                "Tier 3 (VUS / Moderate Potential)": "#f1c40f",
                "Tier 4 (Benign / Tolerated)": "#95a5a6"
            },
            hover_data=["HGVSc_clean", "Consequence"],
            labels={"Offset": "Intron Offset (bp from Junction, signed)", "Splice_Score": "SpliceAI Custom Δ Score"},
            title="Splicing Alteration Landscape vs Intronic Distance"
        )
        fig_splice.add_hline(y=0.50, line_dash="dash", line_color="#e74c3c", annotation_text="High SpliceAI ≥ 0.50")
        fig_splice.add_hline(y=0.20, line_dash="dot", line_color="#f1c40f", annotation_text="Moderate SpliceAI ≥ 0.20")
        fig_splice.update_layout(margin=dict(t=40, b=20, l=20, r=20), height=340)
    else:
        fig_splice = go.Figure()
        fig_splice.update_layout(title="No Splicing Alterations Detected", height=340)

    # 4. Filter High-Priority Candidate Variants for the Clinical Table (Tier 1, Tier 2, and Top Tier 3)
    prio_mask = tier_col.str.contains("Tier 1|Tier 2") | (pd.to_numeric(df.get("VARIANT_PRIORITY_SCORE", 0), errors="coerce").fillna(0) >= 30.0)
    prio_df = df[prio_mask].copy()
    prio_df = prio_df.sort_values("VARIANT_PRIORITY_SCORE", ascending=False)
    
    # Cap table to top 200 candidates to keep HTML responsive
    prio_df = prio_df.head(250)

    table_rows = []
    for _, r in prio_df.iterrows():
        locus = _safe_str(r.get("Locus", "-"), 22)
        hgvsc = _safe_str(r.get("HGVSc", "-"), 30)
        hgvsp = _safe_str(r.get("HGVSp", "-"), 25)
        tier_val = str(r.get("PRIORITY_TIER", "Tier 4"))
        
        tier_badge = "badge-tier4"
        if "Tier 1" in tier_val:
            tier_badge = "badge-tier1"
        elif "Tier 2" in tier_val:
            tier_badge = "badge-tier2"
        elif "Tier 3" in tier_val:
            tier_badge = "badge-tier3"

        score_val = _safe_float(r.get("VARIANT_PRIORITY_SCORE", 0), 1)
        impact_val = _safe_str(r.get("NEW_IMPACT", r.get("IMPACT", "-")), 15)
        conseq_val = _safe_str(r.get("Consequence", "-"), 25)
        af_val = _safe_float(r.get("gnomADv4_AF_grpmax_joint", r.get("AF", 0)), 6)
        
        splice_val = _safe_float(r.get("spliceai_custom_MAX", r.get("spliceAI_MAX", 0)), 2)
        spip_val = _safe_float(r.get("SPiP_prediction", r.get("SPiP", 0)), 2)
        am_val = _safe_float(r.get("am_pathogenicity", 0), 2)
        revel_val = _safe_float(r.get("REVEL_score", 0), 2)
        sv_stat = _safe_str(r.get("SpliceVault_status", "-"), 18)
        cln_val = _safe_str(r.get("CLNSIG", r.get("clinvar_clnsig", "-")), 22)

        # External Links
        loc_parts = str(r.get("Locus", "")).split("-")[0].split(":")
        if len(loc_parts) == 2:
            chr_num, pos_num = loc_parts[0], loc_parts[1]
            gnomad_link = f"https://gnomad.broadinstitute.org/variant/{chr_num}-{pos_num}-{r.get('REF','')}-{r.get('ALT','')}?dataset=gnomad_r4"
        else:
            gnomad_link = "#"

        table_rows.append(f"""
        <tr>
            <td><strong><a href="{gnomad_link}" target="_blank" style="color:#2980b9; text-decoration:none;">{locus}</a></strong></td>
            <td><code>{hgvsc}</code></td>
            <td><code>{hgvsp}</code></td>
            <td><span class="badge {tier_badge}">{tier_val.split(' (')[0]}</span></td>
            <td><strong>{score_val}</strong></td>
            <td><span class="impact-{impact_val.lower()}">{impact_val}</span></td>
            <td>{conseq_val}</td>
            <td>{af_val}</td>
            <td><strong>{splice_val}</strong></td>
            <td>{spip_val}</td>
            <td><strong>{am_val}</strong></td>
            <td>{revel_val}</td>
            <td>{sv_stat}</td>
            <td>{cln_val}</td>
        </tr>
        """)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Build Complete HTML Document
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clinical Prioritization Report: {gene_name} | GRCh38</title>
    <script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
    <style>
        :root {{
            --bg-color: #f4f7f9;
            --card-bg: #ffffff;
            --text-main: #2c3e50;
            --text-muted: #7f8c8d;
            --primary: #2980b9;
            --tier1: #e74c3c;
            --tier2: #e67e22;
            --tier3: #f1c40f;
            --tier4: #2ecc71;
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
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%);
            color: white;
            padding: 28px 32px;
            border-radius: 12px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            margin-bottom: 24px;
        }}
        .header-banner h1 {{
            margin: 0 0 8px 0;
            font-size: 28px;
            font-weight: 700;
        }}
        .meta-tags {{
            display: flex;
            gap: 12px;
            flex-wrap: wrap;
            margin-top: 12px;
        }}
        .meta-tag {{
            background: rgba(255,255,255,0.18);
            padding: 5px 12px;
            border-radius: 20px;
            font-size: 13px;
            font-weight: 500;
        }}
        .kpi-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }}
        .kpi-card {{
            background: var(--card-bg);
            padding: 20px;
            border-radius: 10px;
            border: 1px solid var(--border-color);
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
            text-align: center;
        }}
        .kpi-card.t1 {{ border-top: 5px solid var(--tier1); }}
        .kpi-card.t2 {{ border-top: 5px solid var(--tier2); }}
        .kpi-card.t3 {{ border-top: 5px solid var(--tier3); }}
        .kpi-card.t4 {{ border-top: 5px solid var(--tier4); }}
        .kpi-title {{
            font-size: 13px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            color: var(--text-muted);
            margin-bottom: 6px;
        }}
        .kpi-value {{
            font-size: 32px;
            font-weight: 700;
        }}
        .kpi-value.t1 {{ color: var(--tier1); }}
        .kpi-value.t2 {{ color: var(--tier2); }}
        .kpi-value.t3 {{ color: #d4ac0d; }}
        .kpi-value.t4 {{ color: var(--tier4); }}
        .charts-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(380px, 1fr));
            gap: 20px;
            margin-bottom: 28px;
        }}
        .chart-card {{
            background: var(--card-bg);
            padding: 16px;
            border-radius: 10px;
            border: 1px solid var(--border-color);
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
        }}
        .table-section {{
            background: var(--card-bg);
            padding: 24px;
            border-radius: 12px;
            border: 1px solid var(--border-color);
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
            margin-bottom: 30px;
        }}
        .table-header-flex {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 16px;
            flex-wrap: wrap;
            gap: 12px;
        }}
        .search-box {{
            padding: 8px 16px;
            border-radius: 20px;
            border: 1px solid #cbd5e0;
            font-size: 14px;
            width: 260px;
        }}
        .table-container {{
            overflow-x: auto;
            max-height: 650px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
            text-align: left;
        }}
        th {{
            background-color: #f8fafc;
            color: #475569;
            font-weight: 600;
            padding: 10px 12px;
            border-bottom: 2px solid var(--border-color);
            position: sticky;
            top: 0;
            z-index: 10;
        }}
        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #edf2f7;
            white-space: nowrap;
        }}
        tr:hover {{
            background-color: #f1f5f9;
        }}
        .badge {{
            display: inline-block;
            padding: 3px 8px;
            border-radius: 12px;
            font-size: 11px;
            font-weight: 700;
            color: white;
        }}
        .badge-tier1 {{ background-color: var(--tier1); }}
        .badge-tier2 {{ background-color: var(--tier2); }}
        .badge-tier3 {{ background-color: #f39c12; }}
        .badge-tier4 {{ background-color: var(--tier4); }}
        .impact-high {{ color: #c0392b; font-weight: 700; }}
        .impact-moderate {{ color: #d35400; font-weight: 600; }}
        .impact-low {{ color: #27ae60; }}
        .impact-modifier {{ color: #7f8c8d; }}
        .footer {{
            text-align: center;
            font-size: 12px;
            color: var(--text-muted);
            margin-top: 40px;
        }}
    </style>
</head>
<body>
    <div class="header-banner">
        <h1>🧬 Gene Clinical Prioritization Report: {gene_name}</h1>
        <div>Reference Transcript: <strong>{mane_tx}</strong> | Genome Build: <strong>GRCh38</strong> | Generated: <strong>{timestamp}</strong></div>
        <div class="meta-tags">
            <span class="meta-tag">Gene Constraint (pLI): <strong>{pli_val}</strong></span>
            <span class="meta-tag">ClinGen Priority: <strong>{gene_prio}</strong></span>
            <span class="meta-tag">Total Cohort Variants: <strong>{n_total:,}</strong></span>
            <span class="meta-tag">High Splicing Candidates: <strong>{n_splice_path:,}</strong></span>
            <span class="meta-tag">Consensus Missense Candidates: <strong>{n_missense_path:,}</strong></span>
        </div>
    </div>

    <!-- KPI Cards -->
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
            <small>Secondary research follow-up</small>
        </div>
        <div class="kpi-card t4">
            <div class="kpi-title">Tier 4: Benign / Filtered</div>
            <div class="kpi-value t4">{n_tier4:,}</div>
            <small>Common / non-damaging</small>
        </div>
    </div>

    <!-- Charts -->
    <div class="charts-grid">
        <div class="chart-card">
            {fig_donut.to_html(full_html=False, include_plotlyjs=False)}
        </div>
        <div class="chart-card">
            {fig_missense.to_html(full_html=False, include_plotlyjs=False)}
        </div>
        <div class="chart-card">
            {fig_splice.to_html(full_html=False, include_plotlyjs=False)}
        </div>
    </div>

    <!-- Clinical Triage Table -->
    <div class="table-section">
        <div class="table-header-flex">
            <div>
                <h2 style="margin:0 0 4px 0;">📋 Top Prioritized Candidate Variants ({len(prio_df)} displayed)</h2>
                <small style="color:var(--text-muted);">Includes Tier 1, Tier 2, and high-scoring Tier 3 candidates ranked by multi-omics priority score</small>
            </div>
            <input type="text" id="tableSearch" class="search-box" placeholder="🔍 Search locus, HGVSc, consequence..." onkeyup="filterTable()">
        </div>

        <div class="table-container">
            <table id="clinicalTable">
                <thead>
                    <tr>
                        <th>Locus (GRCh38)</th>
                        <th>HGVSc</th>
                        <th>HGVSp</th>
                        <th>Priority Tier</th>
                        <th>Score</th>
                        <th>NEW_IMPACT</th>
                        <th>Consequence</th>
                        <th>gnomAD AF</th>
                        <th>SpliceAI Δ</th>
                        <th>SPiP</th>
                        <th>AlphaMissense</th>
                        <th>REVEL</th>
                        <th>SpliceVault</th>
                        <th>ClinVar</th>
                    </tr>
                </thead>
                <tbody>
                    {''.join(table_rows)}
                </tbody>
            </table>
        </div>
    </div>

    <script>
        function filterTable() {{
            let input = document.getElementById("tableSearch");
            let filter = input.value.toLowerCase();
            let table = document.getElementById("clinicalTable");
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

    logger.info(f"Successfully generated medical clinical dashboard at: {out_html_path}")
    return True

def main():
    parser = argparse.ArgumentParser(description="Generate Medical Clinical Prioritization Report per Gene")
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
