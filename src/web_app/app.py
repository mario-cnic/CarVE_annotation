#!/usr/bin/env python3
"""
🧬 Universal Clinical Variant Prioritization & Interactive Reporting Web App

Features:
  - Drag-and-Drop file uploader for .xlsx, .parquet (.pq), .csv, .tsv, .vcf, .vcf.gz.
  - Automated ACMG/ClinGen tiering & priority scoring engine (build_newImpact, build_priority_tier).
  - Interactive KPI Metrics & Plotly Analytics charts.
  - Filterable, searchable variant table grid with export options.
  - Live embedded 4-Tab Clinical Prioritization Dashboard.
  - One-click Download Center (.html report, .xlsx table, .pq parquet).
"""

import os
import sys
import tempfile
import io
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as io_plotly
import streamlit as st
import streamlit.components.v1 as components

# Ensure pipeline and shared/utils modules are accessible
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PIPELINE_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, "..", ".."))
SHARED_UTILS = "/home/mruizp/data_lab_PGP/shared/utils/src"

for path in [SHARED_UTILS, os.path.join(PIPELINE_ROOT, "src", "python"), PIPELINE_ROOT]:
    if os.path.exists(path) and path not in sys.path:
        sys.path.insert(0, path)

try:
    from filter_variants import build_newImpact, build_priority_tier, parse_pangolin, parse_spip
    from generate_clinical_prioritization_report import generate_gene_report
except ImportError as e:
    st.error(f"Error importing pipeline modules: {e}")

# Streamlit Page Config
st.set_page_config(
    page_title="Clinical Variant Prioritization Hub",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for Premium Design
st.markdown("""
<style>
    .main-header {
        background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 100%);
        padding: 24px;
        border-radius: 12px;
        color: white;
        margin-bottom: 24px;
        box-shadow: 0 4px 15px rgba(0,0,0,0.1);
    }
    .main-header h1 {
        margin: 0;
        font-size: 28px;
        font-weight: 700;
    }
    .main-header p {
        margin: 6px 0 0 0;
        opacity: 0.85;
        font-size: 14px;
    }
    .metric-card {
        background: #ffffff;
        border-radius: 10px;
        padding: 16px;
        border: 1px solid #e2e8f0;
        box-shadow: 0 2px 8px rgba(0,0,0,0.04);
        text-align: center;
    }
    .metric-val {
        font-size: 26px;
        font-weight: 800;
        color: #0f172a;
    }
    .metric-lbl {
        font-size: 12px;
        font-weight: 600;
        color: #64748b;
        text-transform: uppercase;
    }
    .stDownloadButton > button {
        width: 100%;
        background-color: #2563eb;
        color: white;
        border-radius: 8px;
        font-weight: 600;
        padding: 10px;
    }
</style>
""", unsafe_allow_html=True)


def load_uploaded_file(uploaded_file):
    """Parse uploaded file into a Pandas DataFrame based on file extension."""
    filename = uploaded_file.name.lower()
    
    try:
        if filename.endswith(".parquet") or filename.endswith(".pq"):
            return pd.read_parquet(uploaded_file)
        elif filename.endswith(".xlsx") or filename.endswith(".xls"):
            return pd.read_excel(uploaded_file)
        elif filename.endswith(".csv"):
            return pd.read_csv(uploaded_file)
        elif filename.endswith(".tsv") or filename.endswith(".txt"):
            return pd.read_csv(uploaded_file, sep="\t")
        elif filename.endswith(".vcf") or filename.endswith(".vcf.gz"):
            # Save temporary file to parse with pysam / vcf_parser
            with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz" if filename.endswith(".gz") else ".vcf") as tmp:
                tmp.write(uploaded_file.getvalue())
                tmp_path = tmp.name
            
            try:
                from vcf_parser_pysam import parse_vcf_to_dataframe
                df = parse_vcf_to_dataframe(tmp_path)
            except Exception:
                # Fallback basic VCF header reader
                import gzip
                open_fn = gzip.open if filename.endswith(".gz") else open
                records = []
                with open_fn(tmp_path, "rt") as f:
                    for line in f:
                        if line.startswith("#"):
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) >= 8:
                            records.append({
                                "CHROM": parts[0],
                                "POS": parts[1],
                                "ID": parts[2],
                                "REF": parts[3],
                                "ALT": parts[4],
                                "QUAL": parts[5],
                                "FILTER": parts[6],
                                "INFO": parts[7],
                                "Locus": f"{parts[0]}:{parts[1]}-{parts[3]}-{parts[4]}"
                            })
                df = pd.DataFrame(records)
            finally:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            return df
        else:
            st.error(f"Unsupported file format: {filename}")
            return None
    except Exception as err:
        st.error(f"Failed to parse uploaded file '{uploaded_file.name}': {err}")
        return None


def run_tiering_and_scoring(df):
    """Applies build_newImpact and build_priority_tier to standardise variant classification."""
    # Standardize column mappings if needed
    col_mapping = {
        "gene": "SYMBOL", "Gene": "SYMBOL", "Gene_Name": "SYMBOL",
        "chrom": "CHROM", "pos": "POS", "ref": "REF", "alt": "ALT",
        "locus": "Locus", "LOCUS": "Locus", "hgvsc": "HGVSc", "hgvsp": "HGVSp"
    }
    for old_col, new_col in col_mapping.items():
        if old_col in df.columns and new_col not in df.columns:
            df[new_col] = df[old_col]
            
    if "Locus" not in df.columns and all(c in df.columns for c in ["CHROM", "POS", "REF", "ALT"]):
        df["Locus"] = df["CHROM"].astype(str) + ":" + df["POS"].astype(str) + "-" + df["REF"].astype(str) + "-" + df["ALT"].astype(str)
    elif "Locus" not in df.columns:
        df["Locus"] = [f"var_{i+1}" for i in range(len(df))]

    if "SYMBOL" not in df.columns:
        df["SYMBOL"] = "UNKNOWN"

    # Execute Tiering Engine
    try:
        df = build_newImpact(df)
        df = build_priority_tier(df)
    except Exception as e:
        st.warning(f"Note on tiering execution: {e}")

    return df


def main():
    # Header Banner
    st.markdown("""
    <div class="main-header">
        <h1>🧬 Universal Clinical Variant Prioritization Hub</h1>
        <p>Interactive ACMG/ClinGen tiering, functional splicing analysis, and clinical report generation for human genomic variants.</p>
    </div>
    """, unsafe_allow_html=True)

    # Sidebar Controls & File Upload
    st.sidebar.header("📁 Load Variant Data")
    uploaded_file = st.sidebar.file_uploader(
        "Upload Variant File (.xlsx, .parquet, .csv, .tsv, .vcf)",
        type=["xlsx", "xls", "parquet", "pq", "csv", "tsv", "txt", "vcf", "gz"],
        help="Upload any variant table or annotated VCF file."
    )

    # Demo File Option
    use_demo = st.sidebar.checkbox("Use Demo MYBPC3 Dataset", value=(uploaded_file is None))
    
    df_raw = None
    file_label = ""
    
    if uploaded_file is not None:
        file_label = uploaded_file.name
        df_raw = load_uploaded_file(uploaded_file)
    elif use_demo:
        demo_path = os.path.join(PIPELINE_ROOT, "test_data", "test_run", "results", "MYBPC3.parsed.clean.pq")
        if os.path.exists(demo_path):
            file_label = "MYBPC3 Demo (Parquet)"
            df_raw = pd.read_parquet(demo_path)
        else:
            st.sidebar.warning("Demo dataset file not found.")

    if df_raw is None:
        st.info("👋 Welcome! Please upload a variant file (.xlsx, .parquet, .csv, .vcf) in the sidebar to begin.")
        return

    st.sidebar.success(f"Loaded **{len(df_raw):,}** variants from `{file_label}`")

    # Run Prioritization Engine
    with st.spinner("Calculating ACMG/ClinGen priority tiers and continuous scores..."):
        df = run_tiering_and_scoring(df_raw.copy())

    # Main Tabs
    tab_summary, tab_grid, tab_report, tab_export = st.tabs([
        "📊 Executive Summary",
        "🔍 Interactive Grid",
        "📑 Clinical Dashboard",
        "📥 Export & Downloads"
    ])

    # ---------------- TAB 1: EXECUTIVE SUMMARY ----------------
    with tab_summary:
        st.subheader("📌 Key Clinical Metrics")
        
        tier_col = "PRIORITY_TIER" if "PRIORITY_TIER" in df.columns else None
        tier_counts = df[tier_col].value_counts().to_dict() if tier_col else {}
        
        c1, c2, c3, c4, c5 = st.columns(5)
        with c1:
            st.markdown(f'<div class="metric-card"><div class="metric-val">{len(df):,}</div><div class="metric-lbl">Total Variants</div></div>', unsafe_allow_html=True)
        with c2:
            t1 = tier_counts.get("Tier 1 (Critical Pathogenic Candidate)", 0)
            st.markdown(f'<div class="metric-card"><div class="metric-val" style="color:#dc2626;">{t1:,}</div><div class="metric-lbl">Tier 1 (Critical)</div></div>', unsafe_allow_html=True)
        with c3:
            t2 = tier_counts.get("Tier 2 (Likely Deleterious / Strong Candidate)", 0)
            st.markdown(f'<div class="metric-card"><div class="metric-val" style="color:#ea580c;">{t2:,}</div><div class="metric-lbl">Tier 2 (Strong)</div></div>', unsafe_allow_html=True)
        with c4:
            t3 = tier_counts.get("Tier 3 (VUS / Moderate Potential)", 0)
            st.markdown(f'<div class="metric-card"><div class="metric-val" style="color:#d97706;">{t3:,}</div><div class="metric-lbl">Tier 3 (VUS)</div></div>', unsafe_allow_html=True)
        with c5:
            t4 = tier_counts.get("Tier 4 (Benign / Tolerated)", 0)
            st.markdown(f'<div class="metric-card"><div class="metric-val" style="color:#16a34a;">{t4:,}</div><div class="metric-lbl">Tier 4 (Benign)</div></div>', unsafe_allow_html=True)

        st.markdown("---")
        
        # Interactive Analytics Charts
        g1, g2 = st.columns(2)
        with g1:
            if tier_col and not df[tier_col].empty:
                t_df = df[tier_col].value_counts().reset_index()
                t_df.columns = ["Tier", "Count"]
                color_map = {
                    "Tier 1 (Critical Pathogenic Candidate)": "#dc2626",
                    "Tier 2 (Likely Deleterious / Strong Candidate)": "#ea580c",
                    "Tier 3 (VUS / Moderate Potential)": "#d97706",
                    "Tier 4 (Benign / Tolerated)": "#16a34a"
                }
                fig1 = px.pie(t_df, values="Count", names="Tier", title="Variant Stratification by Priority Tier",
                              color="Tier", color_discrete_map=color_map, hole=0.4)
                st.plotly_chart(fig1, use_container_width=True)
        
        with g2:
            imp_col = "NEW_IMPACT" if "NEW_IMPACT" in df.columns else ("IMPACT" if "IMPACT" in df.columns else None)
            if imp_col and not df[imp_col].empty:
                i_df = df[imp_col].value_counts().reset_index()
                i_df.columns = ["Impact", "Count"]
                fig2 = px.bar(i_df, x="Impact", y="Count", title="Functional Impact Spectrum",
                              color="Impact", color_discrete_sequence=px.colors.qualitative.Set2)
                st.plotly_chart(fig2, use_container_width=True)

    # ---------------- TAB 2: INTERACTIVE GRID ----------------
    with tab_grid:
        st.subheader("🔍 Filter & Search Variant Table")
        
        f1, f2, f3 = st.columns([2, 2, 3])
        with f1:
            if "SYMBOL" in df.columns:
                selected_genes = st.multiselect("Filter Gene(s):", options=sorted(df["SYMBOL"].dropna().unique()))
            else:
                selected_genes = []
        with f2:
            if tier_col:
                selected_tiers = st.multiselect("Filter Priority Tier:", options=list(df[tier_col].unique()))
            else:
                selected_tiers = []
        with f3:
            search_query = st.text_input("Search Locus / HGVSc / HGVSp / Consequence:", value="")

        df_filtered = df.copy()
        if selected_genes:
            df_filtered = df_filtered[df_filtered["SYMBOL"].isin(selected_genes)]
        if selected_tiers:
            df_filtered = df_filtered[df_filtered[tier_col].isin(selected_tiers)]
        if search_query:
            query_lower = search_query.lower()
            match_mask = df_filtered.astype(str).apply(lambda row: row.str.lower().str.contains(query_lower).any(), axis=1)
            df_filtered = df_filtered[match_mask]

        st.caption(f"Showing **{len(df_filtered):,}** of **{len(df):,}** variants")
        st.dataframe(df_filtered, use_container_width=True, height=450)

    # ---------------- TAB 3: CLINICAL DASHBOARD ----------------
    with tab_report:
        st.subheader("📑 4-Tab Clinical Prioritization Dashboard")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_pq = os.path.join(tmpdir, "temp_data.pq")
            tmp_html = os.path.join(tmpdir, "clinical_report.html")
            
            df.to_parquet(tmp_pq)
            try:
                generate_gene_report(tmp_pq, tmp_html)
                if os.path.exists(tmp_html):
                    with open(tmp_html, "r", encoding="utf-8") as hf:
                        html_content = hf.read()
                    
                    st.components.v1.html(html_content, height=850, scrolling=True)
                else:
                    st.error("Report HTML failed to generate.")
            except Exception as ex:
                st.error(f"Error rendering clinical report: {ex}")

    # ---------------- TAB 4: EXPORT & DOWNLOADS ----------------
    with tab_export:
        st.subheader("📥 Export & Download Center")
        st.write("Export your analyzed, prioritized variant dataset in your preferred format:")
        
        exp1, exp2, exp3 = st.columns(3)
        
        # 1. HTML Report Download
        with exp1:
            st.markdown("#### 🌐 Interactive HTML Dashboard")
            st.caption("Self-contained interactive clinical dashboard report with embedded Plotly charts.")
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_pq = os.path.join(tmpdir, "temp_data.pq")
                tmp_html = os.path.join(tmpdir, "clinical_report.html")
                df.to_parquet(tmp_pq)
                generate_gene_report(tmp_pq, tmp_html)
                if os.path.exists(tmp_html):
                    with open(tmp_html, "rb") as f:
                        html_bytes = f.read()
                    st.download_button(
                        label="Download Dashboard (.html)",
                        data=html_bytes,
                        file_name=f"clinical_prioritization_report_{file_label.replace(' ', '_')}.html",
                        mime="text/html"
                    )

        # 2. Excel Table Download
        with exp2:
            st.markdown("#### 📊 Excel Table Report")
            st.caption("Structured Excel workbook formatted with prioritized variant columns.")
            excel_buffer = io.BytesIO()
            with pd.ExcelWriter(excel_buffer, engine="openpyxl") as writer:
                df.to_excel(writer, index=False, sheet_name="Prioritized Variants")
            excel_bytes = excel_buffer.getvalue()
            st.download_button(
                label="Download Table (.xlsx)",
                data=excel_bytes,
                file_name=f"prioritized_variants_{file_label.replace(' ', '_')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

        # 3. Parquet Data Download
        with exp3:
            st.markdown("#### 📦 Parquet Binary Dataset")
            st.caption("High-performance compressed columnar Parquet binary format (.pq).")
            pq_buffer = io.BytesIO()
            df.to_parquet(pq_buffer, index=False)
            pq_bytes = pq_buffer.getvalue()
            st.download_button(
                label="Download Parquet (.pq)",
                data=pq_bytes,
                file_name=f"prioritized_variants_{file_label.replace(' ', '_')}.pq",
                mime="application/octet-stream"
            )


if __name__ == "__main__":
    main()
