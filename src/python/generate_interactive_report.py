#!/usr/bin/env python3
"""
Interactive HTML Report Generator for Variant Annotation Pipeline

Generates a standalone, interactive HTML dashboard containing:
  - Interactive Consequence Breakdown Pie & Bar Charts
  - Interactive 2D Pathogenicity Scatter Plot (REVEL vs AlphaMissense vs SPiP)
  - Interactive Allele Frequency Spectrum & Log-scale distribution
  - Interactive summary table
"""

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["ARROW_IO_THREADS"] = "1"

import sys
import argparse
import logging
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("generate_interactive_report")

def parse_args():
    parser = argparse.ArgumentParser(description="Generate Interactive HTML Variant Dashboard")
    parser.add_argument("--input", required=True, help="Input .pq, .tsv, or .csv file")
    parser.add_argument("--output", required=True, help="Output HTML file path")
    parser.add_argument("--title", default="Variant Annotation Pipeline Dashboard", help="Report title")
    return parser.parse_args()

def main():
    args = parse_args()
    logger.info(f"Loading input file: {args.input}...")
    
    if args.input.endswith('.pq') or args.input.endswith('.parquet'):
        df = pd.read_parquet(args.input)
    elif args.input.endswith('.tsv'):
        df = pd.read_csv(args.input, sep='\t')
    else:
        df = pd.read_csv(args.input)

    try:
        import plotly.express as px
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        logger.warning("Plotly not installed. Generating HTML fallback summary report...")
        generate_html_fallback(df, args.output, args.title)
        return

    logger.info("Building Plotly interactive charts...")
    fig_list = []

    # 1. Consequence Bar Chart
    cons_col = intersect_col(['Consequence', 'consequence', 'VEP_Consequence'], df)
    if cons_col:
        cons_counts = df[cons_col].value_counts().reset_index()
        cons_counts.columns = ['Consequence', 'Count']
        fig_cons = px.bar(
            cons_counts.head(20),
            x='Count',
            y='Consequence',
            orientation='h',
            title='Top Variant Consequences',
            color='Count',
            color_continuous_scale='Viridis'
        )
        fig_cons.update_layout(yaxis={'categoryorder': 'total ascending'})
        fig_list.append(("Variant Consequences", fig_cons.to_html(full_html=False, include_plotlyjs='cdn')))

    # 2. Pathogenicity Scatter Plot (REVEL vs AlphaMissense)
    revel_col = intersect_col(['REVEL_score', 'REVEL', 'revel'], df)
    am_col = intersect_col(['am_pathogenicity', 'AlphaMissense', 'alphamissense'], df)
    spip_col = intersect_col(['SPiP_score', 'SPiP', 'spip_score'], df)

    if revel_col and am_col:
        plot_df = df.copy()
        plot_df['REVEL'] = pd.to_numeric(plot_df[revel_col], errors='coerce')
        plot_df['AlphaMissense'] = pd.to_numeric(plot_df[am_col], errors='coerce')
        if spip_col:
            plot_df['SPiP'] = pd.to_numeric(plot_df[spip_col], errors='coerce').fillna(0.0)
            color_arg = 'SPiP'
        else:
            color_arg = None
        
        plot_df = plot_df.dropna(subset=['REVEL', 'AlphaMissense'])
        if len(plot_df) > 0:
            fig_patho = px.scatter(
                plot_df,
                x='REVEL',
                y='AlphaMissense',
                color=color_arg,
                hover_data=[c for c in ['Gene', 'cDNA', 'HGVSc', 'Consequence'] if c in plot_df.columns],
                title='Pathogenicity Score Correlation (REVEL vs AlphaMissense)',
                color_continuous_scale='Portland'
            )
            fig_patho.add_vline(x=0.75, line_dash="dash", line_color="red")
            fig_patho.add_hline(y=0.80, line_dash="dash", line_color="red")
            fig_list.append(("Pathogenicity Correlation", fig_patho.to_html(full_html=False, include_plotlyjs=False)))

    # 3. AF Spectrum
    af_col = intersect_col(['AF', 'gnomAD_AF', 'gnomAD_AF_joint', 'AF_joint'], df)
    if af_col:
        af_df = df.copy()
        af_df['AF_num'] = pd.to_numeric(af_df[af_col], errors='coerce')
        af_df = af_df[af_df['AF_num'] > 0]
        if len(af_df) > 0:
            fig_af = px.histogram(
                af_df,
                x='AF_num',
                log_x=True,
                nbins=40,
                title='gnomAD Allele Frequency Spectrum (Log Scale)',
                color_discrete_sequence=['#2ecc71']
            )
            fig_list.append(("Allele Frequency Spectrum", fig_af.to_html(full_html=False, include_plotlyjs=False)))

    # Assemble HTML document
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{args.title}</title>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 30px; background-color: #f8f9fa; color: #2c3e50; }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        .metric-card {{ background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); margin-bottom: 25px; }}
        .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-bottom: 30px; }}
        .card {{ background: white; padding: 15px; border-radius: 8px; text-align: center; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }}
        .card-num {{ font-size: 28px; font-weight: bold; color: #3498db; }}
    </style>
</head>
<body>
    <h1>📊 {args.title}</h1>
    <div class="grid">
        <div class="card"><div>Total Variants</div><div class="card-num">{len(df)}</div></div>
        <div class="card"><div>Unique Genes</div><div class="card-num">{df['Gene'].nunique() if 'Gene' in df.columns else 'N/A'}</div></div>
        <div class="card"><div>Columns</div><div class="card-num">{len(df.columns)}</div></div>
    </div>
"""

    for title, div in fig_list:
        html_content += f"""
    <div class="metric-card">
        <h2>{title}</h2>
        {div}
    </div>
"""

    html_content += """
</body>
</html>
"""

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output, "w") as f:
        f.write(html_content)

    logger.info(f"Successfully generated interactive HTML report at: {args.output}")

def intersect_col(candidates, df):
    for c in candidates:
        if c in df.columns:
            return c
    return None

def generate_html_fallback(df, output_path, title):
    html = f"""<!DOCTYPE html>
<html>
<head><title>{title}</title></head>
<body>
    <h1>{title}</h1>
    <p>Total Variants: {len(df)}</p>
    <h3>Columns:</h3>
    <ul>{''.join(f'<li>{c}</li>' for c in df.columns)}</ul>
</body>
</html>
"""
    with open(output_path, "w") as f:
        f.write(html)

if __name__ == "__main__":
    main()
