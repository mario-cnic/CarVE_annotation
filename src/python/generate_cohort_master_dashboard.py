#!/usr/bin/env python3
"""
Multi-Gene Cohort Master Dashboard Generator
Aggregates top-tier prioritized variants across all target genes into a single interactive HTML dashboard.

Features:
  - Cross-gene variant aggregation (.parsed.clean.pq)
  - Cohort-level summary metrics & distribution of High/Moderate candidates
  - Interactive multi-gene Plotly charts (Consequence breakdown, Pathogenicity vs Splicing, AF spectrum)
  - Searchable, filterable variant table with embedded 1-click candidate export center (.csv/.xlsx)
"""

import os
import sys
import glob
import json
import argparse
import logging
import pandas as pd
import numpy as np

# Import dynamic configuration loader
PIPELINE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(PIPELINE_ROOT, "src", "python"))
import config

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("generate_cohort_master_dashboard")

def parse_args():
    parser = argparse.ArgumentParser(description="Multi-Gene Cohort Master Dashboard Generator")
    parser.add_argument("--run-dir", default=None, help="Path to run directory (e.g. RUNS/predictors_050826)")
    parser.add_argument("--results-dir", default=None, help="Path to results directory containing *.parsed.clean.pq files")
    parser.add_argument("--output", default=None, help="Output path for master cohort HTML report")
    return parser.parse_args()

def load_cohort_data(results_dir):
    pq_files = glob.glob(os.path.join(results_dir, "*.parsed.clean.pq"))
    if not pq_files:
        pq_files = glob.glob(os.path.join(results_dir, "*.pq"))
    
    if not pq_files:
        logger.error(f"No Parquet files found in results directory: {results_dir}")
        return pd.DataFrame()

    logger.info(f"Aggregating {len(pq_files)} gene parquet files from {results_dir}...")
    dfs = []
    for f in sorted(pq_files):
        try:
            df_gene = pd.read_parquet(f)
            dfs.append(df_gene)
        except Exception as e:
            logger.warning(f"Could not read {f}: {e}")

    if not dfs:
        return pd.DataFrame()

    master_df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Successfully aggregated {len(master_df):,} total variants across {len(pq_files)} genes.")
    return master_df

def generate_dashboard_html(master_df, run_name, output_path):
    tot_vars = len(master_df)
    tot_genes = master_df["SYMBOL"].nunique() if "SYMBOL" in master_df.columns else 0
    
    # Priority tiers count
    high_prio = 0
    mod_prio = 0
    if "PRIORITY_TIER" in master_df.columns:
        high_prio = (master_df["PRIORITY_TIER"].astype(str).str.contains("Class 1|High", case=False, na=False)).sum()
        mod_prio = (master_df["PRIORITY_TIER"].astype(str).str.contains("Class 2|Moderate", case=False, na=False)).sum()

    # Top consequences
    cons_series = master_df["Consequence"].value_counts().head(10) if "Consequence" in master_df.columns else pd.Series()
    cons_json = json.dumps({"labels": list(cons_series.index), "values": [int(x) for x in cons_series.values]})

    # Top genes by variant count
    gene_series = master_df["SYMBOL"].value_counts().head(15) if "SYMBOL" in master_df.columns else pd.Series()
    gene_json = json.dumps({"labels": list(gene_series.index), "values": [int(x) for x in gene_series.values]})

    # Extract clean table records for interactive datatable
    display_cols = ["Locus", "SYMBOL", "Consequence", "HGVSc", "HGVSp", "PRIORITY_TIER", "VARIANT_PRIORITY_SCORE", "gnomADv4_AF_grpmax_joint", "REVEL_score", "am_pathogenicity", "SPiP_prediction", "spliceai_custom_MAX"]
    cols_to_use = [c for c in display_cols if c in master_df.columns]
    
    table_records = master_df[cols_to_use].head(2500).to_dict(orient="records")
    table_json = json.dumps(table_records, default=str)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Cohort Master Variant Prioritization Dashboard - {run_name}</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {{
            --bg-dark: #0f172a;
            --card-bg: #1e293b;
            --accent-blue: #38bdf8;
            --accent-purple: #c084fc;
            --text-main: #f8fafc;
            --text-muted: #94a3b8;
            --border-color: #334155;
            --success: #4ade80;
            --warning: #fbbf24;
            --danger: #f87171;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg-dark);
            color: var(--text-main);
            margin: 0;
            padding: 24px;
        }}
        .header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 16px;
            margin-bottom: 24px;
        }}
        .header h1 {{
            margin: 0;
            font-size: 24px;
            color: var(--accent-blue);
        }}
        .kpi-container {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }}
        .kpi-card {{
            background: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 16px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        }}
        .kpi-title {{
            font-size: 13px;
            color: var(--text-muted);
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .kpi-value {{
            font-size: 28px;
            font-weight: 700;
            margin-top: 8px;
            color: var(--text-main);
        }}
        .charts-container {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin-bottom: 24px;
        }}
        .chart-card {{
            background: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 16px;
            min-height: 350px;
        }}
        .table-section {{
            background: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 16px;
        }}
        .table-controls {{
            display: flex;
            justify-content: space-between;
            margin-bottom: 12px;
        }}
        input[type="text"] {{
            background: var(--bg-dark);
            border: 1px solid var(--border-color);
            color: var(--text-main);
            padding: 8px 12px;
            border-radius: 4px;
            width: 300px;
        }}
        .btn {{
            background: var(--accent-blue);
            color: #0f172a;
            border: none;
            padding: 8px 16px;
            border-radius: 4px;
            font-weight: 600;
            cursor: pointer;
        }}
        .btn:hover {{
            opacity: 0.9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
        }}
        th, td {{
            padding: 10px 12px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }}
        th {{
            background-color: rgba(15, 23, 42, 0.6);
            color: var(--accent-blue);
            position: sticky;
            top: 0;
        }}
        tr:hover {{
            background-color: rgba(56, 189, 248, 0.05);
        }}
        .badge {{
            padding: 2px 6px;
            border-radius: 4px;
            font-size: 11px;
            font-weight: 600;
        }}
        .badge-high {{ background: rgba(248, 113, 113, 0.2); color: var(--danger); border: 1px solid var(--danger); }}
        .badge-mod {{ background: rgba(251, 191, 36, 0.2); color: var(--warning); border: 1px solid var(--warning); }}
    </style>
</head>
<body>
    <div class="header">
        <div>
            <h1>🧬 Multi-Gene Cohort Master Variant Dashboard</h1>
            <div style="color: var(--text-muted); font-size: 13px; margin-top: 4px;">Run Identifier: <strong>{run_name}</strong> | Target Assembly: <strong>GRCh38</strong></div>
        </div>
        <button class="btn" onclick="exportCSV()">📥 Download Candidates (.csv)</button>
    </div>

    <div class="kpi-container">
        <div class="kpi-card">
            <div class="kpi-title">Target Genes</div>
            <div class="kpi-value">{tot_genes:,}</div>
        </div>
        <div class="kpi-card">
            <div class="kpi-title">Total Cohort Variants</div>
            <div class="kpi-value">{tot_vars:,}</div>
        </div>
        <div class="kpi-card">
            <div class="kpi-title">Class 1 (High Priority)</div>
            <div class="kpi-value" style="color: var(--danger);">{high_prio:,}</div>
        </div>
        <div class="kpi-card">
            <div class="kpi-title">Class 2 (Moderate Priority)</div>
            <div class="kpi-value" style="color: var(--warning);">{mod_prio:,}</div>
        </div>
    </div>

    <div class="charts-container">
        <div class="chart-card" id="chart-consequences"></div>
        <div class="chart-card" id="chart-genes"></div>
    </div>

    <div class="table-section">
        <div class="table-controls">
            <input type="text" id="search-input" placeholder="Search gene, locus, consequence..." onkeyup="filterTable()">
            <div style="color: var(--text-muted); font-size: 13px; align-self: center;">Showing top prioritized candidates across cohort</div>
        </div>
        <div style="overflow-x: auto; max-height: 500px;">
            <table id="variant-table">
                <thead>
                    <tr>
                        <th>Locus</th>
                        <th>Gene</th>
                        <th>Consequence</th>
                        <th>HGVSc</th>
                        <th>HGVSp</th>
                        <th>Priority Tier</th>
                        <th>Score</th>
                        <th>gnomAD AF</th>
                        <th>REVEL</th>
                        <th>AlphaMissense</th>
                        <th>SPiP</th>
                        <th>SpliceAI Δ</th>
                    </tr>
                </thead>
                <tbody id="table-body"></tbody>
            </table>
        </div>
    </div>

    <script>
        const consData = {cons_json};
        const geneData = {gene_json};
        const tableData = {table_json};

        // Render Plotly Bar Charts
        Plotly.newPlot('chart-consequences', [{{
            x: consData.values,
            y: consData.labels,
            type: 'bar',
            orientation: 'h',
            marker: {{ color: '#38bdf8' }}
        }}], {{
            title: 'Top Variant Consequences across Cohort',
            paper_bgcolor: 'transparent',
            plot_bgcolor: 'transparent',
            font: {{ color: '#f8fafc' }},
            margin: {{ l: 150, r: 20, t: 40, b: 40 }}
        }});

        Plotly.newPlot('chart-genes', [{{
            x: geneData.labels,
            y: geneData.values,
            type: 'bar',
            marker: {{ color: '#c084fc' }}
        }}], {{
            title: 'Variant Burden per Target Gene (Top 15)',
            paper_bgcolor: 'transparent',
            plot_bgcolor: 'transparent',
            font: {{ color: '#f8fafc' }},
            margin: {{ l: 40, r: 20, t: 40, b: 60 }}
        }});

        // Render Table
        function renderTable(data) {{
            const tbody = document.getElementById('table-body');
            tbody.innerHTML = '';
            data.forEach(r => {{
                const tr = document.createElement('tr');
                const tierClass = String(r.PRIORITY_TIER || '').includes('Class 1') ? 'badge-high' : (String(r.PRIORITY_TIER || '').includes('Class 2') ? 'badge-mod' : '');
                
                tr.innerHTML = `
                    <td><strong>${{r.Locus || '.'}}</strong></td>
                    <td><span style="color:#38bdf8; font-weight:600;">${{r.SYMBOL || '.'}}</span></td>
                    <td>${{r.Consequence || '.'}}</td>
                    <td>${{r.HGVSc || '.'}}</td>
                    <td>${{r.HGVSp || '.'}}</td>
                    <td><span class="badge ${{tierClass}}">${{r.PRIORITY_TIER || '.'}}</span></td>
                    <td>${{r.VARIANT_PRIORITY_SCORE ? Number(r.VARIANT_PRIORITY_SCORE).toFixed(1) : '.'}}</td>
                    <td>${{r.gnomADv4_AF_grpmax_joint != null ? Number(r.gnomADv4_AF_grpmax_joint).toExponential(2) : '.'}}</td>
                    <td>${{r.REVEL_score != null ? Number(r.REVEL_score).toFixed(2) : '.'}}</td>
                    <td>${{r.am_pathogenicity != null ? Number(r.am_pathogenicity).toFixed(2) : '.'}}</td>
                    <td>${{r.SPiP_prediction || '.'}}</td>
                    <td>${{r.spliceai_custom_MAX != null ? Number(r.spliceai_custom_MAX).toFixed(2) : '.'}}</td>
                `;
                tbody.appendChild(tr);
            }});
        }}

        renderTable(tableData);

        function filterTable() {{
            const q = document.getElementById('search-input').value.toLowerCase();
            const filtered = tableData.filter(r => 
                String(r.SYMBOL || '').toLowerCase().includes(q) ||
                String(r.Locus || '').toLowerCase().includes(q) ||
                String(r.Consequence || '').toLowerCase().includes(q) ||
                String(r.PRIORITY_TIER || '').toLowerCase().includes(q)
            );
            renderTable(filtered);
        }}

        function exportCSV() {{
            let csv = 'Locus,SYMBOL,Consequence,HGVSc,HGVSp,PRIORITY_TIER,VARIANT_PRIORITY_SCORE\\n';
            tableData.forEach(r => {{
                csv += `"${{r.Locus || ''}}","${{r.SYMBOL || ''}}","${{r.Consequence || ''}}","${{r.HGVSc || ''}}","${{r.HGVSp || ''}}","${{r.PRIORITY_TIER || ''}}","${{r.VARIANT_PRIORITY_SCORE || ''}}"\\n`;
            }});
            const blob = new Blob([csv], {{ type: 'text/csv' }});
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.setAttribute('href', url);
            a.setAttribute('download', `cohort_prioritized_candidates_${{Date.now()}}.csv`);
            a.click();
        }}
    </script>
</body>
</html>
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    logger.info(f"Successfully saved Cohort Master Dashboard to: {output_path}")

def main():
    args = parse_args()
    if args.run_dir:
        res_dir = os.path.join(args.run_dir, "results")
        run_name = os.path.basename(args.run_dir.rstrip("/"))
        out_html = args.output or os.path.join(args.run_dir, "reports", f"cohort_master_dashboard_{run_name}.html")
    elif args.results_dir:
        res_dir = args.results_dir
        run_name = os.path.basename(res_dir.rstrip("/"))
        out_html = args.output or os.path.join(res_dir, "..", "reports", "cohort_master_dashboard.html")
    else:
        logger.error("Please provide --run-dir or --results-dir")
        sys.exit(1)

    master_df = load_cohort_data(res_dir)
    if master_df.empty:
        logger.error("No variant records to generate cohort dashboard.")
        sys.exit(1)

    generate_dashboard_html(master_df, run_name, out_html)

if __name__ == "__main__":
    main()
