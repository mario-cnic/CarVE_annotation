#!/usr/bin/env python3
"""
Comprehensive Run Auditor & Quality Control Report Generator

Audits pipeline run output directories:
  1. Checks input vs output variant counts per gene.
  2. Verifies completion status for all 7 pipeline steps.
  3. Audits column count, schema consistency, and annotation completion rates.
  4. Checks log files (_log/ folder) for execution errors or warnings.
  5. Generates a detailed Markdown audit report.
"""

import os
import sys
import glob
import argparse
import logging
import pandas as pd
import numpy as np
from datetime import datetime

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("audit_run_results")

def parse_args():
    parser = argparse.ArgumentParser(description="Audit pipeline run output results")
    parser.add_argument("--run-dir", default="RUNS/predictors_050826", help="Path to run directory")
    parser.add_argument("--raw-dir", default="RUNS/predictors_050826/input/by_gene", help="Path to input raw folder")
    parser.add_argument("--output-report", default=None, help="Path for output markdown report file")
    return parser.parse_args()

def audit_run(run_dir, raw_dir, output_report_path=None):
    run_name = os.path.basename(os.path.normpath(run_dir))
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    tmp_dir = os.path.join(run_dir, "_tmp")
    ann_dir = os.path.join(run_dir, "annotation")
    res_dir = os.path.join(run_dir, "results")
    log_dir = os.path.join(run_dir, "_log")
    rep_dir = os.path.join(run_dir, "reports")
    flt_dir = os.path.join(run_dir, "filtered")

    if output_report_path is None:
        output_report_path = os.path.join(rep_dir, f"audit_report_{run_name}.md")

    # 1. Identify input files
    input_files = sorted(glob.glob(os.path.join(raw_dir, "*.pq")) + glob.glob(os.path.join(raw_dir, "*.parquet")) + glob.glob(os.path.join(raw_dir, "*.tsv")) + glob.glob(os.path.join(raw_dir, "*.csv")))

    if not input_files:
        logger.error(f"No input files found in {raw_dir}")
        return

    logger.info(f"Auditing run '{run_name}' with {len(input_files)} input genes...")

    results = []
    all_gene_columns = {}
    schema_set = None
    column_inconsistencies = []

    key_annotation_cols = [
        'Consequence', 'REVEL_score', 'am_pathogenicity', 'SPiP', 'SPiP_interpretation', 'CADD_PHRED',
        'splicevardb', 'SpliceVault_top_events', 'AF_joint', 'gnomAD_AF_joint'
    ]

    for inp_file in input_files:
        base_name = os.path.basename(inp_file)
        gene_name = base_name.split('.')[0].split('_')[0]

        # Read input counts
        try:
            if inp_file.endswith('.pq') or inp_file.endswith('.parquet'):
                df_inp = pd.read_parquet(inp_file)
            else:
                sep = '\t' if inp_file.endswith('.tsv') else ','
                df_inp = pd.read_csv(inp_file, sep=sep)
            n_input_vars = len(df_inp)
            inp_status = "VALID"
        except Exception as e:
            n_input_vars = 0
            inp_status = f"CORRUPTED ({e})"

        # Check step files
        vcf_tmp = os.path.join(tmp_dir, f"{gene_name}.vcf.gz")
        spip_vcf = os.path.join(ann_dir, f"{gene_name}.annSPiP.vcf.gz")
        if not os.path.exists(spip_vcf):
            spip_vcf = os.path.join(ann_dir, f"{gene_name}.annSPiP.vcf")
        vep_vcf = os.path.join(ann_dir, f"{gene_name}.annVEP.vcf.gz")
        annotated_vcf = os.path.join(ann_dir, f"{gene_name}.annotated.vcf.gz")
        parsed_tsv = os.path.join(res_dir, f"{gene_name}.parsed.tsv")
        final_pq = os.path.join(res_dir, f"{gene_name}.parsed.clean.pq")

        step_status = {
            "vcf_conv": os.path.exists(vcf_tmp),
            "spip": os.path.exists(spip_vcf) and os.path.getsize(spip_vcf) > 0,
            "vep": os.path.exists(vep_vcf) and os.path.getsize(vep_vcf) > 0,
            "mvep": os.path.exists(annotated_vcf) and os.path.getsize(annotated_vcf) > 0,
            "vcf2tsv": os.path.exists(parsed_tsv) and os.path.getsize(parsed_tsv) > 0,
            "final_pq": os.path.exists(final_pq) and os.path.getsize(final_pq) > 0
        }

        all_steps_ok = all(step_status.values())

        # Audit final pq file
        n_out_vars = 0
        n_cols = 0
        pq_status = "MISSING"
        ann_stats = {}

        if step_status["final_pq"]:
            try:
                df_out = pd.read_parquet(final_pq)
                n_out_vars = len(df_out)
                n_cols = len(df_out.columns)
                pq_status = "VALID" if n_out_vars > 0 else "EMPTY"
                all_gene_columns[gene_name] = list(df_out.columns)

                if schema_set is None:
                    schema_set = set(df_out.columns)
                else:
                    diff = schema_set.symmetric_difference(set(df_out.columns))
                    if diff:
                        column_inconsistencies.append((gene_name, list(diff)))

                # Check annotation completion percentages
                for col in key_annotation_cols:
                    found_col = None
                    for c in [col, col.lower(), col.upper()]:
                        if c in df_out.columns:
                            found_col = c
                            break
                    if found_col:
                        non_null = df_out[found_col].notna() & (df_out[found_col].astype(str).str.strip() != '') & (df_out[found_col].astype(str) != 'nan') & (df_out[found_col].astype(str) != '-')
                        pct = (non_null.sum() / max(1, n_out_vars)) * 100
                        ann_stats[col] = f"{pct:.1f}%"
                    else:
                        ann_stats[col] = "N/A"
            except Exception as e:
                pq_status = f"CORRUPTED ({e})"

        # Check logs for errors
        gene_log_dir = os.path.join(log_dir, gene_name)
        log_errors = []
        if os.path.exists(gene_log_dir):
            for err_file in glob.glob(os.path.join(gene_log_dir, "*.err")):
                if os.path.getsize(err_file) > 0:
                    with open(err_file, 'r', errors='ignore') as ef:
                        content = ef.read()
                        if "error" in content.lower() or "exception" in content.lower() or "command not found" in content.lower() or "killed" in content.lower():
                            log_errors.append(os.path.basename(err_file))

        results.append({
            "Gene": gene_name,
            "Input_Vars": n_input_vars,
            "Output_Vars": n_out_vars,
            "Status": "FINISHED" if all_steps_ok and pq_status == "VALID" else ("PENDING" if not os.path.exists(final_pq) else "ERROR"),
            "Columns": n_cols,
            "Final_PQ": pq_status,
            "Steps_OK": f"{sum(step_status.values())}/{len(step_status)}",
            "REVEL_%": ann_stats.get("REVEL_score", "N/A"),
            "AlphaMissense_%": ann_stats.get("am_pathogenicity", "N/A"),
            "SPiP_%": ann_stats.get("SPiP", "N/A"),
            "SPiP_Interp_%": ann_stats.get("SPiP_interpretation", "N/A"),
            "SpliceVarDB_%": ann_stats.get("splicevardb", "N/A"),
            "SpliceVault_%": ann_stats.get("SpliceVault_top_events", "N/A"),
            "Log_Issues": ", ".join(log_errors) if log_errors else "None"
        })

    df_results = pd.DataFrame(results)

    # Build Markdown Report
    n_total = len(df_results)
    n_finished = len(df_results[df_results["Status"] == "FINISHED"])
    n_pending = len(df_results[df_results["Status"] == "PENDING"])
    n_error = len(df_results[df_results["Status"] == "ERROR"])

    total_in_vars = df_results["Input_Vars"].sum()
    total_out_vars = df_results["Output_Vars"].sum()

    md = []
    md.append(f"# 🧪 Pipeline Execution Audit Report: `{run_name}`")
    md.append(f"**Audit Timestamp**: `{timestamp}`  ")
    md.append(f"**Target Build**: `GRCh38`  ")
    md.append(f"**Input Directory**: `{raw_dir}`  ")
    md.append(f"**Results Directory**: `{res_dir}`  \n")

    md.append("---")
    md.append("## 📊 Executive Summary")
    md.append(f"- **Total Genes Input**: `{n_total}`")
    md.append(f"- **Genes Fully Finished**: `{n_finished} / {n_total}` ({(n_finished/n_total)*100:.1f}%)")
    md.append(f"- **Genes Pending**: `{n_pending}`")
    md.append(f"- **Genes with Errors**: `{n_error}`")
    md.append(f"- **Total Input Variants**: `{total_in_vars:,}`")
    md.append(f"- **Total Annotated & Clean Output Variants**: `{total_out_vars:,}`")
    md.append(f"- **Variant Retention Rate**: `{((total_out_vars/total_in_vars)*100):.2f}%`" if total_in_vars > 0 else "N/A")
    md.append("")

    md.append("---")
    md.append("## 📋 Per-Gene Audit Metrics")
    md.append(df_results[["Gene", "Status", "Input_Vars", "Output_Vars", "Columns", "Steps_OK", "Final_PQ", "Log_Issues"]].to_markdown(index=False))
    md.append("")

    md.append("---")
    md.append("## 🧬 Key Annotation Completeness (% Non-Null)")
    md.append(df_results[["Gene", "REVEL_%", "AlphaMissense_%", "SPiP_%", "SPiP_Interp_%", "SpliceVarDB_%", "SpliceVault_%"]].to_markdown(index=False))
    md.append("")

    md.append("---")
    md.append("## 📐 Schema & Column Audit")
    if all_gene_columns:
        col_counts = [len(v) for v in all_gene_columns.values()]
        min_cols, max_cols = min(col_counts), max(col_counts)
        md.append(f"- **Column Count Range**: `{min_cols}` to `{max_cols}` columns per file.")
        if min_cols == max_cols:
            md.append(f"- **Schema Consistency**: ✅ **100% Consistent** across all {n_finished} completed files (`{min_cols}` columns total).")
        else:
            md.append(f"- **Schema Consistency**: ⚠️ **Inconsistency detected** across genes!")
            for g, diff in column_inconsistencies:
                md.append(f"  - `{g}` diff columns: {diff[:5]}...")
    md.append("")

    md.append("---")
    md.append("## ⚠️ Error & Risk Audit")
    if n_error == 0 and n_pending == 0:
        md.append("✅ **No critical errors or missing steps detected.** All input genes processed completely.")
    else:
        if n_pending > 0:
            md.append(f"⚠️ **{n_pending} genes are currently pending completion**.")
        if n_error > 0:
            md.append(f"❌ **{n_error} genes experienced issues during processing**.")
    md.append("")

    report_md_text = "\n".join(md)

    os.makedirs(os.path.dirname(output_report_path), exist_ok=True)
    with open(output_report_path, "w") as f:
        f.write(report_md_text)

    logger.info(f"Audit report saved to: {output_report_path}")
    print("\n" + "="*50)
    print(f"  Pipeline Audit Summary ({run_name})")
    print("="*50)
    print(f"Finished Genes: {n_finished} / {n_total}")
    print(f"Pending Genes : {n_pending}")
    print(f"Error Genes   : {n_error}")
    print(f"Input Vars    : {total_in_vars:,}")
    print(f"Output Vars   : {total_out_vars:,}")
    print(f"Report File   : {output_report_path}")
    print("="*50 + "\n")

    return df_results

if __name__ == "__main__":
    args = parse_args()
    audit_run(args.run_dir, args.raw_dir, args.output_report)
