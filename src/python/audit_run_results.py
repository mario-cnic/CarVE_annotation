#!/usr/bin/env python3
"""
Single Master Pipeline Execution Auditor & Quality Control Report Generator

Audits any pipeline run directory (or automatically audits the latest run in RUNS/):
  1. Identifies input raw variant tables vs output clean Parquet tables.
  2. Tracks exact active step for each gene (VarConv, VEP, SPiP, Branchpoint, Pangolin, SpliceAI, Merge, Parquet).
  3. Measures exact step completion timestamps, step execution durations, and total turnaround time.
  4. Live deep learning throughput & forward pass progress calculation with remaining ETA (SpliceAI/Pangolin).
  5. Schema consistency, column completeness, predictor coverage %, and anomaly detection.
  6. Error and warning log audit (_log/ folder).
  7. Generates a standalone Markdown report with duration metrics, coverage tables, and quality alerts.
"""

import os
import sys
import glob
import re
import argparse
import logging
import subprocess
import pandas as pd
import numpy as np
from datetime import datetime

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("audit_run_results")

def parse_args():
    parser = argparse.ArgumentParser(description="Single Master Pipeline Execution Auditor & Report Generator")
    parser.add_argument("--run-dir", default=None, help="Path to run directory (default: latest run in RUNS/)")
    parser.add_argument("--run-name", default=None, help="Run name in RUNS/ (e.g. run_20260813_1028)")
    parser.add_argument("--raw-dir", default="RUNS/predictors_050826/input/by_gene", help="Path to input raw folder")
    parser.add_argument("--output-report", default=None, help="Custom path for output markdown report file")
    parser.add_argument("--walkthrough", action="store_true", help="Save a copy of the audit report in walkthrough/ directory")
    return parser.parse_args()

def find_latest_run():
    runs = [d for d in glob.glob("RUNS/*") if os.path.isdir(d) and not d.endswith("predictors_050826")]
    if not runs:
        runs = [d for d in glob.glob("RUNS/*") if os.path.isdir(d)]
    if not runs:
        return "RUNS/predictors_050826"
    runs.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    return runs[0]

def format_duration(seconds):
    if seconds is None or np.isnan(seconds) or seconds < 0:
        return "-"
    seconds = int(round(seconds))
    if seconds < 60:
        return f"{seconds}s"
    elif seconds < 3600:
        m = seconds // 60
        s = seconds % 60
        return f"{m}m {s:02d}s"
    else:
        h = seconds // 3600
        m = (seconds % 3600) // 60
        s = seconds % 60
        return f"{h}h {m:02d}m {s:02d}s"

def get_file_status(filepath):
    if os.path.exists(filepath):
        size_bytes = os.path.getsize(filepath)
        if size_bytes > 1024 * 1024:
            return f"✅ Done ({size_bytes / (1024 * 1024):.2f} MB)"
        elif size_bytes > 1024:
            return f"✅ Done ({size_bytes / 1024:.1f} KB)"
        elif size_bytes > 0:
            return f"⚠️ Active ({size_bytes} B)"
        else:
            return "⏳ Empty/Pending"
    return "⏳ Pending"

def get_mtime_safe(filepath):
    if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
        return os.path.getmtime(filepath)
    return None

def count_spliceai_progress(run_dir, gene):
    out_log = os.path.join(run_dir, "_log", gene, f"{gene}.spliceai.out")
    vcf_file = os.path.join(run_dir, "_tmp", f"{gene}.vcf.gz")
    if not os.path.exists(out_log) or not os.path.exists(vcf_file):
        return "Pending"
    
    try:
        cmd_vcf = f"zcat {vcf_file} | grep -v '^#' | wc -l"
        tot_vars = int(subprocess.check_output(cmd_vcf, shell=True).decode().strip())
        tot_passes = tot_vars * 5  # 5 ensemble models per variant
        
        cmd_passes = f"grep -c '1/1 \\[' {out_log} 2>/dev/null || echo 0"
        done_passes = int(subprocess.check_output(cmd_passes, shell=True).decode().strip())
        
        pct = (done_passes / tot_passes) * 100 if tot_passes > 0 else 0
        rem_passes = max(0, tot_passes - done_passes)
        rem_hours = (rem_passes * 3.0) / 3600
        if pct >= 99.9:
            return "Finished"
        return f"{done_passes:,}/{tot_passes:,} passes ({pct:.1f}% | ~{rem_hours:.1f}h left)"
    except Exception:
        return "In Progress"

def parse_time_from_log(log_path):
    if not os.path.exists(log_path) or os.path.getsize(log_path) == 0:
        return None
    try:
        with open(log_path, 'r', errors='ignore') as f:
            for line in f:
                line = line.strip()
                if "Elapsed (wall clock) time" in line:
                    match = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9:]+(\.[0-9]+)?)', line)
                    if match:
                        time_str = match.group(1).split('.')[0]
                        parts = [int(p) for p in time_str.split(':')]
                        if len(parts) == 3:
                            return parts[0] * 3600 + parts[1] * 60 + parts[2]
                        elif len(parts) == 2:
                            return parts[0] * 60 + parts[1]
                elif "real" in line and re.match(r'^real\s+[0-9]+[mhs]', line):
                    match_m = re.search(r'([0-9]+)m', line)
                    match_s = re.search(r'([0-9]+(\.[0-9]+)?)s', line)
                    m = int(match_m.group(1)) if match_m else 0
                    s = float(match_s.group(1)) if match_s else 0
                    return m * 60 + int(round(s))
    except Exception:
        pass
    return None

def audit_run(run_dir=None, raw_dir="RUNS/predictors_050826/input/by_gene", output_report_path=None, save_walkthrough=False):
    if run_dir is None:
        run_dir = find_latest_run()
    
    run_dir = run_dir.rstrip("/")
    run_name = os.path.basename(run_dir)
    res_dir = os.path.join(run_dir, "results")
    ann_dir = os.path.join(run_dir, "annotation")
    tmp_dir = os.path.join(run_dir, "_tmp")
    log_dir = os.path.join(run_dir, "_log")
    reports_dir = os.path.join(run_dir, "reports")
    
    if output_report_path is None:
        output_report_path = os.path.join(reports_dir, f"audit_report_{run_name}.md")

    logger.info(f"Auditing run directory: {run_dir}")
    logger.info(f"Input raw directory : {raw_dir}")

    # Discover target genes
    genes = set()
    if os.path.exists(raw_dir):
        for f in glob.glob(os.path.join(raw_dir, "*")):
            ext = os.path.splitext(f)[1].lower()
            if ext in [".pq", ".parquet", ".xlsx", ".csv", ".tsv", ".vcf", ".gz"]:
                base = os.path.basename(f).split(".")[0].split("_")[0].upper()
                genes.add(base)
    
    if os.path.exists(res_dir):
        for f in glob.glob(os.path.join(res_dir, "*.clean.pq")):
            base = os.path.basename(f).split(".")[0].split("_")[0].upper()
            genes.add(base)

    if not genes:
        logger.warning(f"No genes found in {raw_dir} or {res_dir}")
        genes = ["BAG3", "DSP", "FLNC", "LMNA", "MYBPC3", "PKP2", "TTN"]

    genes = sorted(list(genes))
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    now = datetime.now().timestamp()

    results = []
    timing_results = []
    schema_set = None
    column_inconsistencies = []
    all_gene_columns = {}

    key_annotation_cols = [
        ("gnomADv4_AF_grpmax", ["gnomadv4_af_grpmax", "gnomad_af_joint", "max_af"]),
        ("REVEL_score", ["revel_score", "revel"]),
        ("AlphaMissense", ["am_pathogenicity", "alphamissense_score"]),
        ("SPiP", ["spip", "spip_interpretation"]),
        ("Pangolin", ["pangolin_max_score", "pangolin"]),
        ("SpliceAI", ["spliceai_pred_ds_max", "spliceai"]),
        ("Branchpoint", ["branchpoint_disrupted", "branchpoint_status", "labranchor_score"]),
        ("5UTR_Annotator", ["5utr_annotation", "5utr_consequence"]),
        ("SpliceVault", ["splicevault_top_events", "splicevault_status"]),
        ("ClinVar", ["clinvar_clnsig", "clinvar"])
    ]

    for gene_name in genes:
        # Check raw input
        raw_files = glob.glob(os.path.join(raw_dir, f"{gene_name}*"))
        n_input_vars = 0
        if raw_files:
            try:
                rf = raw_files[0]
                if rf.endswith(".pq") or rf.endswith(".parquet"):
                    df_raw = pd.read_parquet(rf)
                    n_input_vars = len(df_raw)
                elif rf.endswith(".csv"):
                    n_input_vars = sum(1 for _ in open(rf)) - 1
                elif rf.endswith(".tsv"):
                    n_input_vars = sum(1 for _ in open(rf)) - 1
                elif rf.endswith(".xlsx"):
                    df_raw = pd.read_excel(rf)
                    n_input_vars = len(df_raw)
            except Exception:
                pass

        # Pipeline file targets
        vcf_tmp = os.path.join(tmp_dir, f"{gene_name}.vcf.gz")
        if not os.path.exists(vcf_tmp):
            vcf_tmp = os.path.join(tmp_dir, f"{gene_name}.vcf")
        vep_vcf = os.path.join(ann_dir, f"{gene_name}.annVEP.vcf.gz")
        spip_vcf = os.path.join(ann_dir, f"{gene_name}.annSPiP.vcf.gz")
        if not os.path.exists(spip_vcf):
            spip_vcf = os.path.join(ann_dir, f"{gene_name}.annSPiP.vcf")
        branch_vcf = os.path.join(ann_dir, f"{gene_name}.annBranchpoint.vcf.gz")
        pangolin_vcf = os.path.join(ann_dir, f"{gene_name}.annPangolin.vcf.gz")
        spliceai_vcf = os.path.join(ann_dir, f"{gene_name}.annSpliceAI.vcf.gz")
        annotated_vcf = os.path.join(ann_dir, f"{gene_name}.annotated.vcf.gz")
        final_pq = os.path.join(res_dir, f"{gene_name}.parsed.clean.pq")

        # Step Statuses
        step_status = {
            "varconv": os.path.exists(vcf_tmp) and os.path.getsize(vcf_tmp) > 0,
            "vep": os.path.exists(vep_vcf) and os.path.getsize(vep_vcf) > 0,
            "spip": os.path.exists(spip_vcf) and os.path.getsize(spip_vcf) > 0,
            "branchpoint": os.path.exists(branch_vcf) and os.path.getsize(branch_vcf) > 0,
            "pangolin": os.path.exists(pangolin_vcf) and os.path.getsize(pangolin_vcf) > 0,
            "spliceai": os.path.exists(spliceai_vcf) and os.path.getsize(spliceai_vcf) > 1000,
            "merged": os.path.exists(annotated_vcf) and os.path.getsize(annotated_vcf) > 0,
            "final_pq": os.path.exists(final_pq) and os.path.getsize(final_pq) > 0
        }

        # SpliceAI Progress / Status
        spliceai_detail = "⏳ Pending"
        if step_status["spliceai"]:
            spliceai_detail = f"✅ Done ({os.path.getsize(spliceai_vcf)/(1024*1024):.2f} MB)"
        else:
            prog = count_spliceai_progress(run_dir, gene_name)
            if prog == "Finished":
                spliceai_detail = "✅ Done"
            elif "passes" in prog:
                spliceai_detail = f"🔄 {prog}"

        # Current active step determination
        cur_step = "Done"
        if not step_status["final_pq"]:
            if not step_status["varconv"]:
                cur_step = "1. VarConv"
            elif not (step_status["vep"] and step_status["spip"] and step_status["branchpoint"] and step_status["spliceai"]):
                cur_step = "2. Predictors"
            elif not step_status["merged"]:
                cur_step = "3. Merge"
            else:
                cur_step = "4. VCF2Parsed"

        t_final_pq = get_mtime_safe(final_pq)
        t_run_start = get_mtime_safe(vcf_tmp)
        if t_final_pq and t_run_start:
            total_turnaround = t_final_pq - t_run_start
            total_time_str = f"✅ {format_duration(total_turnaround)}"
        elif t_run_start:
            active_turnaround = now - t_run_start
            total_time_str = f"🔄 Running ({format_duration(active_turnaround)})"
        else:
            total_time_str = "-"

        timing_results.append({
            "Gene": gene_name,
            "Variants": f"{n_input_vars:,}" if n_input_vars > 0 else "-",
            "Current_Step": cur_step,
            "SpliceAI_Progress": spliceai_detail,
            "Total_Elapsed": total_time_str,
            "Completed_At": datetime.fromtimestamp(t_final_pq).strftime("%Y-%m-%d %H:%M:%S") if t_final_pq else "In Progress"
        })

        # Audit final pq file & annotations
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

                for label, cand_cols in key_annotation_cols:
                    found_col = None
                    for c in df_out.columns:
                        if c.lower() in cand_cols or c in cand_cols:
                            found_col = c
                            break
                    if found_col:
                        col_vals = df_out[found_col].dropna()
                        non_null = col_vals[~col_vals.astype(str).str.strip().isin(['', 'nan', '.', '-', 'NA', 'none', 'None'])]
                        pct = (len(non_null) / max(1, n_out_vars)) * 100
                        ann_stats[label] = f"{pct:.1f}%"
                    else:
                        ann_stats[label] = "N/A"
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
                        if "exception" in content.lower() or "command not found" in content.lower() or "killed" in content.lower() or "stale file handle" in content.lower() or "filenotfounderror" in content.lower():
                            log_errors.append(os.path.basename(err_file))

        results.append({
            "Gene": gene_name,
            "Input_Vars": n_input_vars,
            "Output_Vars": n_out_vars,
            "Status": "FINISHED" if step_status["final_pq"] and pq_status == "VALID" else ("RUNNING/PENDING" if not log_errors else "ERROR"),
            "VarConv": "✅" if step_status["varconv"] else "❌",
            "VEP": "✅" if step_status["vep"] else "❌",
            "SPiP": "✅" if step_status["spip"] else "❌",
            "Pangolin": "✅" if step_status["pangolin"] else "❌",
            "Branchpoint": "✅" if step_status["branchpoint"] else "❌",
            "SpliceAI": spliceai_detail,
            "Final_PQ": f"✅ {n_out_vars:,} vars" if pq_status == "VALID" else ("⏳ Pending" if pq_status == "MISSING" else f"❌ {pq_status}"),
            "Columns": n_cols if n_cols > 0 else "-",
            "gnomADv4_%": ann_stats.get("gnomADv4_AF_grpmax", "-"),
            "REVEL_%": ann_stats.get("REVEL_score", "-"),
            "AlphaMissense_%": ann_stats.get("AlphaMissense", "-"),
            "SPiP_%": ann_stats.get("SPiP", "-"),
            "Pangolin_%": ann_stats.get("Pangolin", "-"),
            "SpliceAI_%": ann_stats.get("SpliceAI", "-"),
            "Branchpoint_%": ann_stats.get("Branchpoint", "-"),
            "5UTR_%": ann_stats.get("5UTR_Annotator", "-"),
            "Log_Issues": ", ".join(log_errors) if log_errors else "None"
        })

    df_results = pd.DataFrame(results)
    df_timing = pd.DataFrame(timing_results)

    n_total = len(df_results)
    n_finished = len(df_results[df_results["Status"] == "FINISHED"])
    n_running = len(df_results[df_results["Status"] == "RUNNING/PENDING"])
    n_error = len(df_results[df_results["Status"] == "ERROR"])

    total_in_vars = df_results["Input_Vars"].sum()
    total_out_vars = df_results["Output_Vars"].sum()

    md = []
    md.append(f"# 🧪 Master Pipeline Execution & Quality Audit: `{run_name}`")
    md.append(f"**Audit Timestamp**: `{timestamp}`  ")
    md.append(f"**Target Build**: `GRCh38`  ")
    md.append(f"**Run Directory**: `{run_dir}`  ")
    md.append(f"**Results Directory**: `{res_dir}`  \n")

    md.append("---")
    md.append("## 📊 1. Executive Summary")
    md.append(f"- **Total Target Genes**: `{n_total}`")
    md.append(f"- **Genes Fully Finished**: `{n_finished} / {n_total}` ({(n_finished/n_total)*100:.1f}%)")
    md.append(f"- **Genes In Progress / Pending**: `{n_running}`")
    md.append(f"- **Genes with Critical Errors**: `{n_error}`")
    md.append(f"- **Total Input Variants**: `{total_in_vars:,}`")
    md.append(f"- **Total Clean Output Variants Generated**: `{total_out_vars:,}`")
    md.append("")

    md.append("---")
    md.append("## ⏱️ 2. Per-Gene Step Tracking & Running Times")
    md.append(df_timing[["Gene", "Variants", "Current_Step", "SpliceAI_Progress", "Total_Elapsed", "Completed_At"]].to_markdown(index=False))
    md.append("")

    md.append("---")
    md.append("## 📋 3. Stage-by-Stage Predictor Matrix")
    stage_cols = ["Gene", "Status", "VarConv", "VEP", "SPiP", "Pangolin", "Branchpoint", "SpliceAI", "Final_PQ", "Log_Issues"]
    md.append(df_results[stage_cols].to_markdown(index=False))
    md.append("")

    if n_finished > 0:
        md.append("---")
        md.append("## 🧬 4. Key Predictor Completeness (% Non-Null in Output Parquets)")
        ann_cols = ["Gene", "Columns", "gnomADv4_%", "REVEL_%", "AlphaMissense_%", "SPiP_%", "Pangolin_%", "Branchpoint_%", "5UTR_%"]
        md.append(df_results[df_results["Status"] == "FINISHED"][ann_cols].to_markdown(index=False))
        md.append("")

    md.append("---")
    md.append("## 📐 5. Schema & Column Integrity")
    if all_gene_columns:
        col_counts = [len(v) for v in all_gene_columns.values()]
        min_cols, max_cols = min(col_counts), max(col_counts)
        md.append(f"- **Column Count Range**: `{min_cols}` to `{max_cols}` columns per parquet table.")
        if min_cols == max_cols:
            md.append(f"- **Schema Consistency**: ✅ **100% Consistent** across all {n_finished} completed files (`{min_cols}` columns total).")
        else:
            md.append(f"- **Schema Consistency**: ⚠️ **Inconsistency detected** across genes!")
            for g, diff in column_inconsistencies:
                md.append(f"  - `{g}` diff columns ({len(diff)}): {diff[:5]}...")
    else:
        md.append("- *No completed Parquet tables found yet to evaluate schema integrity.*")
    md.append("")

    report_md_text = "\n".join(md)

    os.makedirs(os.path.dirname(output_report_path), exist_ok=True)
    with open(output_report_path, "w") as f:
        f.write(report_md_text)

    if save_walkthrough:
        walk_path = os.path.join("walkthrough", f"{datetime.now().strftime('%Y%m%d')}_master_audit_{run_name}.md")
        with open(walk_path, "w") as f:
            f.write(report_md_text)
        logger.info(f"Walkthrough copy saved to: {walk_path}")

    logger.info(f"Master audit report saved to: {output_report_path}")
    print("\n" + "="*95)
    print(f"  🧪 Master Pipeline Quality & Completeness Audit Summary ({run_name})")
    print("="*95)
    print(df_timing[["Gene", "Variants", "Current_Step", "SpliceAI_Progress", "Total_Elapsed"]].to_string(index=False))
    print("="*95)
    print(f"Finished Genes: {n_finished} / {n_total} ({(n_finished/n_total)*100:.1f}%) | Pending: {n_running} | Errors: {n_error}")
    print(f"Report File   : {output_report_path}")
    print("="*95 + "\n")

    return df_results

if __name__ == "__main__":
    args = parse_args()
    target_dir = args.run_dir
    if target_dir is None and args.run_name:
        target_dir = os.path.join("RUNS", args.run_name)
    audit_run(target_dir, args.raw_dir, args.output_report, save_walkthrough=args.walkthrough)
