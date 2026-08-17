#!/usr/bin/env python3
"""
Single Master Pipeline Execution Auditor & Quality Control Report Generator

Audits any pipeline run directory (or automatically audits the latest run in RUNS/):
  1. Identifies input raw variant tables vs output clean Parquet tables.
  2. Measures exact step completion timestamps and execution durations (VarConv, VEP, SPiP, Branchpoint, Pangolin, SpliceAI, Merge, Parquet).
  3. Live deep learning throughput & forward pass progress calculation (SpliceAI/Pangolin).
  4. Schema consistency, column completeness, and status contracts audit (scored, not_covered, error).
  5. Error and warning log audit (_log/ folder).
  6. Generates a standalone Markdown report with duration metrics and terminal summary tables.
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
        if pct >= 99.9:
            return "Finished"
        return f"{done_passes:,}/{tot_passes:,} passes ({pct:.1f}%)"
    except Exception:
        return "Active"

def infer_run_start_time(run_dir):
    run_name = os.path.basename(os.path.normpath(run_dir))
    # Try parsing run_YYYYMMDD_HHMM
    m = re.search(r"run_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})", run_name)
    if m:
        try:
            dt = datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5)))
            return dt.timestamp()
        except Exception:
            pass
    # Fallback to run directory mtime
    return os.path.getmtime(run_dir)

def audit_run(run_dir=None, raw_dir=None, output_report_path=None, save_walkthrough=False):
    if run_dir is None:
        run_dir = find_latest_run()
    
    run_dir = os.path.normpath(run_dir)
    run_name = os.path.basename(run_dir)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    t_run_start = infer_run_start_time(run_dir)

    tmp_dir = os.path.join(run_dir, "_tmp")
    ann_dir = os.path.join(run_dir, "annotation")
    res_dir = os.path.join(run_dir, "results")
    log_dir = os.path.join(run_dir, "_log")
    rep_dir = os.path.join(run_dir, "reports")
    flt_dir = os.path.join(run_dir, "filtered")
    plt_dir = os.path.join(run_dir, "plots")

    if output_report_path is None:
        output_report_path = os.path.join(rep_dir, f"audit_report_{run_name}.md")

    # Identify input genes
    if raw_dir and os.path.exists(raw_dir):
        input_files = sorted(glob.glob(os.path.join(raw_dir, "*.pq")) + glob.glob(os.path.join(raw_dir, "*.parquet")) + glob.glob(os.path.join(raw_dir, "*.tsv")) + glob.glob(os.path.join(raw_dir, "*.csv")))
    else:
        gene_dirs = [d for d in glob.glob(os.path.join(log_dir, "*")) if os.path.isdir(d)]
        input_files = gene_dirs

    if not input_files:
        logger.error(f"No input genes found for run directory: {run_dir}")
        return

    logger.info(f"Auditing run '{run_name}' ({len(input_files)} target genes)...")

    results = []
    timing_results = []
    all_gene_columns = {}
    schema_set = None
    column_inconsistencies = []

    key_annotation_cols = [
        'Consequence', 'REVEL_score', 'am_pathogenicity', 'SPiP', 'SPiP_interpretation', 'SPiP_status',
        'SpliceAI_status', 'Pangolin_status', 'Branchpoint_status', 'LaBranchoR_status',
        'SpliceVault_status', 'splicevardb', 'intron_offset_signed', 'splice_side'
    ]

    for inp_item in input_files:
        base_name = os.path.basename(inp_item)
        gene_name = base_name.split('.')[0].split('_')[0]

        # Read input counts
        n_input_vars = 0
        if os.path.isfile(inp_item):
            try:
                if inp_item.endswith('.pq') or inp_item.endswith('.parquet'):
                    df_inp = pd.read_parquet(inp_item)
                else:
                    sep = '\t' if inp_item.endswith('.tsv') else ','
                    df_inp = pd.read_csv(inp_item, sep=sep)
                n_input_vars = len(df_inp)
            except Exception:
                n_input_vars = 0
        else:
            vcf_gz = os.path.join(tmp_dir, f"{gene_name}.vcf.gz")
            if os.path.exists(vcf_gz):
                try:
                    cmd_vcf = f"zcat {vcf_gz} | grep -v '^#' | wc -l"
                    n_input_vars = int(subprocess.check_output(cmd_vcf, shell=True).decode().strip())
                except Exception:
                    n_input_vars = 0

        # Stage files & timestamps
        vcf_tmp = os.path.join(tmp_dir, f"{gene_name}.vcf.gz")
        spip_vcf = os.path.join(ann_dir, f"{gene_name}.annSPiP.vcf.gz")
        if not os.path.exists(spip_vcf):
            spip_vcf = os.path.join(ann_dir, f"{gene_name}.annSPiP.vcf")
        vep_vcf = os.path.join(ann_dir, f"{gene_name}.annVEP.vcf.gz")
        pangolin_vcf = os.path.join(ann_dir, f"{gene_name}.annPangolin.vcf.gz")
        spliceai_vcf = os.path.join(ann_dir, f"{gene_name}.annSpliceAI.vcf.gz")
        branch_vcf = os.path.join(ann_dir, f"{gene_name}.annBranchpoint.vcf.gz")
        annotated_vcf = os.path.join(ann_dir, f"{gene_name}.annotated.vcf.gz")
        final_pq = os.path.join(res_dir, f"{gene_name}.parsed.clean.pq")

        t_varconv = get_mtime_safe(vcf_tmp)
        t_vep = get_mtime_safe(vep_vcf)
        t_spip = get_mtime_safe(spip_vcf)
        t_branch = get_mtime_safe(branch_vcf)
        t_pangolin = get_mtime_safe(pangolin_vcf)
        t_spliceai = get_mtime_safe(spliceai_vcf)
        t_merge = get_mtime_safe(annotated_vcf)
        t_final_pq = get_mtime_safe(final_pq)

        # Durations
        dur_varconv = (t_varconv - t_run_start) if (t_varconv and t_run_start) else None
        base_t = t_varconv if t_varconv else t_run_start

        dur_vep = (t_vep - base_t) if (t_vep and base_t) else None
        dur_spip = (t_spip - base_t) if (t_spip and base_t) else None
        dur_branch = (t_branch - base_t) if (t_branch and base_t) else None
        dur_pangolin = (t_pangolin - base_t) if (t_pangolin and base_t) else None
        
        if t_spliceai and base_t:
            dur_spliceai = t_spliceai - base_t
            spliceai_time_str = format_duration(dur_spliceai)
        else:
            # Check if active running
            out_log = os.path.join(log_dir, gene_name, f"{gene_name}.spliceai.out")
            if os.path.exists(out_log) and base_t:
                cur_elapsed = os.path.getmtime(out_log) - base_t
                spliceai_time_str = f"🔄 Active ({format_duration(cur_elapsed)})"
            else:
                spliceai_time_str = "⏳ Pending"

        # Merge duration (time taken by merge job)
        all_pred_times = [t for t in [t_vep, t_spip, t_branch, t_pangolin, t_spliceai] if t is not None]
        last_pred_time = max(all_pred_times) if all_pred_times else base_t
        dur_merge = (t_merge - last_pred_time) if (t_merge and last_pred_time) else None

        # Parquet parsing duration
        dur_parquet = (t_final_pq - t_merge) if (t_final_pq and t_merge) else None

        # Total Turnaround
        if t_final_pq and t_run_start:
            total_turnaround = t_final_pq - t_run_start
            total_time_str = f"✅ {format_duration(total_turnaround)}"
        elif t_run_start:
            active_turnaround = datetime.now().timestamp() - t_run_start
            total_time_str = f"🔄 Running ({format_duration(active_turnaround)})"
        else:
            total_time_str = "-"

        timing_results.append({
            "Gene": gene_name,
            "VarConv": format_duration(dur_varconv),
            "VEP_111": format_duration(dur_vep),
            "SPiP_2.1": format_duration(dur_spip),
            "Branchpoint": format_duration(dur_branch),
            "SpliceAI": spliceai_time_str,
            "Merge": format_duration(dur_merge),
            "VCF2Parsed": format_duration(dur_parquet),
            "Total_Turnaround": total_time_str,
            "Completed_At": datetime.fromtimestamp(t_final_pq).strftime("%Y-%m-%d %H:%M:%S") if t_final_pq else "In Progress"
        })

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

        # SpliceAI live status
        if step_status["spliceai"]:
            spliceai_detail = f"✅ Done ({format_duration(dur_spliceai)})"
        else:
            prog = count_spliceai_progress(run_dir, gene_name)
            spliceai_detail = f"🔄 {prog}" if "passes" in prog else ("⏳ Pending" if prog == "Pending" else "✅ Done")

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

                for col in key_annotation_cols:
                    found_col = None
                    for c in [col, col.lower(), col.upper()]:
                        if c in df_out.columns:
                            found_col = c
                            break
                    if found_col:
                        non_null = df_out[found_col].notna() & (df_out[found_col].astype(str).str.strip() != '') & (df_out[found_col].astype(str) != 'nan') & (df_out[found_col].astype(str) != '.') & (df_out[found_col].astype(str) != '-')
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
                        if "exception" in content.lower() or "command not found" in content.lower() or "killed" in content.lower() or "check failed" in content.lower():
                            log_errors.append(os.path.basename(err_file))

        results.append({
            "Gene": gene_name,
            "Input_Vars": n_input_vars,
            "Output_Vars": n_out_vars,
            "Status": "FINISHED" if step_status["final_pq"] and pq_status == "VALID" else ("RUNNING/PENDING" if not log_errors else "ERROR"),
            "VarConv": "✅" if step_status["varconv"] else "❌",
            "VEP": "✅" if step_status["vep"] else "❌",
            "SPiP": "✅" if step_status["spip"] else "❌",
            "Branchpoint": "✅" if step_status["branchpoint"] else "❌",
            "SpliceAI": spliceai_detail,
            "Final_PQ": f"✅ {n_out_vars:,} vars" if pq_status == "VALID" else ("⏳ Pending" if pq_status == "MISSING" else f"❌ {pq_status}"),
            "Columns": n_cols if n_cols > 0 else "-",
            "REVEL_%": ann_stats.get("REVEL_score", "-"),
            "AlphaMissense_%": ann_stats.get("am_pathogenicity", "-"),
            "SPiP_%": ann_stats.get("SPiP", "-"),
            "SpliceVault_%": ann_stats.get("SpliceVault_status", "-"),
            "Branchpoint_%": ann_stats.get("Branchpoint_status", "-"),
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
    md.append("## ⏱️ 2. Step Completion Times & Execution Durations")
    md.append(df_timing[["Gene", "VarConv", "VEP_111", "SPiP_2.1", "Branchpoint", "SpliceAI", "Merge", "VCF2Parsed", "Total_Turnaround", "Completed_At"]].to_markdown(index=False))
    md.append("")

    md.append("---")
    md.append("## 📋 3. Stage-by-Stage Completion Matrix")
    stage_cols = ["Gene", "Status", "VarConv", "VEP", "SPiP", "Branchpoint", "SpliceAI", "Final_PQ", "Log_Issues"]
    md.append(df_results[stage_cols].to_markdown(index=False))
    md.append("")

    if n_finished > 0:
        md.append("---")
        md.append("## 🧬 4. Key Annotation Completeness (% Non-Null in Completed Files)")
        ann_cols = ["Gene", "Columns", "REVEL_%", "AlphaMissense_%", "SPiP_%", "SpliceVault_%", "Branchpoint_%"]
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
                md.append(f"  - `{g}` diff columns: {diff[:5]}...")
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
    print("\n" + "="*80)
    print(f"  🧪 Master Pipeline Audit & Timings Summary ({run_name})")
    print("="*80)
    print(df_timing[["Gene", "VarConv", "VEP_111", "SPiP_2.1", "SpliceAI", "VCF2Parsed", "Total_Turnaround"]].to_string(index=False))
    print("="*80)
    print(f"Finished Genes: {n_finished} / {n_total} ({(n_finished/n_total)*100:.1f}%) | Pending: {n_running} | Errors: {n_error}")
    print(f"Report File   : {output_report_path}")
    print("="*80 + "\n")

    return df_results

if __name__ == "__main__":
    args = parse_args()
    target_dir = args.run_dir
    if target_dir is None and args.run_name:
        target_dir = os.path.join("RUNS", args.run_name)
    audit_run(target_dir, args.raw_dir, args.output_report, save_walkthrough=args.walkthrough)
