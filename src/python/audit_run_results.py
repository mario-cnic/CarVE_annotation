#!/usr/bin/env python3
"""
Single Master Pipeline Execution & Quality Control Auditor

Audits genomic variant annotation pipeline runs:
  1. Identifies input raw variant tables vs output clean Parquet tables.
  2. Tracks exact active step per gene (VarConv, Predictors, Merge, VCF2Parsed, Done).
  3. Live deep learning throughput & forward pass progress calculation with remaining ETA (SpliceAI/Pangolin).
  4. Predictor status matrix (VarConv, VEP, SPiP, Pangolin, Branchpoint, SpliceAI, Final PQ).
  5. Schema consistency, column completeness, and key predictor coverage %.
  6. Error and warning log audit (_log/ folder).
  7. Saves the master markdown report ONLY inside the target run's reports/ folder:
     RUNS/<RUN_NAME>/reports/audit_report_<RUN_NAME>.md
"""

import os
import sys

# Prevent CWD sys.path precedence from shadowing site-packages C extension modules (numpy/pandas)
sys.path = [p for p in sys.path if p not in ("", ".")]

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
    parser.add_argument("--run-dir", default=None, help="Path to run directory (e.g. RUNS/run_20260813_1028)")
    parser.add_argument("--run-name", default=None, help="Run name inside RUNS/ (e.g. run_20260813_1028)")
    parser.add_argument("--raw-dir", default=None, help="Path to input raw files folder (optional)")
    parser.add_argument("--output-report", default=None, help="Custom path for output markdown report file")
    parser.add_argument("--all-runs", "--all", action="store_true", help="Audit all run folders found in RUNS/")
    return parser.parse_args()

def discover_available_runs():
    runs = [d for d in glob.glob("RUNS/*") if os.path.isdir(d) and not d.endswith("predictors_050826")]
    if not runs:
        runs = [d for d in glob.glob("RUNS/*") if os.path.isdir(d)]
    runs.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    return runs

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

def get_mtime_safe(filepath):
    if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
        return os.path.getmtime(filepath)
    return None

def count_vcf_variants(vcf_path):
    if not os.path.exists(vcf_path) or os.path.getsize(vcf_path) == 0:
        return 0
    try:
        if vcf_path.endswith(".gz"):
            if os.path.exists(vcf_path + ".tbi") or os.path.exists(vcf_path + ".csi"):
                for b_bin in ["/data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools", "/home/mruizp/apps/miniforge3/envs/genomics/bin/bcftools", "bcftools"]:
                    try:
                        res = subprocess.check_output(f"{b_bin} index -n '{vcf_path}' 2>/dev/null", shell=True).decode().strip()
                        if res.isdigit() and int(res) > 0:
                            return int(res)
                    except Exception:
                        pass
            cmd = f"zgrep -c '^[^#]' '{vcf_path}' 2>/dev/null || echo 0"
        else:
            cmd = f"grep -c '^[^#]' '{vcf_path}' 2>/dev/null || echo 0"
        res = subprocess.check_output(cmd, shell=True).decode().strip()
        return int(res) if res.isdigit() else 0
    except Exception:
        return 0

def check_vcf_completeness(vcf_path, expected_vars):
    if not os.path.exists(vcf_path) or os.path.getsize(vcf_path) == 0:
        return {"status": "PENDING", "count": 0, "pct": 0.0, "is_complete": False}
    
    n_vars = count_vcf_variants(vcf_path)
    pct = (n_vars / expected_vars) * 100 if expected_vars > 0 else 0
    is_complete = (n_vars >= expected_vars) if expected_vars > 0 else (n_vars > 0)
    
    return {
        "status": "DONE" if is_complete else "INCOMPLETE",
        "count": n_vars,
        "pct": pct,
        "is_complete": is_complete,
        "size_mb": os.path.getsize(vcf_path) / (1024 * 1024)
    }

def check_spliceai_status(run_dir, gene, n_input_vars, vcf_tmp_mtime=None):
    ann_dir = os.path.join(run_dir, "annotation")
    tmp_dir = os.path.join(run_dir, "_tmp")
    log_dir = os.path.join(run_dir, "_log", gene)
    
    final_vcf_gz = os.path.join(ann_dir, f"{gene}.annSpliceAI.vcf.gz")
    final_vcf = os.path.join(ann_dir, f"{gene}.annSpliceAI.vcf")
    out_log = os.path.join(log_dir, f"{gene}.spliceai.out")
    vcf_tmp = os.path.join(tmp_dir, f"{gene}.vcf.gz")
    if not os.path.exists(vcf_tmp):
        vcf_tmp = os.path.join(tmp_dir, f"{gene}.vcf")

    log_mtime = get_mtime_safe(out_log)
    target_vcf = final_vcf_gz if os.path.exists(final_vcf_gz) else (final_vcf if os.path.exists(final_vcf) else None)
    vcf_mtime = get_mtime_safe(target_vcf) if target_vcf else None

    # Check parallel chunking directory
    chunk_dirs = glob.glob(os.path.join(tmp_dir, f"{gene}_spliceai_chunks*")) + \
                 glob.glob(os.path.join(tmp_dir, f"spliceai_chunks_{gene}*"))
    if chunk_dirs:
        cdir = chunk_dirs[0]
        chunks_all = glob.glob(os.path.join(cdir, "chunk_*.vcf.gz"))
        chunks_done = glob.glob(os.path.join(cdir, "chunk_*.annSpliceAI.vcf.gz"))
        if chunks_all and len(chunks_done) < len(chunks_all):
            pct = (len(chunks_done) / len(chunks_all)) * 100
            return {"status": "IN_PROGRESS", "detail": f"🔄 Chunking: {len(chunks_done)}/{len(chunks_all)} chunks ({pct:.1f}%)", "is_complete": False}

    # Check forward-pass log activity
    passes_detail = None
    done_passes = 0
    tot_passes = n_input_vars * 5
    if os.path.exists(out_log) and n_input_vars > 0:
        try:
            cmd_passes = f"grep -c '1/1 \\[' '{out_log}' 2>/dev/null || echo 0"
            done_passes = int(subprocess.check_output(cmd_passes, shell=True).decode().strip())
            pct_passes = (done_passes / tot_passes) * 100 if tot_passes > 0 else 0
            rem_passes = max(0, tot_passes - done_passes)
            rem_hours = (rem_passes * 0.35) / 3600
            if done_passes < tot_passes:
                passes_detail = f"🔄 {done_passes:,}/{tot_passes:,} passes ({pct_passes:.1f}% | ~{rem_hours:.1f}h left)"
        except Exception:
            pass

    # Check VCF record completeness
    if target_vcf:
        vcf_chk = check_vcf_completeness(target_vcf, n_input_vars)
        if vcf_chk["is_complete"]:
            return {
                "status": "DONE",
                "detail": f"✅ Done ({vcf_chk['size_mb']:.2f} MB)",
                "is_complete": True,
                "count": vcf_chk["count"]
            }
        else:
            # File is partial / currently being written
            if passes_detail:
                detail = f"{passes_detail} [VCF: {vcf_chk['count']:,}/{n_input_vars:,} ({vcf_chk['pct']:.1f}%)]"
            else:
                detail = f"🔄 Writing: {vcf_chk['count']:,}/{n_input_vars:,} vars ({vcf_chk['pct']:.1f}%)"
            return {
                "status": "IN_PROGRESS",
                "detail": detail,
                "is_complete": False,
                "count": vcf_chk["count"]
            }

    if passes_detail:
        return {"status": "IN_PROGRESS", "detail": passes_detail, "is_complete": False}

    return {"status": "PENDING", "detail": "⏳ Pending", "is_complete": False}

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

def audit_run(run_dir=None, raw_dir=None, output_report_path=None):
    if run_dir is None:
        runs = discover_available_runs()
        run_dir = runs[0] if runs else "RUNS/run_20260813_1028"
    
    run_dir = run_dir.rstrip("/")
    run_name = os.path.basename(run_dir)
    res_dir = os.path.join(run_dir, "results")
    ann_dir = os.path.join(run_dir, "annotation")
    tmp_dir = os.path.join(run_dir, "_tmp")
    log_dir = os.path.join(run_dir, "_log")
    reports_dir = os.path.join(run_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    if not output_report_path:
        output_report_path = os.path.join(reports_dir, f"audit_report_{run_name}.md")

    logger.info(f"Auditing run directory: {run_dir}")
    if raw_dir:
        logger.info(f"Input raw directory : {raw_dir}")

    # Discover target genes comprehensively
    genes = set()
    if raw_dir and os.path.exists(raw_dir):
        for f in glob.glob(os.path.join(raw_dir, "*")):
            ext = os.path.splitext(f)[1].lower()
            if ext in [".pq", ".parquet", ".xlsx", ".csv", ".tsv", ".vcf", ".gz"]:
                base = os.path.basename(f).split(".")[0].split("_")[0]
                if base and not base.startswith("chunk"):
                    genes.add(base)
    
    if os.path.exists(tmp_dir):
        for f in glob.glob(os.path.join(tmp_dir, "*.vcf*")):
            base = os.path.basename(f).split(".")[0].split("_")[0]
            if base and not base.startswith("chunk"):
                genes.add(base)

    if os.path.exists(res_dir):
        for f in glob.glob(os.path.join(res_dir, "*.clean.pq")):
            base = os.path.basename(f).split(".")[0].split("_")[0]
            if base:
                genes.add(base)

    if os.path.exists(ann_dir):
        for f in glob.glob(os.path.join(ann_dir, "*.annVEP.vcf.gz")):
            base = os.path.basename(f).split(".")[0].split("_")[0]
            if base:
                genes.add(base)

    if os.path.exists(log_dir):
        for d in glob.glob(os.path.join(log_dir, "*")):
            if os.path.isdir(d):
                genes.add(os.path.basename(d))

    if not genes:
        logger.warning(f"No genes found in {run_dir}")
        genes = ["BAG3", "DSP", "FLNC", "LMNA", "MYBPC3", "PKP2", "TTN"]

    genes = sorted(list(genes))
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    now = datetime.now().timestamp()

    results = []
    timing_results = []
    version_lineage = []
    schema_set = None
    column_inconsistencies = []
    all_gene_columns = {}

    key_annotation_cols = [
        ("gnomADv4_AF_grpmax", ["gnomadv4_af_grpmax_joint", "gnomadv4_af_grpmax", "gnomad_af_joint", "max_af"]),
        ("REVEL_score", ["revel_score", "revel"]),
        ("AlphaMissense", ["am_pathogenicity", "alphamissense_score", "alphamissense"]),
        ("SPiP", ["spip_prediction", "spip", "spip_interpretation"]),
        ("Pangolin", ["pangolin_max_score", "pangolin_largest_delta", "pangolin"]),
        ("SpliceAI", ["spliceai_custom_max", "spliceai_max", "spliceai_pred_ds_max", "spliceai"]),
        ("Branchpoint", ["branchpoint_disrupted", "branchpoint_status", "labranchor_score"]),
        ("5UTR_Annotator", ["5utr_consequence", "5utr_annotation"]),
        ("SpliceVault", ["splicevault_status", "splicevault_top_events"]),
        ("ClinVar", ["clinvar_clnsig", "clinvar"])
    ]

    for gene_name in genes:
        # Check raw input
        n_input_vars = 0
        if raw_dir and os.path.exists(raw_dir):
            raw_files = glob.glob(os.path.join(raw_dir, f"{gene_name}*"))
            if raw_files:
                try:
                    rf = raw_files[0]
                    if rf.endswith(".pq") or rf.endswith(".parquet"):
                        df_raw = pd.read_parquet(rf)
                        n_input_vars = len(df_raw)
                    elif rf.endswith(".csv") or rf.endswith(".tsv"):
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

        if n_input_vars == 0 and os.path.exists(vcf_tmp):
            n_input_vars = count_vcf_variants(vcf_tmp)

        vep_vcf = os.path.join(ann_dir, f"{gene_name}.annVEP.vcf.gz")
        spip_vcf = os.path.join(ann_dir, f"{gene_name}.annSPiP.vcf.gz")
        if not os.path.exists(spip_vcf):
            spip_vcf = os.path.join(ann_dir, f"{gene_name}.annSPiP.vcf")
        branch_vcf = os.path.join(ann_dir, f"{gene_name}.annBranchpoint.vcf.gz")
        pangolin_vcf = os.path.join(ann_dir, f"{gene_name}.annPangolin.vcf.gz")
        spliceai_vcf = os.path.join(ann_dir, f"{gene_name}.annSpliceAI.vcf.gz")
        if not os.path.exists(spliceai_vcf):
            spliceai_vcf = os.path.join(ann_dir, f"{gene_name}.annSpliceAI.vcf")
        annotated_vcf = os.path.join(ann_dir, f"{gene_name}.annotated.vcf.gz")
        final_pq = os.path.join(res_dir, f"{gene_name}.parsed.clean.pq")

        # Step Statuses
        t_varconv = get_mtime_safe(vcf_tmp)
        t_vep = get_mtime_safe(vep_vcf)
        t_spip = get_mtime_safe(spip_vcf)
        t_branch = get_mtime_safe(branch_vcf)
        t_pangolin = get_mtime_safe(pangolin_vcf)
        t_spliceai = get_mtime_safe(spliceai_vcf)
        t_merged = get_mtime_safe(annotated_vcf)
        t_final_pq = get_mtime_safe(final_pq)

        # Check completeness for each predictor
        vep_chk = check_vcf_completeness(vep_vcf, n_input_vars)
        spip_chk = check_vcf_completeness(spip_vcf, n_input_vars)
        branch_chk = check_vcf_completeness(branch_vcf, n_input_vars)
        pangolin_chk = check_vcf_completeness(pangolin_vcf, n_input_vars)
        spliceai_info = check_spliceai_status(run_dir, gene_name, n_input_vars, t_varconv)
        spliceai_detail = spliceai_info["detail"]
        spliceai_is_complete = spliceai_info.get("is_complete", False)

        # Step Status Booleans (Strict Record Count Validation)
        step_status = {
            "varconv": os.path.exists(vcf_tmp) and os.path.getsize(vcf_tmp) > 0,
            "vep": vep_chk["count"] >= n_input_vars if n_input_vars > 0 else vep_chk["is_complete"],
            "spip": spip_chk["is_complete"],
            "branchpoint": branch_chk["is_complete"],
            "pangolin": pangolin_chk["is_complete"],
            "spliceai": spliceai_is_complete,
            "merged": os.path.exists(annotated_vcf) and os.path.getsize(annotated_vcf) > 1000,
            "final_pq": os.path.exists(final_pq) and os.path.getsize(final_pq) > 0
        }

        # Predictor status strings for the matrix
        def format_pred_status(chk, exists):
            if not exists:
                return "❌ Not created"
            if chk["is_complete"]:
                return "✅ Done"
            return f"⚠️ {chk['count']:,}/{n_input_vars:,} ({chk['pct']:.1f}%)"

        pred_vep_str = format_pred_status(vep_chk, os.path.exists(vep_vcf))
        pred_spip_str = format_pred_status(spip_chk, os.path.exists(spip_vcf))
        pred_branch_str = format_pred_status(branch_chk, os.path.exists(branch_vcf))
        pred_pangolin_str = format_pred_status(pangolin_chk, os.path.exists(pangolin_vcf))

        # Check for active re-run vs older artifact version
        # If any predictor or active log is newer than t_final_pq, or if spliceai is live:
        is_re_running = False
        older_artifact_note = None

        if step_status["final_pq"] and t_final_pq:
            newer_predictors = []
            if t_vep and t_vep > t_final_pq:
                newer_predictors.append("VEP")
            if t_spliceai and t_spliceai > t_final_pq:
                newer_predictors.append("SpliceAI")
            if not spliceai_is_complete:
                newer_predictors.append("SpliceAI (Incomplete/Writing)")
            if t_pangolin and t_pangolin > t_final_pq:
                newer_predictors.append("Pangolin")
            
            # Check gene log mtimes
            gene_log_dir = os.path.join(log_dir, gene_name)
            if os.path.exists(gene_log_dir):
                for lf in glob.glob(os.path.join(gene_log_dir, "*.out")):
                    tl = get_mtime_safe(lf)
                    if tl and tl > t_final_pq + 60: # 1 min buffer
                        log_step = os.path.basename(lf).split(".")[1]
                        if log_step not in ["plot", "report", "filter", "clinrep"] and log_step not in [p.lower() for p in newer_predictors]:
                            newer_predictors.append(f"{log_step.upper()} (Active Log)")

            if newer_predictors or (not spliceai_is_complete):
                is_re_running = True
                t_old_str = datetime.fromtimestamp(t_final_pq).strftime("%Y-%m-%d %H:%M")
                try:
                    df_old = pd.read_parquet(final_pq)
                    n_old_rows = len(df_old)
                except Exception:
                    n_old_rows = 0
                older_artifact_note = f"Older artifact exists ({n_old_rows:,} vars, {t_old_str}) ➔ New run in progress: {', '.join(newer_predictors)}"
                version_lineage.append({
                    "Gene": gene_name,
                    "Previous_Artifact": f"{n_old_rows:,} records ({t_old_str})",
                    "Active_Execution_Status": f"🔄 Re-running: {', '.join(newer_predictors)}",
                    "SpliceAI_Progress": spliceai_detail
                })

        # Determine current active step
        if not step_status["varconv"]:
            cur_step = "1. VarConv"
        elif not (step_status["vep"] and step_status["spip"] and step_status["branchpoint"] and step_status["spliceai"]):
            cur_step = "2. Predictors"
        elif is_re_running or (t_merged and t_final_pq and t_merged > t_final_pq):
            cur_step = "3. Merge / VCF2Parsed"
        elif not step_status["merged"]:
            cur_step = "3. Merge"
        elif not step_status["final_pq"]:
            cur_step = "4. VCF2Parsed"
        else:
            cur_step = "Done"

        # Timing calculation
        t_run_start = t_varconv
        if step_status["final_pq"] and not is_re_running and t_final_pq and t_run_start:
            total_turnaround = t_final_pq - t_run_start
            total_time_str = f"✅ {format_duration(total_turnaround)}"
            completed_at = datetime.fromtimestamp(t_final_pq).strftime("%Y-%m-%d %H:%M:%S")
        elif t_run_start:
            active_turnaround = now - t_run_start
            total_time_str = f"🔄 Running ({format_duration(active_turnaround)})"
            completed_at = f"🔄 Active Re-run" if is_re_running else "In Progress"
        else:
            total_time_str = "-"
            completed_at = "Pending"

        timing_results.append({
            "Gene": gene_name,
            "Variants": f"{n_input_vars:,}" if n_input_vars > 0 else "-",
            "Current_Step": cur_step,
            "SpliceAI_Progress": spliceai_detail,
            "Total_Elapsed": total_time_str,
            "Completed_At": completed_at
        })

        # Check logs for issues
        gene_log_dir = os.path.join(log_dir, gene_name)
        log_issues = []
        if os.path.exists(gene_log_dir):
            for err_file in glob.glob(os.path.join(gene_log_dir, "*.err")):
                if os.path.getsize(err_file) > 0:
                    with open(err_file, 'r', errors='ignore') as ef:
                        txt = ef.read()
                        if any(err_word in txt for err_word in ["ERROR", "Error", "Traceback", "Exception", "fatal", "Fatal", "Killed"]):
                            log_issues.append(os.path.basename(err_file))

        # Check Parquet Completeness & Schemas
        coverage = {}
        n_out_vars = 0
        if step_status["final_pq"]:
            try:
                df_out = pd.read_parquet(final_pq)
                n_out_vars = len(df_out)
                cols_lower = {c.lower(): c for c in df_out.columns}
                all_gene_columns[gene_name] = set(df_out.columns)

                if schema_set is None:
                    schema_set = set(df_out.columns)
                else:
                    if set(df_out.columns) != schema_set:
                        diff = schema_set.symmetric_difference(set(df_out.columns))
                        column_inconsistencies.append((gene_name, len(df_out.columns), list(diff)[:5]))

                for std_name, candidate_cols in key_annotation_cols:
                    matched = None
                    for cand in candidate_cols:
                        if cand.lower() in cols_lower:
                            matched = cols_lower[cand.lower()]
                            break
                    if matched:
                        non_null_cnt = df_out[matched].notna().sum()
                        pct = (non_null_cnt / n_out_vars) * 100 if n_out_vars > 0 else 0
                        coverage[std_name] = f"{pct:.1f}%"
                    else:
                        coverage[std_name] = "N/A"
            except Exception as e:
                log_issues.append(f"ParquetReadError: {str(e)[:30]}")

        if is_re_running:
            final_pq_status = f"🔄 Re-running ({n_out_vars:,} vars in older version)"
            gene_status = "RE-RUNNING"
        elif step_status["final_pq"]:
            final_pq_status = f"✅ {n_out_vars:,} vars"
            gene_status = "FINISHED"
        else:
            final_pq_status = "⏳ Pending"
            gene_status = "RUNNING/PENDING"

        results.append({
            "Gene": gene_name,
            "Status": gene_status,
            "VarConv": "✅" if step_status["varconv"] else "⏳",
            "VEP": pred_vep_str,
            "SPiP": pred_spip_str,
            "Pangolin": pred_pangolin_str,
            "Branchpoint": pred_branch_str,
            "SpliceAI": spliceai_detail,
            "Final_PQ": final_pq_status,
            "Log_Issues": ", ".join(log_issues) if log_issues else "None",
            "Total_Columns": len(all_gene_columns.get(gene_name, [])),
            "Coverage": coverage,
            "Older_Artifact": older_artifact_note
        })

    # Summary KPIs
    tot_genes = len(genes)
    finished_genes = sum(1 for r in results if r["Status"] == "FINISHED")
    pending_genes = tot_genes - finished_genes
    genes_with_errors = sum(1 for r in results if r["Log_Issues"] != "None")
    tot_input_vars = sum(int(t["Variants"].replace(",", "")) for t in timing_results if t["Variants"] != "-")
    
    tot_clean_vars = 0
    for r in results:
        m = re.search(r'([0-9,]+)\s+vars', r["Final_PQ"])
        if m:
            tot_clean_vars += int(m.group(1).replace(",", ""))

    # Build Markdown Report Content
    df_timing = pd.DataFrame(timing_results)
    df_stage = pd.DataFrame([{
        "Gene": r["Gene"],
        "Status": r["Status"],
        "VarConv": r["VarConv"],
        "VEP": r["VEP"],
        "SPiP": r["SPiP"],
        "Pangolin": r["Pangolin"],
        "Branchpoint": r["Branchpoint"],
        "SpliceAI": r["SpliceAI"],
        "Final_PQ": r["Final_PQ"],
        "Log_Issues": r["Log_Issues"]
    } for r in results])

    cov_rows = []
    for r in results:
        if r["Coverage"]:
            c_dict = {"Gene": r["Gene"], "Columns": r["Total_Columns"]}
            c_dict.update(r["Coverage"])
            cov_rows.append(c_dict)
    df_cov = pd.DataFrame(cov_rows) if cov_rows else pd.DataFrame()

    md = []
    md.append(f"# 🧪 Master Pipeline Execution & Quality Audit: `{run_name}`")
    md.append(f"**Audit Timestamp**: `{timestamp}`  ")
    md.append(f"**Target Build**: `GRCh38`  ")
    md.append(f"**Run Directory**: `{run_dir}`  ")
    md.append(f"**Results Directory**: `{res_dir}`  \n")
    md.append("---\n")
    md.append("## 📊 1. Executive Summary")
    md.append(f"- **Total Target Genes**: `{tot_genes}`")
    md.append(f"- **Genes Fully Finished**: `{finished_genes} / {tot_genes}` ({finished_genes/tot_genes*100:.1f}%)")
    md.append(f"- **Genes In Progress / Pending**: `{pending_genes}`")
    md.append(f"- **Genes with Critical Errors**: `{genes_with_errors}`")
    if tot_input_vars > 0:
        md.append(f"- **Total Input Variants**: `{tot_input_vars:,}`")
    md.append(f"- **Total Clean Output Variants Generated**: `{tot_clean_vars:,}`\n")
    md.append("---\n")
    md.append("## ⏱️ 2. Per-Gene Step Tracking & Running Times")
    md.append(df_timing.to_markdown(index=False))
    md.append("\n---\n")
    md.append("## 📋 3. Stage-by-Stage Predictor Matrix")
    md.append(df_stage.to_markdown(index=False))
    md.append("\n---\n")
    md.append("## 🧬 4. Key Predictor Completeness (% Non-Null in Output Parquets)")
    if not df_cov.empty:
        md.append(df_cov.to_markdown(index=False))
    else:
        md.append("_No output Parquet files generated yet to assess completeness._")
    md.append("\n---\n")
    md.append("## 📐 5. Schema & Column Integrity")
    if all_gene_columns:
        col_counts = [len(v) for v in all_gene_columns.values()]
        md.append(f"- **Column Count Range**: `{min(col_counts)}` to `{max(col_counts)}` columns per parquet table.")
        if not column_inconsistencies:
            md.append("- **Schema Consistency**: ✅ **100% Unified Schema** across all finished genes!")
        else:
            md.append("- **Schema Consistency**: ⚠️ **Inconsistency detected** across genes!")
            for g, count, diff_cols in column_inconsistencies[:5]:
                md.append(f"  - `{g}` diff columns ({len(diff_cols)}): {diff_cols}")
    else:
        md.append("_No output Parquet files available to assess schema integrity._")
    md.append("\n")

    if version_lineage:
        df_lineage = pd.DataFrame(version_lineage)
        md.append("---\n")
        md.append("## 🔄 6. Active Re-runs & Older Artifact Lineage")
        md.append("> [!IMPORTANT]")
        md.append("> The following genes have active ongoing executions or newly generated predictor data that postdate older parquet artifacts:")
        md.append(df_lineage.to_markdown(index=False))
        md.append("\n")


    report_content = "\n".join(md)

    # Save Markdown Report strictly in run reports folder
    with open(output_report_path, "w") as f:
        f.write(report_content)
    logger.info(f"Master audit report saved to: {output_report_path}")

    # Print Terminal Summary
    print("\n" + "=" * 95)
    print(f"  🧪 Master Pipeline Quality & Completeness Audit Summary ({run_name})")
    print("=" * 95)
    print(df_timing.to_string(index=False))
    print("=" * 95)
    print(f"Finished Genes: {finished_genes} / {tot_genes} ({finished_genes/tot_genes*100:.1f}%) | Pending: {pending_genes} | Errors: {genes_with_errors}")
    print(f"Report Saved  : {output_report_path}")
    print("=" * 95 + "\n")

    return df_timing

audit_single_run = audit_run

def main():
    args = parse_args()
    
    if args.all_runs:
        available_runs = discover_available_runs()
        if not available_runs:
            logger.error("No run folders found in RUNS/")
            sys.exit(1)
        for rdir in available_runs:
            audit_single_run(run_dir=rdir, raw_dir=args.raw_dir)
        sys.exit(0)

    target_run = None
    if args.run_dir:
        target_run = args.run_dir
    elif args.run_name:
        target_run = os.path.join("RUNS", args.run_name)
    else:
        available_runs = discover_available_runs()
        if available_runs:
            target_run = available_runs[0]
        else:
            logger.error("No run directories found in RUNS/")
            sys.exit(1)

    audit_single_run(run_dir=target_run, raw_dir=args.raw_dir, output_report_path=args.output_report)

if __name__ == "__main__":
    main()
