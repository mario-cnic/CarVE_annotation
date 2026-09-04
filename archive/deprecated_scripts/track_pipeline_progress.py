#!/usr/bin/env python3
"""
Comprehensive Multi-Pipeline Progress Tracker and Health Auditor
Targeting:
  - RUNS/run_20260813_1028 (Cardiomyopathy Benchmark Cohort - 7 Core Genes)
  - RUNS/run_more_genes_20260817_1143 (Extended Cardiomyopathy & Channelopathy Cohort - 118 Genes)

Audits:
  1. Input variants vs output clean Parquet records per gene.
  2. Granular step-by-step progress (VarConv, VEP, SPiP, Pangolin, Branchpoint, SpliceAI, Merge, Parsed PQ).
  3. SpliceAI deep learning forward-pass progress and ETA calculation.
  4. Schema consistency (column counts, key predictor coverage %).
  5. Error/warning diagnostics from _log/ subdirectories.
  6. Overall completion percentages and active bottleneck detection.
"""

import os
import sys
import glob
import re
import subprocess
from datetime import datetime
import pandas as pd
import numpy as np

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
            cmd = f"zcat '{vcf_path}' | grep -v '^#' | wc -l"
        else:
            cmd = f"grep -v '^#' '{vcf_path}' | wc -l"
        return int(subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL).decode().strip())
    except Exception:
        return 0

def check_spliceai_status(run_dir, gene):
    ann_dir = os.path.join(run_dir, "annotation")
    tmp_dir = os.path.join(run_dir, "_tmp")
    log_dir = os.path.join(run_dir, "_log", gene)
    
    final_vcf_gz = os.path.join(ann_dir, f"{gene}.annSpliceAI.vcf.gz")
    final_vcf = os.path.join(ann_dir, f"{gene}.annSpliceAI.vcf")
    
    if os.path.exists(final_vcf_gz) and os.path.getsize(final_vcf_gz) > 1000:
        sz_mb = os.path.getsize(final_vcf_gz) / (1024 * 1024)
        return {"status": "DONE", "detail": f"✅ Done ({sz_mb:.2f} MB gz)"}
    
    if os.path.exists(final_vcf) and os.path.getsize(final_vcf) > 1000:
        sz_mb = os.path.getsize(final_vcf) / (1024 * 1024)
        return {"status": "DONE_RAW", "detail": f"✅ Done raw ({sz_mb:.2f} MB uncompressed)"}
    
    # Check parallel chunking directory
    chunk_dirs = glob.glob(os.path.join(tmp_dir, f"{gene}_spliceai_chunks*")) + \
                 glob.glob(os.path.join(tmp_dir, f"spliceai_chunks_{gene}*"))
    if chunk_dirs:
        cdir = chunk_dirs[0]
        chunks_all = glob.glob(os.path.join(cdir, "chunk_*.vcf.gz"))
        chunks_done = glob.glob(os.path.join(cdir, "chunk_*.annSpliceAI.vcf.gz"))
        if chunks_all:
            pct = (len(chunks_done) / len(chunks_all)) * 100
            return {"status": "IN_PROGRESS", "detail": f"🔄 Chunking: {len(chunks_done)}/{len(chunks_all)} chunks ({pct:.1f}%)"}

    # Check stdout logs for pass counts
    out_log = os.path.join(log_dir, f"{gene}.spliceai.out")
    vcf_tmp = os.path.join(tmp_dir, f"{gene}.vcf.gz")
    if not os.path.exists(vcf_tmp):
        vcf_tmp = os.path.join(tmp_dir, f"{gene}.vcf")
        
    if os.path.exists(out_log) and os.path.exists(vcf_tmp):
        try:
            tot_vars = count_vcf_variants(vcf_tmp)
            tot_passes = tot_vars * 5  # 5 models
            cmd_passes = f"grep -c '1/1 \\[' '{out_log}' 2>/dev/null || echo 0"
            done_passes = int(subprocess.check_output(cmd_passes, shell=True).decode().strip())
            
            pct = (done_passes / tot_passes) * 100 if tot_passes > 0 else 0
            rem_passes = max(0, tot_passes - done_passes)
            rem_hours = (rem_passes * 3.0) / 3600
            if pct >= 99.9 or (tot_vars > 0 and done_passes >= tot_passes):
                return {"status": "DONE", "detail": "✅ Finished passes"}
            elif done_passes > 0:
                return {"status": "IN_PROGRESS", "detail": f"🔄 {done_passes:,}/{tot_passes:,} passes ({pct:.1f}% | ~{rem_hours:.1f}h left)"}
        except Exception:
            pass

    # Check if process is active or queued
    if os.path.exists(out_log) and os.path.getsize(out_log) > 0:
        return {"status": "IN_PROGRESS", "detail": "🔄 Active Inference"}

    return {"status": "PENDING", "detail": "⏳ Pending"}

def audit_single_run(run_dir):
    run_name = os.path.basename(run_dir.rstrip("/"))
    ann_dir = os.path.join(run_dir, "annotation")
    res_dir = os.path.join(run_dir, "results")
    tmp_dir = os.path.join(run_dir, "_tmp")
    log_dir = os.path.join(run_dir, "_log")
    
    # Discover genes
    genes = set()
    for f in glob.glob(os.path.join(tmp_dir, "*.vcf*")):
        base = os.path.basename(f).split(".")[0].split("_")[0]
        if base and not base.startswith("chunk"):
            genes.add(base)
    for f in glob.glob(os.path.join(res_dir, "*.clean.pq")):
        base = os.path.basename(f).split(".")[0].split("_")[0]
        genes.add(base)
    for f in glob.glob(os.path.join(ann_dir, "*.annVEP.vcf.gz")):
        base = os.path.basename(f).split(".")[0].split("_")[0]
        genes.add(base)
    for d in glob.glob(os.path.join(log_dir, "*")):
        if os.path.isdir(d):
            genes.add(os.path.basename(d))
            
    genes = sorted(list(genes))
    
    gene_records = []
    schema_cols = {}
    
    key_annotation_cols = [
        ("gnomADv4_AF_grpmax", ["gnomadv4_af_grpmax_joint", "gnomadv4_af_grpmax", "gnomad_af_joint", "max_af"]),
        ("REVEL_score", ["revel_score", "revel"]),
        ("AlphaMissense", ["am_pathogenicity", "alphamissense_score", "alphamissense"]),
        ("SPiP", ["spip_prediction", "spip", "spip_interpretation"]),
        ("Pangolin", ["pangolin_max_score", "pangolin_largest_delta", "pangolin"]),
        ("SpliceAI", ["spliceai_custom_max", "spliceai_max", "spliceai_pred_ds_max"]),
        ("Branchpoint", ["branchpoint_disrupted", "branchpoint_status", "labranchor_score"]),
        ("5UTR_Annotator", ["5utr_consequence", "5utr_annotation"]),
        ("SpliceVault", ["splicevault_status", "splicevault_top_events"]),
        ("MaxEntScan", ["maxentscan_diff", "maxentscan_ref"]),
        ("Priority_Tier", ["priority_tier"]),
        ("Priority_Score", ["variant_priority_score"])
    ]
    
    for gene in genes:
        vcf_tmp = os.path.join(tmp_dir, f"{gene}.vcf.gz")
        if not os.path.exists(vcf_tmp):
            vcf_tmp = os.path.join(tmp_dir, f"{gene}.vcf")
            
        vep_vcf = os.path.join(ann_dir, f"{gene}.annVEP.vcf.gz")
        spip_vcf = os.path.join(ann_dir, f"{gene}.annSPiP.vcf.gz")
        if not os.path.exists(spip_vcf):
            spip_vcf = os.path.join(ann_dir, f"{gene}.annSPiP.vcf")
        branch_vcf = os.path.join(ann_dir, f"{gene}.annBranchpoint.vcf.gz")
        pangolin_vcf = os.path.join(ann_dir, f"{gene}.annPangolin.vcf.gz")
        spliceai_vcf = os.path.join(ann_dir, f"{gene}.annSpliceAI.vcf.gz")
        if not os.path.exists(spliceai_vcf):
            spliceai_vcf = os.path.join(ann_dir, f"{gene}.annSpliceAI.vcf")
        merged_vcf = os.path.join(ann_dir, f"{gene}.annotated.vcf.gz")
        final_pq = os.path.join(res_dir, f"{gene}.parsed.clean.pq")
        final_tsv = os.path.join(res_dir, f"{gene}.parsed.tsv")
        
        has_varconv = os.path.exists(vcf_tmp) and os.path.getsize(vcf_tmp) > 0
        n_input_vars = count_vcf_variants(vcf_tmp) if has_varconv else 0
        
        has_vep = os.path.exists(vep_vcf) and os.path.getsize(vep_vcf) > 1000
        has_spip = os.path.exists(spip_vcf) and os.path.getsize(spip_vcf) > 100
        has_branch = os.path.exists(branch_vcf) and os.path.getsize(branch_vcf) > 100
        has_pangolin = os.path.exists(pangolin_vcf) and os.path.getsize(pangolin_vcf) > 1000
        
        spliceai_info = check_spliceai_status(run_dir, gene)
        has_spliceai = spliceai_info["status"] in ["DONE", "DONE_RAW"]
        
        has_merged = os.path.exists(merged_vcf) and os.path.getsize(merged_vcf) > 1000
        has_final_pq = os.path.exists(final_pq) and os.path.getsize(final_pq) > 0
        
        n_output_records = 0
        n_cols = 0
        ann_coverages = {}
        
        if has_final_pq:
            try:
                df = pd.read_parquet(final_pq)
                n_output_records = len(df)
                n_cols = len(df.columns)
                schema_cols[gene] = list(df.columns)
                
                for label, cands in key_annotation_cols:
                    found = None
                    for c in df.columns:
                        if c.lower() in cands or c in cands:
                            found = c
                            break
                    if found:
                        valid = df[found].dropna()
                        valid = valid[~valid.astype(str).str.strip().isin(['', 'nan', '.', '-', 'NA', 'none', 'None', 'not_covered'])]
                        pct = (len(valid) / max(1, n_output_records)) * 100
                        ann_coverages[label] = f"{pct:.1f}%"
                    else:
                        ann_coverages[label] = "❌ Missing"
            except Exception as e:
                ann_coverages["ERROR"] = str(e)
                
        # Determine active stage
        if has_final_pq:
            cur_stage = "✅ Complete"
        elif not has_varconv:
            cur_stage = "1. VarConv"
        elif not (has_vep and has_spip and has_branch and has_spliceai):
            sub_pending = []
            if not has_vep: sub_pending.append("VEP")
            if not has_spip: sub_pending.append("SPiP")
            if not has_branch: sub_pending.append("Branch")
            if not has_spliceai: sub_pending.append("SpliceAI")
            cur_stage = f"2. Predictors ({'/'.join(sub_pending)})"
        elif not has_merged:
            cur_stage = "3. Merge"
        else:
            cur_stage = "4. VCF2Parsed"
            
        # Log inspection
        g_log = os.path.join(log_dir, gene)
        log_errors = []
        if os.path.exists(g_log):
            for ef in glob.glob(os.path.join(g_log, "*.err")):
                if os.path.getsize(ef) > 0:
                    with open(ef, "r", errors="ignore") as f:
                        c = f.read()
                        if any(kw in c.lower() for kw in ["exception", "command not found", "killed", "stale file handle", "filenotfounderror", "segmentation fault"]):
                            log_errors.append(os.path.basename(ef))

        gene_records.append({
            "Gene": gene,
            "Input_Vars": n_input_vars,
            "Output_Records": n_output_records,
            "Stage": cur_stage,
            "VarConv": "✅" if has_varconv else "❌",
            "VEP": "✅" if has_vep else "❌",
            "SPiP": "✅" if has_spip else "❌",
            "Pangolin": "✅" if has_pangolin else "❌",
            "Branchpoint": "✅" if has_branch else "❌",
            "SpliceAI": spliceai_info["detail"],
            "Merged": "✅" if has_merged else "❌",
            "Final_PQ": f"✅ {n_output_records:,} rows ({n_cols} cols)" if has_final_pq else "⏳ Pending",
            "gnomADv4": ann_coverages.get("gnomADv4_AF_grpmax", "-"),
            "REVEL": ann_coverages.get("REVEL_score", "-"),
            "AlphaMissense": ann_coverages.get("AlphaMissense", "-"),
            "SPiP_cov": ann_coverages.get("SPiP", "-"),
            "Branch_cov": ann_coverages.get("Branchpoint", "-"),
            "SpliceAI_cov": ann_coverages.get("SpliceAI", "-"),
            "5UTR_cov": ann_coverages.get("5UTR_Annotator", "-"),
            "MaxEnt_cov": ann_coverages.get("MaxEntScan", "-"),
            "Tier_cov": ann_coverages.get("Priority_Tier", "-"),
            "Score_cov": ann_coverages.get("Priority_Score", "-"),
            "Log_Issues": ", ".join(log_errors) if log_errors else "None"
        })
        
    df_all = pd.DataFrame(gene_records)
    return {
        "run_name": run_name,
        "run_dir": run_dir,
        "df": df_all,
        "schema_cols": schema_cols,
        "n_genes": len(genes),
        "n_completed": len(df_all[df_all["Stage"] == "✅ Complete"]),
        "total_in_vars": df_all["Input_Vars"].sum(),
        "total_out_records": df_all["Output_Records"].sum()
    }

def main():
    target_runs = [
        "/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_20260813_1028",
        "/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/RUNS/run_more_genes_20260817_1143"
    ]
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S CEST")
    report_lines = []
    
    report_lines.append("# 📊 Comprehensive Multi-Pipeline Progress & Execution Audit")
    report_lines.append(f"**Audit Execution Timestamp**: `{timestamp}`  ")
    report_lines.append("**Target Environment**: Local Sandbox & HPC SGE Cluster State  ")
    report_lines.append("**Genome Build**: `GRCh38` (Ensembl v111 / GATK Bundle v0)  \n")
    report_lines.append("---")
    
    audit_summaries = []
    
    for rdir in target_runs:
        if not os.path.exists(rdir):
            continue
        audit = audit_single_run(rdir)
        audit_summaries.append(audit)
        
        df = audit["df"]
        n_tot = audit["n_genes"]
        n_comp = audit["n_completed"]
        pct_comp = (n_comp / n_tot) * 100 if n_tot > 0 else 0
        
        report_lines.append(f"## 🚀 Pipeline Run: `{audit['run_name']}`")
        report_lines.append(f"- **Path**: `{audit['run_dir']}`")
        report_lines.append(f"- **Progress**: `{n_comp} / {n_tot}` genes complete (**{pct_comp:.1f}%**)")
        report_lines.append(f"- **Total Input Variants Processed**: `{audit['total_in_vars']:,}`")
        report_lines.append(f"- **Total Output Consequence Records**: `{audit['total_out_records']:,}`")
        
        # Breakdown of stages
        stage_counts = df["Stage"].value_counts().to_dict()
        report_lines.append(f"- **Stage Breakdown**:")
        for stg, cnt in stage_counts.items():
            report_lines.append(f"  - `{stg}`: **{cnt}** genes")
            
        report_lines.append("\n### Per-Gene Execution Status Matrix")
        summary_cols = ["Gene", "Input_Vars", "Stage", "VarConv", "VEP", "SPiP", "Pangolin", "Branchpoint", "SpliceAI", "Merged", "Final_PQ", "Log_Issues"]
        report_lines.append(df[summary_cols].to_markdown(index=False))
        report_lines.append("")
        
        # Predictor coverage for completed genes
        df_comp = df[df["Stage"] == "✅ Complete"]
        if len(df_comp) > 0:
            report_lines.append("### Key Predictor Coverage in Completed Parquet Files")
            cov_cols = ["Gene", "Output_Records", "gnomADv4", "REVEL", "AlphaMissense", "SPiP_cov", "Branch_cov", "SpliceAI_cov", "5UTR_cov", "MaxEnt_cov", "Tier_cov", "Score_cov"]
            report_lines.append(df_comp[cov_cols].to_markdown(index=False))
            report_lines.append("")
            
        # Schema integrity
        if audit["schema_cols"]:
            col_lens = {g: len(cols) for g, cols in audit["schema_cols"].items()}
            min_c = min(col_lens.values())
            max_c = max(col_lens.values())
            report_lines.append(f"### Schema Consistency")
            if min_c == max_c:
                report_lines.append(f"- ✅ **100% Consistent Schema**: All {len(col_lens)} parquet files have exactly **`{min_c}` columns**.")
            else:
                report_lines.append(f"- ⚠️ **Schema Variance**: Columns range from `{min_c}` to `{max_c}` across genes.")
                for g, l in col_lens.items():
                    if l != max_c:
                        report_lines.append(f"  - `{g}`: {l} cols")
        report_lines.append("\n" + "="*80 + "\n")

    # Generate Cross-Run Executive Overview
    overview_lines = [
        "## 📌 Executive Cross-Pipeline Summary",
        f"| Pipeline Run | Total Genes | Completed | In Progress / Pending | Total Input Variants | Clean Output Records | Schema Status |",
        f"| :--- | :---: | :---: | :---: | :---: | :---: | :---: |"
    ]
    for s in audit_summaries:
        n_p = s["n_genes"] - s["n_completed"]
        pct = (s["n_completed"] / s["n_genes"]) * 100 if s["n_genes"] > 0 else 0
        col_lens = [len(v) for v in s["schema_cols"].values()] if s["schema_cols"] else [0]
        schema_st = f"✅ Unified ({col_lens[0]} cols)" if len(set(col_lens)) == 1 and col_lens[0] > 0 else f"⚠️ {min(col_lens)}-{max(col_lens)} cols"
        overview_lines.append(
            f"| **`{s['run_name']}`** | `{s['n_genes']}` | **`{s['n_completed']}`** ({pct:.1f}%) | `{n_p}` | `{s['total_in_vars']:,}` | `{s['total_out_records']:,}` | {schema_st} |"
        )
    overview_lines.append("")

    full_report = "\n".join(report_lines[:5] + overview_lines + report_lines[5:])
    
    # Save report to walkthrough
    walkthrough_out = "/home/mruizp/data_lab_PGP/pipelines/annotation_pipeline_new/walkthrough/20260819_dual_pipeline_execution_progress_audit.md"
    with open(walkthrough_out, "w") as f:
        f.write(full_report)
        
    print(f"Audit successfully generated and saved to: {walkthrough_out}")

if __name__ == "__main__":
    main()
