#!/usr/bin/env python3
"""
eda_parsed_pq.py
Exploratory Data Analysis (EDA) utility for ".parsed.clean.pq" Parquet files.
Analyzes transcript mappings, variant classes, sources, allele frequencies, and pathogenicity scores.
Produces a detailed Markdown report.
"""

import os
import glob
import re
import argparse
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Perform EDA on pipeline Parquet results.')
    parser.add_argument('--input', type=str, default='results',
                        help='Input file path (.pq) or directory containing .parsed.clean.pq files.')
    parser.add_argument('--output', type=str, default='results/parquet_eda_report.md',
                        help='Output path for the Markdown report (set to "stdout" to print to terminal).')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output.')
    parser.add_argument('--raw-dir', type=str, default=None,
                        help='Raw directory containing input files to audit execution status and coverage.')
    parser.add_argument('--log-dir', type=str, default=None,
                        help='Log directory containing cluster execution logs (defaults to _log/<raw_folder_name>).')
    return parser.parse_args()

def get_gene_pipeline_status(gene, log_dir):
    """
    Scans the gene's log directory to determine why it failed.
    Returns (status, last_failed_step, error_message).
    """
    if not log_dir or not os.path.exists(log_dir):
        return "Failed", "Unknown", "Log directory not found"
        
    gene_log_dir = os.path.join(log_dir, gene)
    if not os.path.exists(gene_log_dir):
        return "Failed / Unstarted", "Unknown", "No log directory found (job likely unstarted or pending in queue)"
        
    pipeline_report_path = os.path.join(gene_log_dir, "pipeline_report.md")
    failed_step = "Unknown"
    exit_code = "Unknown"
    
    if os.path.exists(pipeline_report_path):
        try:
            with open(pipeline_report_path, "r") as f:
                content = f.read()
            
            # Find last failed step in pipeline_report.md
            sections = content.split("### Step: ")
            for section in sections[1:]:
                lines = section.split("\n")
                step_name = lines[0].strip()
                for l in lines:
                    if "Failed" in l:
                        match = re.search(r"Failed\s*\(Exit\s*Code\s*(\d+)\)", l)
                        failed_step = step_name
                        exit_code = match.group(1) if match else "Error"
        except Exception:
            pass
            
    # Try to find the exact error message from the most recent non-empty err file
    err_msg = "No error message captured"
    try:
        err_files = glob.glob(os.path.join(gene_log_dir, "*.err"))
        if err_files:
            err_files.sort(key=os.path.getmtime, reverse=True)
            for err_file in err_files:
                if os.path.getsize(err_file) > 0:
                    with open(err_file, "r") as f:
                        lines = [l.strip() for l in f.readlines() if l.strip()]
                    if lines:
                        err_msg = lines[-1]
                        if "JOB HAD KILL SIGNAL" in err_msg:
                            err_msg = "Job was killed by SGE (OOM, wall-time limit, or qdel)"
                        elif "Dropped because it is disabled" in err_msg or "temporarily not available" in err_msg:
                            if len(lines) > 1:
                                err_msg = lines[-2]
                        break
    except Exception:
        pass
        
    if failed_step != "Unknown":
        status_str = f"Failed at {failed_step} (Exit Code {exit_code})"
    else:
        status_str = "Failed"
        
    return status_str, failed_step, err_msg

def make_ascii_histogram(series, bins=10):
    """Generates an ASCII-based distribution bar chart for allele frequencies."""
    if series.empty:
        return "  No allele frequency data available."
    counts, edges = np.histogram(series, bins=bins)
    max_count = max(counts) if max(counts) > 0 else 1
    max_width = 30
    hist_lines = []
    for i in range(bins):
        width = int((counts[i] / max_count) * max_width)
        bar = "█" * width + "░" * (max_width - width)
        hist_lines.append(f"  [{edges[i]:.2e} - {edges[i+1]:.2e}]: {bar} {counts[i]:,}")
    return "\n".join(hist_lines)

def analyze_single_file(file_path, verbose=False):
    if verbose:
        print(f"Reading {file_path}...")
    
    try:
        df = pd.read_parquet(file_path)
    except Exception as e:
        return {"error": f"Failed to read file: {e}"}
        
    n_rows, n_cols = df.shape
    file_size_kb = os.path.getsize(file_path) / 1024
    
    # Extract Genes
    genes = set()
    for col in ['SYMBOL', 'SYMBOL_VEP', 'GEN', 'gene']:
        if col in df.columns:
            genes.update(df[col].dropna().unique())
    gene_list = sorted([str(g) for g in genes if g])
    
    # Extract Target Transcript from Feature column (VEP filtered transcript)
    target_transcripts = set()
    if 'Feature' in df.columns:
        for val in df['Feature'].dropna().unique():
            for part in str(val).split('&'):
                target_transcripts.add(part)
    target_list = sorted([str(t) for t in target_transcripts if t])
    
    # Extract Associated Transcripts (split from Ensembl_transcriptid dbNSFP annotations)
    associated_transcripts = set()
    for col in ['Ensembl_transcriptid', 'ensembl_transcript_id']:
        if col in df.columns:
            for val in df[col].dropna().unique():
                for part in str(val).split('&'):
                    associated_transcripts.add(part)
    associated_list = sorted([str(t) for t in associated_transcripts if t])
    
    # Extract RefSeq Transcripts (split NM IDs)
    refseq_transcripts = set()
    for col in ['HGVSc', 'HGVSc_VEP', 'Refseq_transcript', 'NM']:
        if col in df.columns:
            vals = df[col].dropna().unique()
            for v in vals:
                for part in str(v).split('&'):
                    match = re.search(r'(NM_\d+(?:\.\d+)?)', part)
                    if match:
                        refseq_transcripts.add(match.group(1))
                    elif re.match(r'^NM_\d+', part):
                        refseq_transcripts.add(part)
    nm_list = sorted(list(refseq_transcripts))
    
    # Variant Class distribution
    variant_classes = {}
    for col in ['VARIANT_CLASS', 'TIPO_DE_VARIANTE', 'variant_class']:
        if col in df.columns:
            variant_classes = df[col].value_counts().to_dict()
            break
            
    # Source Origin distribution
    source_origins = {}
    for col in ['SOURCE_ORIGIN', 'source_origin']:
        if col in df.columns:
            source_origins = df[col].value_counts().to_dict()
            break
            
    # Impact level distribution
    impacts = {}
    for col in ['IMPACT', 'IMPACT_VEP', 'impact']:
        if col in df.columns:
            impacts = df[col].value_counts().to_dict()
            break
            
    # Impact breakdown counts for Summary table
    impact_counts = {'HIGH': 0, 'MODERATE': 0, 'LOW': 0, 'MODIFIER': 0}
    if impacts:
        for k, v in impacts.items():
            k_upper = str(k).upper()
            impact_counts[k_upper] = v
            
    # Splicing predictions (SPiP)
    spip_summary = {}
    for col in ['SPiP_prediction', 'SPiP_interpretation', 'spip_prediction', 'SPiP']:
        if col in df.columns:
            preds = df[col].dropna()
            if not preds.empty:
                spip_summary = preds.value_counts().head(5).to_dict()
            break
            
    af_col = 'gnomADv4_AF_grpmax_joint' if 'gnomADv4_AF_grpmax_joint' in df.columns else None
                
    af_stats = {}
    if af_col:
        freqs = pd.to_numeric(df[af_col], errors='coerce').dropna()
        if not freqs.empty:
            af_stats = {
                "column": af_col,
                "mean": freqs.mean(),
                "median": freqs.median(),
                "max": freqs.max(),
                "min": freqs.min(),
                "missing_pct": (df[af_col].isna().sum() / n_rows) * 100,
                "raw_values": freqs
            }
            
            # Mean frequency by origin
            if 'SOURCE_ORIGIN' in df.columns:
                origin_means = {}
                for origin, sub_df in df.groupby('SOURCE_ORIGIN'):
                    sub_freqs = pd.to_numeric(sub_df[af_col], errors='coerce').dropna()
                    if not sub_freqs.empty:
                        origin_means[str(origin)] = sub_freqs.mean()
                af_stats["by_origin_mean"] = origin_means
    # Additional population genetics metrics
    pop_stats = {}
    
    # 1. FAF (Filter Allele Frequency)
    faf_col = None
    for col in ['gnomADv4_fafmax_faf95_max_joint', 'gnomADv4_faf95_joint', 'gnomADv4_faf95_joint_nfe']:
        if col in df.columns:
            faf_col = col
            break
            
    if faf_col:
        faf_vals = pd.to_numeric(df[faf_col], errors='coerce').dropna()
        if not faf_vals.empty:
            pop_stats["faf"] = {
                "column": faf_col,
                "max": faf_vals.max(),
                "non_zero_count": int((faf_vals > 0).sum()),
                "total": len(faf_vals)
            }
            
    # 2. NFE vs Non-NFE Allele Counts and Frequencies
    ac_joint_col = 'gnomADv4_AC_joint'
    an_joint_col = 'gnomADv4_AN_joint'
    ac_nfe_col = 'gnomADv4_AC_joint_nfe'
    an_nfe_col = 'gnomADv4_AN_joint_nfe'
    
    if all(col in df.columns for col in [ac_joint_col, an_joint_col, ac_nfe_col, an_nfe_col]):
        ac_joint = pd.to_numeric(df[ac_joint_col], errors='coerce').fillna(0)
        an_joint = pd.to_numeric(df[an_joint_col], errors='coerce').fillna(0)
        ac_nfe = pd.to_numeric(df[ac_nfe_col], errors='coerce').fillna(0)
        an_nfe = pd.to_numeric(df[an_nfe_col], errors='coerce').fillna(0)
        
        ac_non_nfe = ac_joint - ac_nfe
        an_non_nfe = an_joint - an_nfe
        ac_non_nfe = ac_non_nfe.clip(lower=0)
        an_non_nfe = an_non_nfe.clip(lower=0)
        
        sum_ac_nfe = ac_nfe.sum()
        sum_an_nfe = an_nfe.sum()
        sum_ac_non_nfe = ac_non_nfe.sum()
        sum_an_nfe_total = an_non_nfe.sum()
        
        af_nfe_avg = (sum_ac_nfe / sum_an_nfe) if sum_an_nfe > 0 else 0
        af_non_nfe_avg = (sum_ac_non_nfe / sum_an_nfe_total) if sum_an_nfe_total > 0 else 0
        
        pop_stats["nfe_comparison"] = {
            "nfe_ac": int(sum_ac_nfe),
            "nfe_an": int(sum_an_nfe),
            "nfe_af": af_nfe_avg,
            "non_nfe_ac": int(sum_ac_non_nfe),
            "non_nfe_an": int(sum_an_nfe_total),
            "non_nfe_af": af_non_nfe_avg
        }
                
    # Pathogenicity Predictors: REVEL & AlphaMissense
    predictors = {}
    
    # REVEL
    revel_col = None
    for col in ['REVEL', 'REVEL_score', 'revel_score', 'revel']:
        if col in df.columns:
            revel_col = col
            break
    if revel_col:
        rev_scores = pd.to_numeric(df[revel_col], errors='coerce').dropna()
        if not rev_scores.empty:
            predictors["REVEL"] = {
                "mean": rev_scores.mean(),
                "max": rev_scores.max(),
                "high_count": int((rev_scores > 0.5).sum()),
                "total": len(rev_scores)
            }
            
    # AlphaMissense
    am_col = None
    for col in ['AlphaMissense_score', 'am_score', 'AlphaMissense']:
        if col in df.columns:
            am_col = col
            break
    if am_col:
        am_scores = pd.to_numeric(df[am_col], errors='coerce').dropna()
        if not am_scores.empty:
            predictors["AlphaMissense"] = {
                "mean": am_scores.mean(),
                "max": am_scores.max(),
                "high_count": int((am_scores > 0.564).sum()),
                "total": len(am_scores)
            }

    return {
        "file_name": os.path.basename(file_path),
        "file_size_kb": file_size_kb,
        "rows": n_rows,
        "cols": n_cols,
        "genes": gene_list,
        "target_transcripts": target_list,
        "associated_transcripts": associated_list,
        "refseq": nm_list,
        "variant_classes": variant_classes,
        "source_origins": source_origins,
        "impacts": impacts,
        "impact_counts": impact_counts,
        "spip": spip_summary,
        "af": af_stats,
        "pop_stats": pop_stats,
        "predictors": predictors
    }

def format_markdown_report(results, audit_info=None):
    lines = []
    lines.append("# Parquet Outlier & Annotation EDA Report")
    lines.append(f"*Generated automatically on genomic datasets.*")
    lines.append("")
    
    if audit_info:
        lines.append("## 📋 Pipeline Execution & Coverage Audit")
        lines.append("")
        lines.append(f"- **Total Genes in Input Cohort**: {audit_info['total']}")
        lines.append(f"- **Successfully Processed (Parquet Generated)**: {audit_info['success_count']} / {audit_info['total']} ({audit_info['success_pct']:.1f}%)")
        lines.append(f"- **Failed / Missing Genes**: {audit_info['failed_count']} / {audit_info['total']} ({audit_info['failed_pct']:.1f}%)")
        lines.append("")
        
        if audit_info['failed_genes']:
            lines.append("### ❌ Details of Failed / Missing Genes")
            lines.append("")
            lines.append("| Gene | Pipeline Status | Primary Error Message / Rationale |")
            lines.append("| :--- | :--- | :--- |")
            for g in sorted(audit_info['failed_genes']):
                status_str, _, err_msg = audit_info['failed_details'][g]
                lines.append(f"| **{g}** | {status_str} | `{err_msg}` |")
            lines.append("")
        lines.append("---")
        lines.append("")
    
    lines.append("> [!IMPORTANT]")
    lines.append("> **Transcript Filtering Guide**:")
    lines.append("> - Always filter by transcript using the **`Feature`** column. This column contains the specific, single transcript ID for which the row's VEP annotations (like `Consequence` and `IMPACT`) were filtered.")
    lines.append("> - Do **NOT** filter using `Ensembl_transcriptid` or `ensembl_transcript_id` columns from database sources like dbNSFP, as they contain `&`-concatenated lists of all transcripts that overlap the variant coordinates, which causes identical variant lists to be returned across different isoforms.")
    lines.append("")
    
    lines.append("## 📊 Executive Summary Table")
    lines.append("")
    lines.append("| File Name | Genes | Rows | Target Transcript | HIGH | MODERATE | LOW | MODIFIER | Associated ENSTs | Size (KB) |")
    lines.append("| :--- | :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |")
    for r in results:
        if "error" in r:
            lines.append(f"| {r['file_name']} | *Error* | - | - | - | - | - | - | - | - |")
            continue
        genes_str = ", ".join(r['genes'][:3]) + ("..." if len(r['genes']) > 3 else "")
        target_str = ", ".join(r['target_transcripts'][:2]) + ("..." if len(r['target_transcripts']) > 2 else "")
        if not target_str:
            target_str = "None"
        
        ic = r['impact_counts']
        high_c = ic.get('HIGH', 0)
        mod_c = ic.get('MODERATE', 0)
        low_c = ic.get('LOW', 0)
        modifier_c = ic.get('MODIFIER', 0)
        assoc_count = len(r['associated_transcripts'])
        
        lines.append(f"| {r['file_name']} | {genes_str} | {r['rows']:,} | {target_str} | {high_c:,} | {mod_c:,} | {low_c:,} | {modifier_c:,} | {assoc_count} | {r['file_size_kb']:.1f} |")
    lines.append("")
    
    for r in results:
        if "error" in r:
            lines.append(f"### ❌ Error analyzing {r['file_name']}")
            lines.append(f"```{r['error']}```\n")
            continue
            
        lines.append(f"## 🔍 File Details: `{r['file_name']}`")
        lines.append("")
        
        lines.append("### 🧬 Transcripts & Gene Mapping")
        lines.append(f"- **Gene Symbol(s)**: {', '.join(r['genes']) if r['genes'] else 'Unknown'}")
        lines.append(f"- **Target Transcript (from VEP `Feature`)**: {', '.join(r['target_transcripts']) if r['target_transcripts'] else 'None'}")
        
        if r['associated_transcripts']:
            if len(r['associated_transcripts']) > 12:
                assoc_str = ", ".join(r['associated_transcripts'][:12]) + f" ... (+{len(r['associated_transcripts'])-12} more)"
            else:
                assoc_str = ", ".join(r['associated_transcripts'])
            lines.append(f"- **Unique Associated Transcripts (from dbNSFP/annotations)**: {assoc_str} (Total: {len(r['associated_transcripts'])})")
        else:
            lines.append("- **Unique Associated Transcripts (from dbNSFP/annotations)**: None")
            
        lines.append(f"- **RefSeq Transcript(s)**: {', '.join(r['refseq']) if r['refseq'] else 'None mapped'}")
        lines.append("")
        
        lines.append("### 🏷️ Variant Classes & Classifications")
        lines.append("")
        
        def make_dict_table(data_dict, col1, col2):
            table = [f"| {col1} | {col2} |", "| :--- | :---: |"]
            for k, v in sorted(data_dict.items(), key=lambda x: x[1], reverse=True):
                k_safe = str(k).replace('|', '\\|')
                table.append(f"| {k_safe} | {v:,} |")
            return "\n".join(table)
            
        if r['variant_classes']:
            lines.append("**Variant Classes**:")
            lines.append(make_dict_table(r['variant_classes'], "Class", "Variants"))
            lines.append("")
            
        if r['source_origins']:
            lines.append("**Source Origins (HiC vs. gnomAD)**:")
            lines.append(make_dict_table(r['source_origins'], "Source Origin", "Variants"))
            lines.append("")
            
        if r['impacts']:
            lines.append("**VEP Predicted Splicing/Translation Impacts**:")
            lines.append(make_dict_table(r['impacts'], "Impact Severity", "Variants"))
            lines.append("")

        if r['spip']:
            lines.append("### ✂️ SPiP Splicing Predictions (Top 5)")
            lines.append(make_dict_table(r['spip'], "Splicing Consequence / Score", "Variants"))
            lines.append("")
            
        if r['af']:
            lines.append("### 📈 Allele Frequency Characteristics")
            lines.append(f"Analysis performed using column **`{r['af']['column']}`**:")
            lines.append(f"- **Frequency Range**: {r['af']['min']:.3e} to {r['af']['max']:.3e}")
            lines.append(f"- **Mean Frequency**: {r['af']['mean']:.4e}")
            lines.append(f"- **Median Frequency**: {r['af']['median']:.4e}")
            lines.append(f"- **Missing/Absent Frequency**: {r['af']['missing_pct']:.1f}% of variants")
            lines.append("")
            
            # Add FAF details if present
            if "faf" in r['pop_stats']:
                faf = r['pop_stats']['faf']
                lines.append(f"**Filtering Allele Frequency (FAF95)** (using `{faf['column']}`):")
                lines.append(f"- **Max FAF**: {faf['max']:.4e}")
                lines.append(f"- **Variants with non-zero FAF**: {faf['non_zero_count']:,} / {faf['total']:,}")
                lines.append("")
                
            # Add NFE vs Non-NFE details if present
            if "nfe_comparison" in r['pop_stats']:
                nfe = r['pop_stats']['nfe_comparison']
                lines.append("**NFE (Non-Finnish European) vs. Non-NFE Allele Stats**:")
                lines.append("| Population Group | Total Allele Count (AC) | Total Allele Number (AN) | Calculated Overall AF |")
                lines.append("| :--- | :---: | :---: | :---: |")
                lines.append(f"| **NFE** | {nfe['nfe_ac']:,} | {nfe['nfe_an']:,} | {nfe['nfe_af']:.4e} |")
                lines.append(f"| **Non-NFE** | {nfe['non_nfe_ac']:,} | {nfe['non_nfe_an']:,} | {nfe['non_nfe_af']:.4e} |")
                lines.append("")
            
            # ASCII Histogram
            lines.append("#### 📊 Allele Frequency Distribution (ASCII)")
            lines.append("```text")
            lines.append(make_ascii_histogram(r['af']['raw_values']))
            lines.append("```")
            lines.append("")
            
            if "by_origin_mean" in r['af'] and r['af']["by_origin_mean"]:
                lines.append("**Mean Allele Frequency by Cohort Source**:")
                lines.append("| Origin Source | Mean Frequency |")
                lines.append("| :--- | :---: |")
                for origin, mean_val in r['af']["by_origin_mean"].items():
                    origin_safe = str(origin).replace('|', '\\|')
                    lines.append(f"| {origin_safe} | {mean_val:.4e} |")
                lines.append("")
                
        if r['predictors']:
            lines.append("### 🩺 Pathogenicity Scores Summary")
            lines.append("")
            lines.append("| Predictor Tool | Mean Score | Max Score | High-Risk Variants | Coverage |")
            lines.append("| :--- | :---: | :---: | :---: | :---: |")
            for pred_name, p in r['predictors'].items():
                threshold_info = ">0.5" if pred_name == "REVEL" else ">0.564 (Pathogenic)"
                lines.append(f"| {pred_name} | {p['mean']:.3f} | {p['max']:.3f} | {p['high_count']:,} ({threshold_info}) | {p['total']:,} variants |")
            lines.append("")
            
        lines.append("---")
        lines.append("")
        
    return "\n".join(lines)

def main():
    args = parse_args()
    
    files_to_analyze = []
    if os.path.isdir(args.input):
        pattern = os.path.join(args.input, "**", "*.parsed.clean.pq")
        files_to_analyze = glob.glob(pattern, recursive=True)
        if not files_to_analyze:
            pattern2 = os.path.join(args.input, "**", "*.pq")
            pattern3 = os.path.join(args.input, "**", "*.parquet")
            files_to_analyze = glob.glob(pattern2, recursive=True) + glob.glob(pattern3, recursive=True)
    elif os.path.isfile(args.input):
        files_to_analyze = [args.input]
        
    if not files_to_analyze:
        print(f"Error: No Parquet files found at input: {args.input}")
        return
        
    print(f"Found {len(files_to_analyze)} Parquet file(s) for EDA.")
    
    results = []
    processed_genes = set()
    for f in sorted(files_to_analyze):
        res = analyze_single_file(f, args.verbose)
        results.append(res)
        if "error" not in res:
            processed_genes.update(res['genes'])
            
    # Audit pipeline status if raw-dir is provided
    audit_info = None
    if args.raw_dir and os.path.exists(args.raw_dir):
        # Scan for raw genes
        raw_genes = set()
        for f in os.listdir(args.raw_dir):
            ext = os.path.splitext(f)[1].lower()
            if ext in ['.parquet', '.pq', '.xlsx', '.csv'] and f != 'test_write.parquet':
                base_name = os.path.splitext(os.path.basename(f))[0]
                gene = base_name.split('_')[0].strip().upper()
                raw_genes.add(gene)
                
        # Find failed/missing genes
        failed_genes = raw_genes - processed_genes
        failed_details = {}
        for g in failed_genes:
            failed_details[g] = get_gene_pipeline_status(g, args.log_dir)
            
        total = len(raw_genes) if raw_genes else 1
        success_count = len(processed_genes)
        failed_count = len(failed_genes)
        
        audit_info = {
            "total": total,
            "success_count": success_count,
            "success_pct": (success_count / total) * 100,
            "failed_count": failed_count,
            "failed_pct": (failed_count / total) * 100,
            "failed_genes": failed_genes,
            "failed_details": failed_details
        }
        
    report_md = format_markdown_report(results, audit_info)
    
    if args.output == "stdout":
        print("\n" + report_md)
    else:
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        with open(args.output, "w") as f:
            f.write(report_md)
        print(f"EDA Report successfully generated and saved to: {args.output}")

if __name__ == '__main__':
    main()
