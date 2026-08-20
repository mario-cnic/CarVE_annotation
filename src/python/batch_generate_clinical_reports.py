#!/usr/bin/env python3
"""
High-Throughput Parallel Clinical Prioritization Dashboard Generator
"""
import os
import sys
import glob
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

# Ensure repository root is on sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from src.python.generate_clinical_prioritization_report import generate_gene_report

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("batch_report_generator")

def process_single_gene(pq_file):
    run_dir = os.path.dirname(os.path.dirname(pq_file))
    rep_dir = os.path.join(run_dir, "reports")
    os.makedirs(rep_dir, exist_ok=True)
    gname = os.path.basename(pq_file).split(".")[0]
    out_html = os.path.join(rep_dir, f"{gname}_clinical_prioritization_report.html")
    try:
        generate_gene_report(pq_file, out_html)
        return (True, gname, out_html)
    except Exception as e:
        return (False, gname, str(e))

def main():
    runs = ["RUNS/run_20260813_1028", "RUNS/run_more_genes_20260817_1143", "RUNS/predictors_050826"]
    pq_files = []
    for run in runs:
        res_dir = os.path.join(run, "results")
        found = sorted(glob.glob(os.path.join(res_dir, "*.parsed.clean.pq")))
        pq_files.extend(found)
        logger.info(f"Found {len(found)} gene parquets in {run}")

    n_cpus = min(8, os.cpu_count() or 4)
    logger.info(f"Launching parallel generation for total {len(pq_files)} gene reports across {n_cpus} CPU workers...")
    
    success_count = 0
    fail_count = 0
    with ProcessPoolExecutor(max_workers=n_cpus) as executor:
        futures = {executor.submit(process_single_gene, pq): pq for pq in pq_files}
        for future in as_completed(futures):
            ok, gname, res = future.result()
            if ok:
                success_count += 1
                logger.info(f"[{success_count}/{len(pq_files)}] Completed {gname}")
            else:
                fail_count += 1
                logger.error(f"Failed {gname}: {res}")

    logger.info(f"BATCH COMPLETE: Successfully generated {success_count}/{len(pq_files)} clinical dashboards ({fail_count} failed).")

if __name__ == "__main__":
    main()
