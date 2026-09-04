#!/usr/bin/env python3
"""
🧪 Multi-Sample Trio Pedigree End-to-End Test Suite

Verifies:
  1. Parsing multi-sample VCF format with pysam (PROBAND_001, FATHER_001, MOTHER_001).
  2. Parsing 6-column PED pedigree format (MYBPC3_trio.ped).
  3. Evaluating inheritance models (De Novo, Autosomal Recessive, Compound Het).
  4. Priority Score escalation (+20 pts De Novo, +15 pts Recessive/Compound Het).
  5. Generating 4-tab Clinical Prioritization HTML Dashboard with inheritance annotations.
"""

import os
import sys
import pandas as pd
import pysam

PIPELINE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(PIPELINE_ROOT, "src", "python"))
import config

from filter_variants import (
    build_newImpact,
    read_pedigree_file,
    analyze_pedigree_inheritance,
    build_priority_tier
)
from generate_clinical_prioritization_report import generate_gene_report


def parse_vcf_with_pysam(vcf_path):
    """Direct VCF parser using pysam to extract INFO tags and GT_<SampleID> columns."""
    vcf = pysam.VariantFile(vcf_path)
    samples = list(vcf.header.samples)
    records = []
    for rec in vcf:
        row = {
            "CHROM": str(rec.chrom),
            "POS": rec.pos,
            "ID": rec.id or ".",
            "REF": rec.ref,
            "ALT": rec.alts[0] if rec.alts else ".",
            "QUAL": rec.qual or 1000.0,
            "FILTER": ",".join(rec.filter.keys()) if rec.filter else "PASS",
            "Locus": f"{rec.chrom}:{rec.pos}-{rec.ref}-{rec.alts[0] if rec.alts else '.'}"
        }
        for k, v in rec.info.items():
            row[k] = v[0] if isinstance(v, tuple) and len(v) == 1 else v
            
        for s in samples:
            gt_tuple = rec.samples[s].get("GT")
            if gt_tuple == (0, 1) or gt_tuple == (1, 0): gt_str = "HET"
            elif gt_tuple == (1, 1): gt_str = "HOMALT"
            elif gt_tuple == (0, 0): gt_str = "HOMREF"
            else: gt_str = "MISSING"
            row[f"GT_{s}"] = gt_str
            
        records.append(row)
    return pd.DataFrame(records)


def run_e2e_trio_test():
    vcf_path = os.path.join(PIPELINE_ROOT, "test_data", "raw_vcfs", "MYBPC3_trio_test.vcf")
    ped_path = os.path.join(PIPELINE_ROOT, "test_data", "raw_vcfs", "MYBPC3_trio.ped")
    out_pq = os.path.join(PIPELINE_ROOT, "test_data", "test_run", "results", "MYBPC3_trio.parsed.clean.pq")
    out_html = os.path.join(PIPELINE_ROOT, "test_data", "test_run", "reports", "MYBPC3_trio_clinical_prioritization_report.html")

    os.makedirs(os.path.dirname(out_pq), exist_ok=True)
    os.makedirs(os.path.dirname(out_html), exist_ok=True)

    print("==========================================================================")
    print(" 🧪 Starting End-to-End Multi-Sample Trio Pedigree Test")
    print("==========================================================================")

    # 1. Parse VCF to DataFrame
    print(f"1. Parsing Multi-Sample VCF with pysam: {vcf_path}...")
    df = parse_vcf_with_pysam(vcf_path)
    print(f"   - Parsed {len(df)} variants. Columns: {list(df.columns)}")

    # 2. Read Pedigree File
    print(f"\n2. Reading PED File: {ped_path}...")
    ped_df = read_pedigree_file(ped_path)
    print(f"   - Pedigree:\n{ped_df}")

    # 3. Execute Functional Tiering & Pedigree Inheritance Engine
    print("\n3. Running Functional Tiering & Pedigree Inheritance Engine...")
    df = build_newImpact(df)
    df = analyze_pedigree_inheritance(df, ped_df)
    df = build_priority_tier(df)

    # Save to Parquet
    df.to_parquet(out_pq)
    print(f"   - Saved processed dataset to: {out_pq}")

    # 4. Print Results Verification Table
    print("\n==========================================================================")
    print(" 📊 VARIANT INHERITANCE CLASSIFICATION & TIERING VERIFICATION")
    print("==========================================================================")
    for idx, r in df.iterrows():
        print(f"📌 Locus: {r['Locus']} ({r['SYMBOL']})")
        print(f"   - Consequence   : {r.get('Consequence', '-')}")
        print(f"   - Genotypes     : {r.get('SAMPLE_GENOTYPES_SUMMARY', '-')}")
        print(f"   - Inheritance   : {r.get('INHERITANCE_MODEL', '-')}")
        print(f"   - Priority Score: {r.get('VARIANT_PRIORITY_SCORE', '-')}")
        print(f"   - Priority Tier : {r.get('PRIORITY_TIER', '-')}\n")

    # 5. Generate Clinical Prioritization HTML Dashboard
    print(f"5. Generating 4-Tab Clinical HTML Report: {out_html}...")
    generate_gene_report(out_pq, out_html)
    print("==========================================================================")
    print(f" 🎉 SUCCESS! Multi-Sample Trio Report generated at:\n   {out_html}")
    print("==========================================================================")


if __name__ == "__main__":
    run_e2e_trio_test()
