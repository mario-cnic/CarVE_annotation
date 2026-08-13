import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import audit_run_results as arr

def test_audit_run(tmp_path):
    run_dir = tmp_path / "RUNS" / "test_run"
    raw_dir = run_dir / "input" / "by_gene"
    res_dir = run_dir / "results"
    
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(res_dir, exist_ok=True)

    # Create dummy raw input parquet
    raw_pq = raw_dir / "MYBPC3.pq"
    df_raw = pd.DataFrame({'CHROM': ['11'], 'POS': [47353297]})
    df_raw.to_parquet(raw_pq)

    # Create dummy final output clean parquet
    res_pq = res_dir / "MYBPC3.parsed.clean.pq"
    df_res = pd.DataFrame({
        'CHROM': ['11'],
        'POS': [47353297],
        'REF': ['G'],
        'ALT': ['A'],
        'Consequence': ['missense_variant'],
        'gnomAD_AF_joint': [0.0001],
        'REVEL_score': [0.85],
        'am_pathogenicity': [0.90],
        'SPiP': [0.10],
        'SPiP_interpretation': ['OK'],
        'splicevardb': ['OK'],
        'SpliceVault_top_events': ['None']
    })
    df_res.to_parquet(res_pq)

    df_results = arr.audit_run(str(run_dir), str(raw_dir))
    assert df_results is not None
    assert len(df_results) == 1
    assert df_results.iloc[0]['Gene'] == 'MYBPC3'
