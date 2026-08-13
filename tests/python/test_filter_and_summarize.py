import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import filter_and_summarize as fas

def test_filter_dataframe():
    df = pd.DataFrame({
        'CHROM': ['11', '11', '14'],
        'POS': [47353297, 47353400, 23354000],
        'gnomAD_AF_joint': [0.0001, 0.005, 0.01],
        'REVEL': [0.85, 0.50, None],
        'AlphaMissense': [0.90, 0.30, 0.95],
        'SPiP_score': [0.10, 0.80, 0.05],
        'CADD_PHRED': [25.0, 15.0, 30.0],
        'Consequence': ['missense_variant', 'splice_donor_variant', 'stop_gained']
    })

    # Test AF filter
    res_af, stats_af = fas.filter_dataframe(df, max_af=0.001)
    assert len(res_af) == 1
    assert res_af.iloc[0]['POS'] == 47353297

    # Test REVEL filter
    res_revel, stats_revel = fas.filter_dataframe(df, min_revel=0.75)
    # Retains >= 0.75 or NA
    assert len(res_revel) == 2

    # Test Consequence filter
    res_cons, _ = fas.filter_dataframe(df, consequences="stop_gained,splice_donor_variant")
    assert len(res_cons) == 2

def test_load_data(tmp_path):
    tsv_file = tmp_path / "data.tsv"
    df_orig = pd.DataFrame({'GENE': ['MYBPC3'], 'POS': [100]})
    df_orig.to_csv(tsv_file, sep='\t', index=False)

    df_loaded = fas.load_data(str(tsv_file))
    assert len(df_loaded) == 1
    assert df_loaded.iloc[0]['GENE'] == 'MYBPC3'
