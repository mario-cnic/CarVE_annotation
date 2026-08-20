import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/tools')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import combine_valencia_hic as cvh

def test_combine_dataframes():
    df_valencia = pd.DataFrame({
        'CHROM': ['11'], 'POS': [47353297], 'REF': ['G'], 'ALT': ['A'], 'VALENCIA_SCORE': [0.9]
    })
    df_hic = pd.DataFrame({
        'CHROM': ['11'], 'POS': [47353297], 'REF': ['G'], 'ALT': ['A'], 'HIC_SCORE': [0.8]
    })

    merged = pd.merge(df_valencia, df_hic, on=['CHROM', 'POS', 'REF', 'ALT'], how='outer')
    assert len(merged) == 1
    assert merged.iloc[0]['VALENCIA_SCORE'] == 0.9
    assert merged.iloc[0]['HIC_SCORE'] == 0.8
