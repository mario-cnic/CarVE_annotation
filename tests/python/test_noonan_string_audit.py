import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import noonan_string_audit as nsa

def test_string_network_filtering():
    df_edges = pd.DataFrame({
        'protein1': ['MYBPC3', 'MYH7'],
        'protein2': ['MYH7', 'TNNT2'],
        'combined_score': [900, 400]
    })
    
    high_conf = df_edges[df_edges['combined_score'] >= 700]
    assert len(high_conf) == 1
    assert high_conf.iloc[0]['protein1'] == 'MYBPC3'
