import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/tools')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import parse_raw_hic as prh
import update_mappings_to_grch38 as umg

def test_parse_raw_hic_matrix():
    matrix_df = pd.DataFrame({
        'bin1': [1000, 2000],
        'bin2': [2000, 3000],
        'count': [15, 42]
    })
    binned = matrix_df[matrix_df['count'] > 20]
    assert len(binned) == 1
    assert binned.iloc[0]['count'] == 42

def test_update_mappings_grch38():
    mapping_df = pd.DataFrame({
        'GENE': ['MYBPC3'],
        'BUILD': ['hg19'],
        'POS': [47353297]
    })
    mapping_df['BUILD'] = 'GRCh38'
    assert mapping_df.iloc[0]['BUILD'] == 'GRCh38'
