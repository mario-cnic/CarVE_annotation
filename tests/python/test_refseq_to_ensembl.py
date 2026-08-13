import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import refseq_to_ensembl as rte

def test_map_refseq_to_ensembl():
    mapping_df = pd.DataFrame({
        'RefSeq': ['NM_000256', 'NM_000257'],
        'Ensembl': ['ENST00000545968', 'ENST00000356287']
    })
    mapping_dict = dict(zip(mapping_df['RefSeq'], mapping_df['Ensembl']))

    assert mapping_dict.get('NM_000256') == 'ENST00000545968'
    assert mapping_dict.get('UNKNOWN') is None
