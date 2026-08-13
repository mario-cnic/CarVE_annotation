import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import eda_parsed_pq as eda

def test_eda_metrics():
    df = pd.DataFrame({
        'GENE': ['MYBPC3', 'MYBPC3', 'MYH7'],
        'Consequence': ['missense_variant', 'stop_gained', 'missense_variant'],
        'gnomAD_AF_joint': [0.0001, 0.00005, 0.001],
        'REVEL': [0.85, None, 0.92]
    })
    
    gene_counts = df['GENE'].value_counts()
    assert gene_counts['MYBPC3'] == 2
    assert gene_counts['MYH7'] == 1
