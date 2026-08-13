import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

def test_filter_gene_transcript_missing(tmp_path):
    mapping_file = tmp_path / "mapping.txt"
    df = pd.DataFrame({
        'Gen': ['MYBPC3', 'UNKNOWN_GENE'],
        'RefSeq': ['NM_000256', 'N/A'],
        'Ensembl': ['ENST00000545968', 'N/A']
    })
    df.to_csv(mapping_file, index=False)

    df_read = pd.read_csv(mapping_file, keep_default_na=False)
    valid_mask = (df_read['RefSeq'] != 'N/A') & (df_read['Ensembl'] != 'N/A') & (df_read['RefSeq'] != '')
    filtered_df = df_read[valid_mask]

    assert len(filtered_df) == 1
    assert filtered_df.iloc[0]['Gen'] == 'MYBPC3'
