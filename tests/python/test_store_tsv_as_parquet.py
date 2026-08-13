import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import store_tsv_as_parquet as stap

def test_store_tsv_as_parquet(tmp_path):
    tsv_file = tmp_path / "input.tsv"
    df_orig = pd.DataFrame({'CHROM': ['11'], 'POS': [100], 'GENE': ['MYBPC3']})
    df_orig.to_csv(tsv_file, sep='\t', index=False)

    stap.batch_convert_tsv_to_parquet(folder_path=str(tmp_path), output_folder="compressed")

    output_parquet = tmp_path / "compressed" / "input.parquet"
    assert output_parquet.exists()
    df_read = pd.read_parquet(output_parquet)
    assert len(df_read) == 1
    assert df_read.iloc[0]['GENE'] == 'MYBPC3'
