import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/tools')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import validate_output as vo

def test_check_file_status(tmp_path):
    # Test valid CSV
    valid_csv = tmp_path / "valid.csv"
    df = pd.DataFrame({'a': [1, 2]})
    df.to_csv(valid_csv, index=False)
    assert vo.check_file_status(str(valid_csv)) == "VALID"

    # Test empty CSV with header
    empty_csv = tmp_path / "empty.csv"
    pd.DataFrame(columns=['a', 'b']).to_csv(empty_csv, index=False)
    assert vo.check_file_status(str(empty_csv)) == "EMPTY"

    # Test corrupted file
    corrupt_file = tmp_path / "corrupt.pq"
    corrupt_file.write_text("not a parquet file")
    assert vo.check_file_status(str(corrupt_file)) == "CORRUPTED"
