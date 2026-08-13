import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import split_input_by_gene as sibg

def test_split_input_by_gene(tmp_path):
    multi_gene_file = tmp_path / "multi_gene.pq"
    df = pd.DataFrame({
        'GENE': ['MYBPC3', 'MYBPC3', 'MYH7'],
        'CHROM': ['11', '11', '14'],
        'POS': [47353297, 47353400, 23354000],
        'REF': ['G', 'G', 'C'],
        'ALT': ['A', 'A', 'T']
    })
    df.to_parquet(multi_gene_file)

    out_dir = tmp_path / "split_output"

    # Simulate command line arguments
    sys.argv = ['split_input_by_gene.py', '--input', str(multi_gene_file), '--output-dir', str(out_dir)]
    sibg.main()

    assert (out_dir / "MYBPC3.pq").exists()
    assert (out_dir / "MYH7.pq").exists()

    df_mybpc3 = pd.read_parquet(out_dir / "MYBPC3.pq")
    assert len(df_mybpc3) == 2
