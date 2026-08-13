import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import generate_interactive_report as gir

def test_generate_report(tmp_path):
    pq_path = tmp_path / "sample.pq"
    df = pd.DataFrame({
        'CHROM': ['11', '14'],
        'POS': [47353297, 23354000],
        'REF': ['G', 'C'],
        'ALT': ['A', 'T'],
        'GENE': ['MYBPC3', 'MYH7'],
        'Consequence': ['missense_variant', 'stop_gained'],
        'gnomAD_AF_joint': [0.0001, 0.00005],
        'REVEL': [0.85, 0.95],
        'AlphaMissense': [0.90, 0.98],
        'SPiPv2_1_probability': [0.10, 0.05],
        'CADD_PHRED': [25.0, 35.0]
    })
    df.to_parquet(pq_path)

    out_html = tmp_path / "report.html"
    sys.argv = ['generate_interactive_report.py', '--input', str(pq_path), '--output', str(out_html)]
    gir.main()

    assert out_html.exists()
    content = out_html.read_text()
    assert "<html>" in content or "<!DOCTYPE html>" in content or "plotly" in content.lower()
