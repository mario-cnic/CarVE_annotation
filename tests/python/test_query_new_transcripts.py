import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import query_new_transcripts as qnt

def test_extract_refseq_from_raw(tmp_path):
    raw_file = tmp_path / "ACADVL_data.csv"
    df = pd.DataFrame({'cDNA': ['NM_000018.3:c.100A>G', 'c.200C>T']})
    df.to_csv(raw_file, index=False)

    refseq = qnt.extract_refseq_from_raw(str(raw_file))
    assert refseq == "NM_000018.3"

def test_resolve_gene_with_df_local():
    df_local = pd.DataFrame({
        'Gene name': ['MYBPC3'],
        'Ensembl Canonical': [1.0],
        'Transcript stable ID': ['ENST00000545968'],
        'RefSeq match transcript (MANE Select)': ['NM_000256']
    })

    refseq, enst = qnt.resolve_gene('MYBPC3', raw_file_path=None, df_local=df_local)
    assert refseq == 'NM_000256'
    assert enst == 'ENST00000545968'
