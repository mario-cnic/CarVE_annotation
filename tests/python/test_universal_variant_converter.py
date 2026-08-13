import sys
import os
import pytest
import pandas as pd

# Add src/python to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import universal_variant_converter as uvc

def test_clean_string_for_vcf():
    assert uvc.clean_string_for_vcf(None) == "."
    assert uvc.clean_string_for_vcf(float('nan')) == "."
    assert uvc.clean_string_for_vcf("  hello;world=123  ") == "hello,world:123"
    assert uvc.clean_string_for_vcf("María_Pérez") == "Maria_Perez"

def test_normalize_chrom():
    assert uvc.normalize_chrom(None) is None
    assert uvc.normalize_chrom("chr11") == "11"
    assert uvc.normalize_chrom("11") == "11"
    assert uvc.normalize_chrom("chrX") == "X"
    assert uvc.normalize_chrom("chrM") == "MT"
    assert uvc.normalize_chrom("MT") == "MT"

def test_compute_signed_intron_offset():
    # Exonic
    assert uvc.compute_signed_intron_offset("c.526C>T") == (0, "exonic")
    assert uvc.compute_signed_intron_offset("") == (0, "exonic")
    assert uvc.compute_signed_intron_offset(None) == (0, "exonic")

    # Donor (+)
    assert uvc.compute_signed_intron_offset("c.123+5G>A") == (5, "donor")
    assert uvc.compute_signed_intron_offset("NM_000257.3:c.100+12A>G") == (12, "donor")

    # Acceptor (-)
    assert uvc.compute_signed_intron_offset("c.456-12C>T") == (-12, "acceptor")
    assert uvc.compute_signed_intron_offset("c.10-2A>G") == (-2, "acceptor")

def test_convert_coordinates():
    df = pd.DataFrame({
        'CHROM': ['chr11', '14'],
        'POS': [47353297, 23354000],
        'REF': ['G', 'C'],
        'ALT': ['A', 'T'],
        'METADATA': ['sample1', 'sample2']
    })
    id_cols = ['CHROM', 'POS', 'REF', 'ALT']
    variants, errors = uvc.convert_coordinates(df, id_cols, input_build="GRCh38")
    assert len(variants) == 2
    assert len(errors) == 0
    assert variants[0][0] == '11'
    assert variants[0][1] == 47353297
    assert variants[0][3] == 'G'
    assert variants[0][4] == 'A'
    assert 'METADATA=sample1' in variants[0][7]

def test_load_input_dataframe(tmp_path):
    csv_path = tmp_path / "test.csv"
    df_orig = pd.DataFrame({'CHROM': ['11'], 'POS': [100], 'REF': ['A'], 'ALT': ['G']})
    df_orig.to_csv(csv_path, index=False)

    df_loaded = uvc.load_input_dataframe(str(csv_path))
    assert len(df_loaded) == 1
    assert df_loaded.iloc[0]['CHROM'] == 11

def test_write_vcf(tmp_path):
    out_vcf = tmp_path / "out.vcf"
    variants = [("11", 47353297, ".", "G", "A", ".", "PASS", "GENE=MYBPC3")]
    errors = []
    uvc.write_vcf(str(out_vcf), variants, errors, ["GENE"])
    assert out_vcf.exists()
    content = out_vcf.read_text()
    assert "##assembly=GRCh38" in content
    assert "11\t47353297" in content
