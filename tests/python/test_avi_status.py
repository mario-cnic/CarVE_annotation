"""Synthetic checks for AVI_PHRED/AVI_RAW/AVI_status derived from the VEP AVI plugin columns."""
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'src/python'))
from parsing.modules.predictors import parse_avi


def test_status_per_variant_class_and_value():
    df = pd.DataFrame({
        'VARIANT_CLASS': ['SNV', 'SNV', 'SNV', 'deletion', 'insertion', 'SNV', 'SNV', None],
        'AVI_PHRED': ['25.1', None, '.', None, None, '0', 'x&y', None],
        'AVI_RAW': ['1.5', None, '.', None, None, '-1.2', '0.3', None],
    })
    result = parse_avi(df.copy())
    assert result.AVI_status.tolist() == [
        'scored', 'not_covered', 'not_covered', 'not_applicable', 'not_applicable', 'scored', 'error', 'not_applicable']
    np.testing.assert_allclose(result.AVI_PHRED.tolist(), [25.1, np.nan, np.nan, np.nan, np.nan, 0.0, np.nan, np.nan])
    np.testing.assert_allclose(result.AVI_RAW.tolist(), [1.5, np.nan, np.nan, np.nan, np.nan, -1.2, np.nan, np.nan])


def test_numeric_input_and_missing_raw_column():
    df = pd.DataFrame({'VARIANT_CLASS': ['SNV', 'SNV'], 'AVI_PHRED': [59.63756, np.nan]})
    result = parse_avi(df)
    assert result.AVI_status.tolist() == ['scored', 'not_covered']
    assert result.AVI_PHRED.tolist()[0] == 59.63756
    assert result.AVI_RAW.isna().all()


def test_absent_columns_are_not_annotated():
    for df in (pd.DataFrame({'VARIANT_CLASS': ['SNV', 'deletion']}), pd.DataFrame({'AVI_PHRED': ['20', '30']})):
        result = parse_avi(df)
        assert result.AVI_status.tolist() == ['not_annotated', 'not_annotated']
        assert result.AVI_PHRED.isna().all() and result.AVI_RAW.isna().all()


def test_vep_plugin_param_and_whitelist():
    vep = (REPO / 'modules/local/vep.nf').read_text()
    assert re.search(r'--plugin AVI,file=\$\{params\.vep_avi\},raw=1\b', vep)
    cols = (REPO / 'resources/all_but_old_gnomad_vep_cols.txt').read_text().split()
    assert 'AVI_PHRED' in cols and 'AVI_RAW' in cols
    assert re.search(r'^\s*vep_avi\s*=', (REPO / 'nextflow.config').read_text(), re.M)
    assert 'vep_avi' in (REPO / 'run_annotate_vcf.sh').read_text()
