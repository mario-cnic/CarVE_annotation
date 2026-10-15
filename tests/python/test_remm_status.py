"""Synthetic checks for ReMM_score/ReMM_status derived from the VEP --custom ReMM column."""
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / 'src/python'))
from parsing.modules.predictors import parse_remm


def test_status_per_variant_class():
    df = pd.DataFrame({
        'VARIANT_CLASS': ['SNV', 'SNV', 'SNV', 'deletion', 'insertion', 'substitution', 'SNV', 'SNV'],
        'ReMM': ['0.998689', None, '.', None, '0.5', '0.4', '0', '0.1&0.2'],
    }, index=range(5, 13))
    result = parse_remm(df.copy())
    assert result.ReMM_status.tolist() == [
        'scored', 'not_covered', 'not_covered', 'not_applicable', 'not_applicable', 'not_applicable', 'scored', 'error']
    np.testing.assert_allclose(result.ReMM_score.tolist(), [0.998689, np.nan, np.nan, np.nan, np.nan, np.nan, 0.0, np.nan])
    assert result.index.equals(df.index)
    assert result.ReMM.isna().tolist() == [False, True, False, True, True, True, False, False]
    pd.testing.assert_series_equal(result.VARIANT_CLASS, df.VARIANT_CLASS)


def test_numeric_column_and_multiallelic_alts_share_the_position_score():
    df = pd.DataFrame({'VARIANT_CLASS': ['SNV', 'SNV'], 'ReMM': [0.27834, 0.27834]})
    result = parse_remm(df)
    assert result.ReMM_status.tolist() == ['scored', 'scored']
    assert result.ReMM_score.tolist() == [0.27834, 0.27834]


def test_missing_column_is_not_annotated_not_zero():
    result = parse_remm(pd.DataFrame({'VARIANT_CLASS': ['SNV', 'deletion']}))
    assert result.ReMM_status.tolist() == ['not_annotated', 'not_annotated']
    assert result.ReMM_score.isna().all()


def test_missing_variant_class_is_not_annotated_not_scored():
    result = parse_remm(pd.DataFrame({'ReMM': ['0.9', '0.3']}))
    assert result.ReMM_status.tolist() == ['not_annotated', 'not_annotated']
    assert result.ReMM_score.isna().all()


def test_vep_custom_is_exact_bigwig_and_whitelisted():
    vep = (REPO / 'modules/local/vep.nf').read_text()
    assert re.search(r'--custom file=\$\{params\.vep_remm\},short_name=ReMM,format=bigwig,type=exact\b', vep)
    assert 'ReMM' in (REPO / 'resources/all_but_old_gnomad_vep_cols.txt').read_text().split()
    assert re.search(r'^\s*vep_remm\s*=', (REPO / 'nextflow.config').read_text(), re.M)
    assert 'vep_remm' in (REPO / 'run_annotate_vcf.sh').read_text()
