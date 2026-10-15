"""Synthetic regressions for two-source SpliceAI export and score selection."""
import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/python'))
from parsing.modules import splicing, scoring
from parsing.modules.acmg import build_acmg_criteria


def test_two_sources_remain_separate_and_raw_columns_preserved():
    df = pd.DataFrame({'SpliceAI_pred_DS_AG': [.2, None, 0],
                       'spliceai_custom_DS_AG': [.8, .9, .7],
                       'spliceai_custom_MAX': [.8, .9, .7],
                       'SpliceAI': ['raw1', 'raw2', 'raw3']})
    before = df.copy(deep=True)
    result = splicing.build_spliceMAX(df, ['SpliceAI_pred_DS_AG', 'spliceai_custom_DS_AG'])
    pd.testing.assert_frame_equal(result[before.columns], before)
    assert result.spliceAI_MAX.iloc[0] == .2
    assert pd.isna(result.spliceAI_MAX.iloc[1])
    assert result.spliceAI_MAX.iloc[2] == 0
    assert result.SpliceAI_status.tolist() == ['scored', 'not_covered', 'scored']
    selected, methods = splicing._extract_spliceai_details(result)
    assert selected.tolist() == [.8, .9, .7]
    assert methods.str.startswith('Custom').all()


def test_selection_and_zero_priority_per_row():
    df = pd.DataFrame({'spliceai_custom_MAX': [0, None, .2, None],
                       'spliceAI_MAX': [.9, 0, .8, None]}, index=[10, 20, 30, 40])
    score, method = splicing._extract_spliceai_details(df)
    assert score.index.equals(df.index)
    assert score.iloc[:3].tolist() == [0, 0, .2]
    assert pd.isna(score.iloc[3])
    assert method.tolist() == ['Custom Window (-D 4999 / ~10kb)', 'VEP Standard (~50bp)',
                              'Custom Window (-D 4999 / ~10kb)', 'Not Covered']


@pytest.mark.parametrize('value', [None, '.', '', 'bad', -1, 1.1, float('inf')])
def test_invalid_or_missing_custom_falls_back(value):
    score, method = splicing._extract_spliceai_details(pd.DataFrame({'spliceai_custom_MAX': [value], 'spliceAI_MAX': [.4]}))
    assert score.iloc[0] == .4 and method.iloc[0].startswith('VEP')
    score, method = splicing._extract_spliceai_details(pd.DataFrame({'spliceai_custom_MAX': [value], 'spliceAI_MAX': [value]}))
    assert pd.isna(score.iloc[0]) and method.iloc[0] == 'Not Covered'


def test_component_maximum_missing_and_numeric_strings():
    df = pd.DataFrame({'SpliceAI_pred_DS_AG': ['.', '0', '.3', 'bad'],
                       'SpliceAI_pred_DS_AL': [None, None, '.6', float('inf')]})
    result = splicing.build_spliceMAX(df)
    assert pd.isna(result.spliceAI_MAX.iloc[0]) and pd.isna(result.spliceAI_MAX.iloc[3])
    assert result.spliceAI_MAX.iloc[1:3].tolist() == [0, .6]


def test_absent_sources_and_previously_imputed_uncovered_zero():
    for df in [pd.DataFrame(index=[0]), pd.DataFrame({'spliceAI_MAX': [0], 'SpliceAI_status': ['not_covered']})]:
        result = splicing.build_spliceMAX(df)
        score, method = splicing._extract_spliceai_details(result)
        assert pd.isna(score.iloc[0]) and method.iloc[0] == 'Not Covered'
        assert result.SpliceAI_status.iloc[0] == 'not_covered'


def test_unknown_does_not_supply_benign_splicing_evidence():
    df = pd.DataFrame({'spliceai_custom_MAX': [None, 0], 'Consequence': ['synonymous_variant'] * 2,
                       'am_pathogenicity': [.1] * 2, 'REVEL_score': [.1] * 2})
    result = build_acmg_criteria(df)
    assert 'BP4' not in result.ACMG_CRITERIA.iloc[0] and 'BP7' not in result.ACMG_CRITERIA.iloc[0]
    assert 'BP4' in result.ACMG_CRITERIA.iloc[1] and 'BP7' in result.ACMG_CRITERIA.iloc[1]


def test_missing_spliceai_does_not_erase_other_predictor_scores(monkeypatch):
    monkeypatch.setattr(scoring, 'apply_disease_phenotype_curation', lambda data, phenotype: data)
    df = pd.DataFrame({'spliceai_custom_MAX': [None, 0, None], 'spliceAI_MAX': [None] * 3,
                       'SPiP_prediction': [.8, .8, 0], 'Pangolin_max_score': [.3, .3, 0],
                       'QC_STATUS': ['PASS'] * 3})
    result = scoring.build_priority_tier(df)
    assert result.VARIANT_PRIORITY_SCORE.notna().all()
    assert result.VARIANT_PRIORITY_SCORE.iloc[0] == result.VARIANT_PRIORITY_SCORE.iloc[1] == 24
    assert result.VARIANT_PRIORITY_SCORE.iloc[2] == 0
    assert result.SPLICE_MAX_UNIFIED.isna().tolist() == [True, False, True]
    assert result.SPLICEAI_STATUS.tolist() == ['not_covered', 'scored', 'not_covered']
    assert result.SPLICEAI_METHOD.iloc[0] == 'Not Covered'


def test_explicit_error_preserved_without_poisoning_valid_fallback(monkeypatch):
    monkeypatch.setattr(scoring, 'apply_disease_phenotype_curation', lambda data, phenotype: data)
    df = pd.DataFrame({'SpliceAI_status': ['error'] * 2, 'spliceAI_MAX': [None] * 2,
                       'spliceai_custom_MAX': [None, 0], 'QC_STATUS': ['PASS'] * 2})
    result = scoring.build_priority_tier(splicing.build_spliceMAX(df))
    assert result.SpliceAI_status.tolist() == ['error', 'error']
    assert result.SPLICEAI_STATUS.tolist() == ['error', 'scored']
    assert pd.isna(result.SPLICE_MAX_UNIFIED.iloc[0]) and result.SPLICE_MAX_UNIFIED.iloc[1] == 0
