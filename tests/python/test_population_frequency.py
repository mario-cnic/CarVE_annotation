"""Synthetic population-frequency availability and prioritization regressions."""
import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/python'))
from parsing.modules.population import resolve_population_af, GNOMAD_GRPMAX, GNOMAD_JOINT
from parsing.modules.acmg import build_acmg_criteria
from parsing.modules import scoring
from filter_and_summarize import filter_dataframe
from generate_clinical_prioritization_report import _extract_gnomadv4_af_grpmax


def test_resolution_per_row_and_raw_preservation():
    df = pd.DataFrame({GNOMAD_GRPMAX: [0, None, '.', 'bad', -1, 1.1, float('inf'), .002],
                       GNOMAD_JOINT: [.2, .1, None, None, None, None, None, .3],
                       'AF': [0] * 8, 'MAX_AF': [0] * 8})
    before = df.copy(deep=True)
    af, source, status = resolve_population_af(df)
    assert af.iloc[:2].tolist() == [0, .1]
    assert af.iloc[2:7].isna().all()
    assert af.iloc[7] == .002
    assert source.iloc[1] == GNOMAD_JOINT
    assert status.tolist() == ['available_frequency', 'available_frequency', 'missing_frequency'] + ['invalid_frequency'] * 4 + ['available_frequency']
    pd.testing.assert_frame_equal(df, before)


@pytest.mark.parametrize('value', [None, '', '.', '-', 'NA', 'nan'])
def test_missing_is_not_zero(value):
    af, _, status = resolve_population_af(pd.DataFrame({GNOMAD_GRPMAX: [value]}))
    assert af.isna().all()
    assert status.iloc[0] == 'missing_frequency'


def test_legacy_is_only_used_when_explicitly_requested():
    df = pd.DataFrame({'AF': [0], 'MAX_AF': [.1]})
    assert resolve_population_af(df)[0].isna().all()
    assert resolve_population_af(df, 'MAX_AF')[0].iloc[0] == .1
    with pytest.raises(ValueError):
        filter_dataframe(df, max_af=.001)


def test_acmg_missing_zero_common_and_joint_fallback():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0, None, '.'], GNOMAD_JOINT: [None, .2, .1, None], 'MAX_AF': [0] * 4})
    result = build_acmg_criteria(df)
    criteria = result.ACMG_CRITERIA.tolist()
    assert 'PM2' not in criteria[0] and 'BA1' not in criteria[0] and 'BS1' not in criteria[0]
    assert 'PM2_Supporting' in criteria[1]
    assert 'BA1' in criteria[2] and 'PM2' not in criteria[2]
    assert 'PM2' not in criteria[3]


def test_score_has_no_unknown_rarity_bonus(monkeypatch):
    monkeypatch.setattr(scoring, 'apply_disease_phenotype_curation', lambda data, phenotype: data)
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0], 'Consequence': ['missense_variant'] * 2,
                       'REVEL_score': [.8] * 2, 'NEW_IMPACT': ['MODERATE'] * 2,
                       'QC_STATUS': ['PASS'] * 2})
    result = scoring.build_priority_tier(df)
    assert result.VARIANT_PRIORITY_SCORE.iloc[1] - result.VARIANT_PRIORITY_SCORE.iloc[0] == 15
    assert result.POPULATION_AF_STATUS.iloc[0] == 'missing_frequency'
    assert 'PM2' not in result.ACMG_CRITERIA.iloc[0]


def test_filters_retain_unknown_without_zero_imputation():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0, .001, .1, 'bad']})
    result, stats = filter_dataframe(df, max_af=.001)
    assert result.index.tolist() == [0, 1, 2, 4]
    assert result.POPULATION_AF.isna().sum() == 2
    assert stats.AF_unknown_count.iloc[0] == 2
    strict = scoring.filter_freq(df, .001)
    assert strict.index.tolist() == [0, 1, 4]
    assert strict.POPULATION_AF.isna().sum() == 2


def test_report_preserves_zero_and_uses_joint_without_old_af():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0, None], GNOMAD_JOINT: [None, .3, .1], 'MAX_AF': [0] * 3})
    display, af, absent = _extract_gnomadv4_af_grpmax(df)
    assert display.iloc[0] == '-'
    assert display.iloc[1] == '0'
    assert 'joint AF; grpmax undefined' in display.iloc[2]
    assert pd.isna(af.iloc[0]) and not absent
    assert _extract_gnomadv4_af_grpmax(pd.DataFrame({'AF': [0]}))[2]


def test_unknown_cannot_promote_to_tier_one(monkeypatch):
    monkeypatch.setattr(scoring, 'apply_disease_phenotype_curation', lambda data, phenotype: data)
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0], 'Consequence': ['stop_gained'] * 2,
                       'pLI_gene_value': [1.] * 2, 'NEW_IMPACT': ['HIGH'] * 2,
                       'QC_STATUS': ['PASS'] * 2})
    result = scoring.build_priority_tier(df)
    assert not result.PRIORITY_TIER.iloc[0].startswith('Tier 1')
    assert result.PRIORITY_TIER.iloc[1].startswith('Tier 1')
    assert 'population AF unknown' in result.TIER_RATIONALE.iloc[0]


def test_all_unknown_filter_and_availability():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, '.'], 'AF': [0, 0]})
    assert len(scoring.filter_freq(df, .001)) == 2
    audit = scoring.audit_annotation_availability(df)['gnomAD Allele Frequency']
    assert audit['status'] == 'Unknown (no frequency evidence)'
    assert audit['coverage'].startswith('0/2')


def test_explicit_column_used_consistently_for_scoring_and_criteria(monkeypatch):
    monkeypatch.setattr(scoring, 'apply_disease_phenotype_curation', lambda data, phenotype: data)
    df = pd.DataFrame({GNOMAD_GRPMAX: [0], 'chosen_frequency': [.1], 'QC_STATUS': ['PASS']})
    result = scoring.build_priority_tier(df, freq_col='chosen_frequency')
    assert result.POPULATION_AF_SOURCE.iloc[0] == 'chosen_frequency'
    assert 'BA1' in result.ACMG_CRITERIA.iloc[0]
    assert 'PM2' not in result.ACMG_CRITERIA.iloc[0]
