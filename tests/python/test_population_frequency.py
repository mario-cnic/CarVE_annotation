"""Synthetic population-frequency availability and prioritization regressions."""
import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/python'))
from parsing.modules.population import resolve_population_af, ABSENT, GNOMAD_GRPMAX, GNOMAD_JOINT, GNOMAD_AN_JOINT
from parsing.modules.acmg import build_acmg_criteria
from parsing.modules import scoring

DISEASE_GENE_LIST = str(Path(__file__).resolve().parents[1] / 'fixtures/mini_disease_gene_list.csv')


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
    assert status.tolist() == ['available_frequency', 'available_frequency', ABSENT] + ['invalid_frequency'] * 4 + ['available_frequency']
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
        scoring.filter_freq(df, .001)


def test_acmg_absent_zero_common_joint_fallback_and_uncovered_record():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0, None, '.', None],
                       GNOMAD_JOINT: [None, .2, .1, None, None],
                       GNOMAD_AN_JOINT: [None, 1000, 1000, None, 0], 'MAX_AF': [.5] * 5})
    result = build_acmg_criteria(df)
    criteria = result.ACMG_CRITERIA.tolist()
    assert 'PM2_Supporting' in criteria[0] and 'BA1' not in criteria[0] and 'BS1' not in criteria[0]
    assert 'PM2_Supporting' in criteria[1]
    assert 'BA1' in criteria[2] and 'PM2' not in criteria[2]
    assert 'PM2_Supporting' in criteria[3]
    assert 'PM2' not in criteria[4] and 'BA1' not in criteria[4]


def test_absent_requires_gnomad_annotation_and_no_record():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, None, None, 'bad', None],
                       GNOMAD_JOINT: [None, None, .1, None, 'bad'],
                       GNOMAD_AN_JOINT: [None, 0, 100, None, None]})
    af, source, status = resolve_population_af(df)
    assert status.tolist() == [ABSENT, 'missing_frequency', 'available_frequency', 'invalid_frequency', 'invalid_frequency']
    assert pd.isna(af.iloc[0]) and source.iloc[0] == ''
    assert resolve_population_af(df[[GNOMAD_GRPMAX]])[2].iloc[0] == 'missing_frequency'
    assert resolve_population_af(df.drop(columns=GNOMAD_AN_JOINT))[2].iloc[1] == ABSENT
    assert resolve_population_af(df, GNOMAD_JOINT)[2].iloc[0] == ABSENT
    assert resolve_population_af(df.assign(MAX_AF=None), 'MAX_AF')[2].iloc[0] == 'missing_frequency'


def test_score_has_no_unknown_rarity_bonus(monkeypatch):
    monkeypatch.setattr(scoring, 'apply_disease_phenotype_curation', lambda data, gene_list, phenotype: data)
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0], 'Consequence': ['missense_variant'] * 2,
                       'REVEL_score': [.8] * 2, 'NEW_IMPACT': ['MODERATE'] * 2,
                       'QC_STATUS': ['PASS'] * 2})
    result = scoring.build_priority_tier(df, DISEASE_GENE_LIST)
    assert result.VARIANT_PRIORITY_SCORE.iloc[1] - result.VARIANT_PRIORITY_SCORE.iloc[0] == 15
    assert result.POPULATION_AF_STATUS.iloc[0] == 'missing_frequency'
    assert 'PM2' not in result.ACMG_CRITERIA.iloc[0]


def test_filters_retain_unknown_without_zero_imputation():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0, .001, .1, 'bad']})
    strict = scoring.filter_freq(df, .001)
    assert strict.index.tolist() == [0, 1, 4]
    assert strict.POPULATION_AF.isna().sum() == 2


def test_filters_retain_absent_and_dominant_ranks_it_candidate():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, None, None], GNOMAD_JOINT: [None, .1, None],
                       GNOMAD_AN_JOINT: [None, 100, 0],
                       'gene_priority': ['prioritary'] * 3, 'NEW_IMPACT': ['HIGH'] * 3})
    strict = scoring.filter_freq(df, .001)
    assert strict.index.tolist() == [0, 2]
    assert strict.POPULATION_AF.isna().all()
    assert strict.POPULATION_AF_STATUS.tolist() == [ABSENT, 'missing_frequency']
    ranks = scoring.classification_dominant(df, GNOMAD_GRPMAX, .001).Variant_rank.tolist()
    assert ranks == ['Candidate', 'Excluded', 'Unknown frequency']


def test_unknown_cannot_promote_to_tier_one(monkeypatch):
    monkeypatch.setattr(scoring, 'apply_disease_phenotype_curation', lambda data, gene_list, phenotype: data)
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0], 'Consequence': ['stop_gained'] * 2,
                       'pLI_gene_value': [1.] * 2, 'NEW_IMPACT': ['HIGH'] * 2,
                       'QC_STATUS': ['PASS'] * 2})
    result = scoring.build_priority_tier(df, DISEASE_GENE_LIST)
    assert not result.PRIORITY_TIER.iloc[0].startswith('Tier 1')
    assert result.PRIORITY_TIER.iloc[1].startswith('Tier 1')
    assert 'population AF unknown' in result.TIER_RATIONALE.iloc[0]


def test_absent_from_gnomad_earns_pm2_and_can_reach_tier_one(monkeypatch):
    monkeypatch.setattr(scoring, 'apply_disease_phenotype_curation', lambda data, gene_list, phenotype: data)
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0, None], GNOMAD_JOINT: [None, 0, None],
                       GNOMAD_AN_JOINT: [None, 1000, 0], 'Consequence': ['stop_gained'] * 3,
                       'pLI_gene_value': [1.] * 3, 'NEW_IMPACT': ['HIGH'] * 3,
                       'QC_STATUS': ['PASS'] * 3})
    result = scoring.build_priority_tier(df, DISEASE_GENE_LIST)
    assert result.POPULATION_AF_STATUS.tolist() == [ABSENT, 'available_frequency', 'missing_frequency']
    assert pd.isna(result.POPULATION_AF.iloc[0]) and result.POPULATION_AF.iloc[1] == 0
    assert result.VARIANT_PRIORITY_SCORE.iloc[0] == result.VARIANT_PRIORITY_SCORE.iloc[1]
    assert result.VARIANT_PRIORITY_SCORE.iloc[0] - result.VARIANT_PRIORITY_SCORE.iloc[2] == 15
    assert result.PRIORITY_TIER.iloc[0].startswith('Tier 1')
    assert not result.PRIORITY_TIER.iloc[2].startswith('Tier 1')
    assert 'PM2_Supporting' in result.ACMG_CRITERIA.iloc[0]
    assert 'absent from gnomAD v4' in result.SCORE_EXPLANATION.iloc[0]


def test_all_unknown_filter_and_availability():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, '.'], 'AF': [0, 0]})
    assert len(scoring.filter_freq(df, .001)) == 2
    audit = scoring.audit_annotation_availability(df)['gnomAD Allele Frequency']
    assert audit['status'] == 'Unknown (no frequency evidence)'
    assert audit['coverage'].startswith('0/2')
    audit = scoring.audit_annotation_availability(df.assign(**{GNOMAD_JOINT: None}))['gnomAD Allele Frequency']
    assert audit['status'] == 'Available' and audit['coverage'].startswith('2/2')


def test_explicit_column_used_consistently_for_scoring_and_criteria(monkeypatch):
    monkeypatch.setattr(scoring, 'apply_disease_phenotype_curation', lambda data, gene_list, phenotype: data)
    df = pd.DataFrame({GNOMAD_GRPMAX: [0], 'chosen_frequency': [.1], 'QC_STATUS': ['PASS']})
    result = scoring.build_priority_tier(df, DISEASE_GENE_LIST, freq_col='chosen_frequency')
    assert result.POPULATION_AF_SOURCE.iloc[0] == 'chosen_frequency'
    assert 'BA1' in result.ACMG_CRITERIA.iloc[0]
    assert 'PM2' not in result.ACMG_CRITERIA.iloc[0]
