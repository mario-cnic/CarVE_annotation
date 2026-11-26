"""Synthetic population-frequency availability and prioritization regressions."""
import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/python'))
from parsing.modules.population import resolve_population_af, ABSENT, GNOMAD_GRPMAX, GNOMAD_JOINT, GNOMAD_AN_JOINT
from parsing.modules import scoring



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


def test_table_status_absent_zero_joint_fallback_and_uncovered_record():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0, None, '.', None],
                       GNOMAD_JOINT: [None, .2, .1, None, None],
                       GNOMAD_AN_JOINT: [None, 1000, 1000, None, 0], 'MAX_AF': [.5] * 5})
    result = scoring.build_interpretation_inputs(df)
    assert result.POPULATION_AF_STATUS.tolist() == [ABSENT, 'available_frequency', 'available_frequency', ABSENT,
                                                    'missing_frequency']
    assert result.POPULATION_AF.iloc[1] == 0 and result.POPULATION_AF.iloc[2] == .1
    assert result.POPULATION_AF_SOURCE.iloc[2] == GNOMAD_JOINT


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


def test_no_variant_score_tier_or_acmg_columns():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0], 'Consequence': ['missense_variant'] * 2,
                       'REVEL_score': [.8] * 2})
    result = scoring.build_interpretation_inputs(df)
    assert result.POPULATION_AF_STATUS.tolist() == ['missing_frequency', 'available_frequency']
    assert not {'PRIORITY_TIER', 'VARIANT_PRIORITY_SCORE', 'ACMG_CRITERIA', 'TIER_RATIONALE',
                'SCORE_EXPLANATION', 'HPO_MATCH', 'DISEASE_GENE_EVIDENCE', 'ACMG_SF_V3_2',
                'NEW_IMPACT', 'QC_STATUS', 'QC_DETAIL'} & set(result.columns)


def test_filters_retain_unknown_without_zero_imputation():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0, .001, .1, 'bad']})
    strict = scoring.filter_freq(df, .001)
    assert strict.index.tolist() == [0, 1, 4]
    assert strict.POPULATION_AF.isna().sum() == 2


def test_filters_retain_absent():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, None, None], GNOMAD_JOINT: [None, .1, None],
                       GNOMAD_AN_JOINT: [None, 100, 0]})
    strict = scoring.filter_freq(df, .001)
    assert strict.index.tolist() == [0, 2]
    assert strict.POPULATION_AF.isna().all()
    assert strict.POPULATION_AF_STATUS.tolist() == [ABSENT, 'missing_frequency']


def test_absent_from_gnomad_status_and_frequency():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, 0, None], GNOMAD_JOINT: [None, 0, None],
                       GNOMAD_AN_JOINT: [None, 1000, 0], 'Consequence': ['stop_gained'] * 3,
                       'pLI_gene_value': [1.] * 3})
    result = scoring.build_interpretation_inputs(df)
    assert result.POPULATION_AF_STATUS.tolist() == [ABSENT, 'available_frequency', 'missing_frequency']
    assert pd.isna(result.POPULATION_AF.iloc[0]) and result.POPULATION_AF.iloc[1] == 0


def test_all_unknown_filter_and_availability():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, '.'], 'AF': [0, 0]})
    assert len(scoring.filter_freq(df, .001)) == 2
    audit = scoring.audit_annotation_availability(df)['gnomAD Allele Frequency']
    assert audit['status'] == 'Unknown (no frequency evidence)'
    assert audit['coverage'].startswith('0/2')
    audit = scoring.audit_annotation_availability(df.assign(**{GNOMAD_JOINT: None}))['gnomAD Allele Frequency']
    assert audit['status'] == 'Available' and audit['coverage'].startswith('2/2')


def test_explicit_column_used_for_the_table_frequency():
    df = pd.DataFrame({GNOMAD_GRPMAX: [0], 'chosen_frequency': [.1]})
    result = scoring.build_interpretation_inputs(df, freq_col='chosen_frequency')
    assert result.POPULATION_AF_SOURCE.iloc[0] == 'chosen_frequency' and result.POPULATION_AF.iloc[0] == .1


def test_min_absent_an_follows_the_pm2_threshold():
    from parsing.modules.population import min_absent_an
    assert min_absent_an(4e-5) == 74894
    assert min_absent_an(1e-4) == 29958


def test_site_an_takes_minimum_over_positions_and_missing_as_zero():
    from parsing.modules.population import site_allele_number
    df = pd.DataFrame({'gnomADv4AN_exomes': ['1400000', '948&1774&2374', None, ''],
                       'gnomADv4AN_genomes': ['150000', '80000&70000', '152000', None]})
    assert site_allele_number(df).tolist() == [1550000, 70948, 152000, 0]
    assert site_allele_number(df[['gnomADv4AN_genomes']]).tolist() == [150000, 70000, 152000, 0]
    assert site_allele_number(pd.DataFrame({'x': [1]})).isna().all()


def test_absence_needs_enough_site_an():
    from parsing.modules.population import ABSENT_LOW_AN, SITE_AN
    # absent at well-genotyped site, absent at poorly genotyped site, absent where gnomAD has no
    # callable genotype, and a present allele whose site AN is irrelevant
    df = pd.DataFrame({GNOMAD_GRPMAX: [None, None, None, 0], GNOMAD_JOINT: [None, None, None, 0],
                       GNOMAD_AN_JOINT: [None, None, None, 1000],
                       'gnomADv4AN_exomes': ['1400000', '', '', '10'], 'gnomADv4AN_genomes': ['150000', '74893', '', '10'],
                       'Consequence': ['stop_gained'] * 4, 'pLI_gene_value': [1.] * 4})
    af, _, status = resolve_population_af(df)
    assert status.tolist() == [ABSENT, ABSENT_LOW_AN, ABSENT_LOW_AN, 'available_frequency']
    assert af.iloc[:3].isna().all()
    result = scoring.build_interpretation_inputs(df.copy())
    assert result[SITE_AN].tolist() == [1550000, 74893, 0, 20]
    assert result.POPULATION_AF_STATUS.tolist() == [ABSENT, ABSENT_LOW_AN, ABSENT_LOW_AN, 'available_frequency']


def test_absence_without_site_an_annotation_keeps_old_status():
    df = pd.DataFrame({GNOMAD_GRPMAX: [None], GNOMAD_JOINT: [None], GNOMAD_AN_JOINT: [None]})
    assert resolve_population_af(df)[2].tolist() == [ABSENT]


def test_pm2_threshold_option_changes_min_an(monkeypatch):
    from parsing.modules import population
    monkeypatch.setattr(population, 'PM2_ABSENCE_AF', 1e-4)
    df = pd.DataFrame({GNOMAD_GRPMAX: [None], GNOMAD_JOINT: [None], GNOMAD_AN_JOINT: [None],
                       'gnomADv4AN_genomes': ['40000']})
    assert resolve_population_af(df)[2].tolist() == [population.ABSENT]
