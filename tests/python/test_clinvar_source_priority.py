"""Synthetic ClinVar source precedence and same-source review regressions."""
import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/python'))
from parsing.modules.clinvar import resolve_clinvar
from parsing.modules.acmg import _extract_clinvar_status, build_newImpact
from parsing.modules import scoring


def test_custom_authoritative_and_fallback_per_row():
    df = pd.DataFrame({
        'ClinVarVCF_CLNSIG': ['Pathogenic', 'Benign', 'Uncertain_significance', 'Conflicting_classifications_of_pathogenicity', '.', None],
        'clinvar_clnsig': ['Benign', 'Pathogenic', 'Pathogenic', 'Pathogenic', 'Likely_pathogenic', ''],
        'CLIN_SIG': ['benign'] * 6,
        'ClinVarVCF_CLNREVSTAT': ['reviewed_by_expert_panel'] * 6,
        'clinvar_review': ['criteria_provided,_single_submitter'] * 6,
    })
    path, benign = _extract_clinvar_status(df)
    assert path.tolist() == [True, False, False, False, True, False]
    assert benign.tolist() == [False, True, False, False, False, True]
    assert df.CLINVAR_SOURCE.tolist() == ['ClinVarVCF_CLNSIG'] * 4 + ['clinvar_clnsig', 'CLIN_SIG']
    assert df.ClinVar_Stars.iloc[:5].tolist() == [3, 3, 3, 3, 1]
    assert pd.isna(df.ClinVar_Stars.iloc[5])


def test_missing_custom_review_does_not_borrow_old_stars():
    df = pd.DataFrame({'ClinVarVCF_CLNSIG': ['Pathogenic'], 'clinvar_clnsig': ['Pathogenic'],
                       'clinvar_review': ['reviewed_by_expert_panel'], 'ClinVar_Stars': [4]})
    _extract_clinvar_status(df)
    assert pd.isna(df.ClinVar_Stars.iloc[0])
    assert df.CLINVAR_SELECTED_REVIEW.iloc[0] == ''


@pytest.mark.parametrize('review,expected', [
    ('practice_guideline', 4), ('reviewed_by_expert_panel', 3),
    ('criteria_provided,_multiple_submitters,_no_conflicts', 2),
    ('criteria_provided&_multiple_submitters&_no_conflicts', 2),
    ('criteria_provided,_single_submitter', 1),
    ('criteria_provided,_conflicting_classifications', 1),
    ('no_assertion_criteria_provided', 0), ('unrecognized', None), ('', None),
])
def test_review_levels(review, expected):
    df = pd.DataFrame({'ClinVarVCF_CLNSIG': ['Pathogenic'], 'ClinVarVCF_CLNREVSTAT': [review]})
    stars = resolve_clinvar(df)[3].iloc[0]
    assert pd.isna(stars) if expected is None else stars == expected


@pytest.mark.parametrize('classification,path,benign', [
    ('Pathogenic/Likely_pathogenic', True, False), ('Benign&Likely_benign', False, True),
    ('Pathogenic&Benign', False, False), ('Pathogenic/Uncertain_significance', False, False),
    ('Conflicting_interpretations_of_pathogenicity', False, False), ('risk_factor', False, False),
])
def test_ambiguous_or_other_classification_not_promoted(classification,path,benign):
    p,b = _extract_clinvar_status(pd.DataFrame({'ClinVarVCF_CLNSIG': [classification]}))
    assert bool(p.iloc[0]) == path and bool(b.iloc[0]) == benign


def test_scoring_and_impact_use_custom_source(monkeypatch):
    monkeypatch.setattr(scoring, 'apply_disease_phenotype_curation', lambda data, phenotype: data)
    df = pd.DataFrame({'ClinVarVCF_CLNSIG': ['Pathogenic', 'Pathogenic'],
                       'ClinVarVCF_CLNREVSTAT': ['reviewed_by_expert_panel', None],
                       'clinvar_clnsig': ['Benign'] * 2, 'ClinVar_Stars': [4] * 2,
                       'QC_STATUS': ['PASS'] * 2, 'Consequence': ['missense_variant'] * 2,
                       'IMPACT': ['MODERATE'] * 2})
    result = scoring.build_priority_tier(build_newImpact(df))
    assert result.NEW_IMPACT.tolist() == ['HIGH', 'HIGH']
    assert result.VARIANT_PRIORITY_SCORE.iloc[0] - result.VARIANT_PRIORITY_SCORE.iloc[1] == 5
    assert result.CLINVAR_SOURCE.eq('ClinVarVCF_CLNSIG').all()
    assert 'PS1' in result.ACMG_CRITERIA.iloc[0]
    assert pd.isna(result.ClinVar_Stars.iloc[1])


def test_legacy_alias_and_missing_sources():
    df = pd.DataFrame({'CLNSIG': ['Likely_pathogenic', '.'], 'ClinVar_Stars': [3, 3]})
    p,b = _extract_clinvar_status(df)
    assert p.tolist() == [True, False]
    assert df.CLINVAR_SOURCE.tolist() == ['CLNSIG', '']
    assert df.ClinVar_Stars.iloc[0] == 3
    assert pd.isna(df.ClinVar_Stars.iloc[1])
    assert not _extract_clinvar_status(pd.DataFrame(index=[0]))[0].iloc[0]
