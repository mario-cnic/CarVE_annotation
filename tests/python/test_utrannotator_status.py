"""Synthetic checks for UTRAnnotator_status and its NEW_IMPACT effect (Whiffin et al. 2020 high-impact rule)."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/python'))
from parsing.modules.predictors import parse_utrannotator
from parsing.modules.acmg import build_newImpact

GAIN = '5_prime_UTR_premature_start_codon_gain_variant'
STOP_LOSS = '5_prime_UTR_uORF_stop_codon_loss_variant'
FRAMESHIFT = '5_prime_UTR_uORF_frameshift_variant'
START_LOSS = '5_prime_UTR_premature_start_codon_loss_variant'
STOP_GAIN = '5_prime_UTR_uORF_stop_codon_gain_variant'

CASES = [
    (GAIN, 'type=OutOfFrame_oORF:DistanceToCDS=20:KozakStrength=Moderate:KozakContext=ACCATGT', 'high_impact'),
    (GAIN, 'KozakStrength=Strong:type=inFrame_oORF:DistanceToCDS=9', 'high_impact'),
    (GAIN, 'type=uORF:KozakStrength=Strong:DistanceToStop=30', 'other'),
    (GAIN, 'type=OutOfFrame_oORF:KozakStrength=Weak', 'other'),
    (GAIN, 'type=OutOfFrame_oORF:KozakContext=-:KozakStrength=', 'unclassified'),
    (GAIN, 'type=uORF:KozakStrength=Weak&type=OutOfFrame_oORF:KozakStrength=Strong', 'high_impact'),
    (STOP_LOSS, 'AltStop=False:FrameWithCDS=outOfFrame:KozakStrength=Weak:Evidence=True', 'high_impact'),
    (STOP_LOSS, 'AltStop=False:FrameWithCDS=inFrame:KozakStrength=Moderate:Evidence=False', 'high_impact'),
    (STOP_LOSS, 'AltStop=True:FrameWithCDS=outOfFrame:KozakStrength=Strong:Evidence=True', 'other'),
    (STOP_LOSS, 'AltStop=False:FrameWithCDS=outOfFrame:KozakStrength=Weak:Evidence=False', 'other'),
    (STOP_LOSS, 'AltStop=False:FrameWithCDS=outOfFrame:KozakStrength=Weak:Evidence=NA', 'unclassified'),
    (FRAMESHIFT, 'ref_type=uORF:alt_type=OutOfFrame_oORF:KozakStrength=Strong:Evidence=False', 'high_impact'),
    (FRAMESHIFT, 'ref_type=uORF:alt_type=InFrame_oORF:KozakStrength=Weak:Evidence=True', 'high_impact'),
    (FRAMESHIFT, 'ref_type=uORF:alt_type=uORF:KozakStrength=Strong:Evidence=True', 'other'),
    (FRAMESHIFT, 'ref_type=OutOfFrame_oORF:alt_type=InFrame_oORF:KozakStrength=Strong', 'other'),
    (START_LOSS, 'type=OutOfFrame_oORF:KozakStrength=Strong:Evidence=True', 'other'),
    (STOP_GAIN, 'ref_type=OutOfFrame_oORF:KozakStrength=Strong:Evidence=True', 'other'),
    ('5_prime_UTR_unknown_variant', 'type=OutOfFrame_oORF:KozakStrength=Strong', 'unclassified'),
]


def test_uorf_classes_follow_the_high_impact_rule():
    df = pd.DataFrame({
        'Consequence': ['5_prime_UTR_variant'] * len(CASES),
        '5UTR_consequence': [c for c, _, _ in CASES],
        '5UTR_annotation': [a for _, a, _ in CASES],
        'Existing_uORFs': [1.0] * len(CASES),
    })
    assert parse_utrannotator(df).UTRAnnotator_status.tolist() == [s for _, _, s in CASES]


def test_rows_without_a_uorf_change():
    df = pd.DataFrame({
        'Consequence': ['5_prime_UTR_variant', '5_prime_UTR_variant&intron_variant', 'missense_variant'],
        '5UTR_consequence': [None, np.nan, None],
        '5UTR_annotation': [None, None, None],
        'Existing_uORFs': [2.0, np.nan, np.nan],
    }, index=[7, 3, 5])
    result = parse_utrannotator(df.copy())
    assert result.UTRAnnotator_status.tolist() == ['no_uorf_effect', 'not_evaluated', 'not_applicable']
    assert result.index.equals(df.index)
    pd.testing.assert_frame_equal(result[df.columns], df)


def test_missing_plugin_column_is_not_annotated():
    df = pd.DataFrame({'Consequence': ['5_prime_UTR_variant']})
    assert parse_utrannotator(df).UTRAnnotator_status.tolist() == ['not_annotated']


def test_new_impact_moderate_only_for_high_impact_and_never_downgrades():
    df = pd.DataFrame({
        'IMPACT': ['MODIFIER', 'MODIFIER', 'LOW', 'MODIFIER'],
        'Consequence': ['5_prime_UTR_variant'] * 4,
        '5UTR_consequence': [GAIN, GAIN, STOP_LOSS, GAIN],
        '5UTR_annotation': ['type=OutOfFrame_oORF:KozakStrength=Strong',
                            'type=uORF:KozakStrength=Strong',
                            'AltStop=False:KozakStrength=Moderate:Evidence=False',
                            'type=OutOfFrame_oORF:KozakStrength=Strong'],
        'Existing_uORFs': [1.0] * 4,
        'Pangolin_max_score': [np.nan, np.nan, np.nan, 0.9],
    })
    result = build_newImpact(parse_utrannotator(df))
    assert result.NEW_IMPACT.astype(str).tolist() == ['MODERATE', 'MODIFIER', 'MODERATE', 'HIGH']


def test_new_impact_derives_status_when_column_absent():
    df = pd.DataFrame({
        'IMPACT': ['MODIFIER'], 'Consequence': ['5_prime_UTR_variant'],
        '5UTR_consequence': [GAIN], '5UTR_annotation': ['type=inFrame_oORF:KozakStrength=Moderate'],
    })
    result = build_newImpact(df)
    assert result.NEW_IMPACT.astype(str).tolist() == ['MODERATE']
    assert 'UTRAnnotator_status' not in result.columns


def test_linked_entries_are_graded_under_their_own_consequence():
    # UTRAnnotatorPGP output: several consequences, each entry tagged with the one it belongs to
    gain = f'consequence={GAIN}:KozakStrength=Strong:type=OutOfFrame_oORF'
    stop_loss_weak = f'consequence={STOP_LOSS}:AltStop=False:Evidence=False:KozakStrength=Weak'
    start_loss = f'consequence={START_LOSS}:Evidence=True:KozakStrength=Strong:type=OutOfFrame_oORF'
    rows = [
        (f'{GAIN}&{STOP_LOSS}', f'{gain}&{stop_loss_weak}', 'high_impact'),
        (f'{START_LOSS}&{STOP_LOSS}', f'{start_loss}&{stop_loss_weak}', 'other'),
        (f'{GAIN},{START_LOSS}', f'{start_loss},{gain}', 'high_impact'),
        (f'{START_LOSS}&{STOP_LOSS}', start_loss, 'unclassified'),
    ]
    df = pd.DataFrame({
        'Consequence': ['5_prime_UTR_variant'] * len(rows),
        '5UTR_consequence': [c for c, _, _ in rows],
        '5UTR_annotation': [a for _, a, _ in rows],
        'Existing_uORFs': [1.0] * len(rows),
    })
    assert parse_utrannotator(df).UTRAnnotator_status.tolist() == [s for _, _, s in rows]


def test_unlinked_multi_consequence_row_keeps_previous_grading():
    # Output of the shared UTRAnnotator (no consequence= key): every entry graded under every term
    df = pd.DataFrame({
        'Consequence': ['5_prime_UTR_variant'],
        '5UTR_consequence': [f'{GAIN}&{STOP_LOSS}'],
        '5UTR_annotation': ['type=OutOfFrame_oORF:KozakStrength=Strong'],
        'Existing_uORFs': [1.0],
    })
    assert parse_utrannotator(df).UTRAnnotator_status.tolist() == ['high_impact']
