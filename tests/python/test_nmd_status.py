"""Synthetic checks for NMD_STATUS derived from the VEP NMD plugin flag."""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/python'))
from parsing.modules.predictors import parse_nmd


def test_status_per_consequence_and_cds_position():
    df = pd.DataFrame({
        'Consequence': ['stop_gained', 'frameshift_variant', 'stop_gained', 'splice_donor_variant',
                        'splice_acceptor_variant&frameshift_variant', 'stop_gained', 'missense_variant',
                        'frameshift_variant&NMD_transcript_variant', 'stop_gained'],
        'NMD': ['NMD_escaping_variant', None, '', None, 'NMD_escaping_variant', None, None, '-', None],
        'CDS_position': ['1200', '45-46', '300', None, '880', '?-12', '500', '77', '-'],
    }, index=range(10, 19))
    result = parse_nmd(df.copy())
    assert result.NMD_STATUS.tolist() == [
        'escape_predicted', 'nmd_predicted', 'nmd_predicted', 'not_evaluated',
        'escape_predicted', 'not_evaluated', 'not_applicable', 'nmd_predicted', 'not_evaluated']
    assert result.index.equals(df.index)
    pd.testing.assert_frame_equal(result[df.columns], df)


def test_missing_plugin_column_is_not_annotated_not_nmd():
    df = pd.DataFrame({'Consequence': ['stop_gained', 'splice_donor_variant', 'synonymous_variant'],
                       'CDS_position': ['10', None, '20']})
    assert parse_nmd(df).NMD_STATUS.tolist() == ['not_annotated', 'not_annotated', 'not_applicable']


def test_escape_flag_outside_lof_scope_is_ignored():
    df = pd.DataFrame({'Consequence': ['missense_variant'], 'NMD': ['NMD_escaping_variant'], 'CDS_position': ['10']})
    assert parse_nmd(df).NMD_STATUS.tolist() == ['not_applicable']
