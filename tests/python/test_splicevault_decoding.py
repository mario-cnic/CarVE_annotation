"""SpliceVault cryptic-site positions are relative to the annotated splice site, along the transcript (PRED-7)."""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python/parsing")))
from modules import splicing as sp  # noqa: E402


def test_position_is_counted_from_the_annotated_site_on_each_strand():
    plus = sp.decode_splicevault_prediction("Top1:CD:-4:2:out-of-frame", 1000, "+")
    minus = sp.decode_splicevault_prediction("Top1:CD:-4:2:out-of-frame", 1000, "-")
    assert "4 bp upstream of the annotated splice site (pos:996)" in plus
    assert "4 bp upstream of the annotated splice site (pos:1004)" in minus     # 5'->3' on the minus strand
    assert "downstream" in sp.decode_splicevault_prediction("Top1:CA:7:1:inframe", 1000, "+")


def test_without_the_annotated_site_no_position_is_made_up():
    out = sp.decode_splicevault_prediction("Top1:CD:12:3:inframe", None, "+")
    assert "12 bp downstream of the annotated splice site" in out and "pos:" not in out
    assert "Cryptic Site at the annotated splice site" in sp.decode_splicevault_prediction("Top1:CD:0:3:inframe", 50, "+")


def test_parse_uses_site_pos_not_the_variant_position():
    df = pd.DataFrame({"SpliceVault_top_events": ["Top1:CD:-4:2:out-of-frame&Top2:ES:3:1:inframe", None],
                       "SpliceVault_site_pos": [2000, None], "SpliceVault_site_sample_count": [10, None],
                       "POS": [1990, 50], "STRAND": [-1, 1]})
    out = sp.parse_splicevault(df)
    assert "(pos:2004)" in out["SpliceVault_Predictions_Decoded"].iloc[0]
    assert "Exon 3 skipping" in out["SpliceVault_Predictions_Decoded"].iloc[0]
    assert out["SpliceVault_status"].tolist() == ["aberrant_event_detected", "not_covered"]
