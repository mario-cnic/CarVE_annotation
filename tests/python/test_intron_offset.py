"""intron_offset_signed / splice_side from HGVS (PRED-8)."""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python/parsing")))
from modules import splicing as sp  # noqa: E402

CASES = [
    # unchanged by the fix
    ("ENST00000355349.4:c.100+5G>A", None, (5, "donor")),
    ("ENST00000355349.4:c.101-2A>G", None, (-2, "acceptor")),
    ("ENST00000355349.4:c.100G>A", None, (0, "exonic")),
    ("", None, (0, "exonic")),
    # PRED-8: UTR introns, non-coding transcripts, range starting in the exon
    ("ENST00000355349.4:c.-45+3G>A", None, (3, "donor")),        # was acceptor -3 (the leading '-')
    ("ENST00000355349.4:c.-45-3G>A", None, (-3, "acceptor")),
    ("ENST00000355349.4:c.*12+2T>C", None, (2, "donor")),
    ("ENST00000355349.4:c.*12-7T>C", None, (-7, "acceptor")),
    ("ENST00000469530.1:n.123+5G>A", None, (5, "donor")),
    ("ENST00000355349.4:c.99_100+2del", None, (2, "donor")),
    ("ENST00000355349.4:c.-45_-44del", None, (0, "exonic")),
    # HGVS_OFFSET is VEP's 3' shift of the notation, not an intron position
    ("ENST00000355349.4:c.-45delA", 2, (0, "exonic")),            # was acceptor -2
    ("ENST00000355349.4:c.120dupT", 1, (0, "exonic")),            # was donor 1
]


@pytest.mark.parametrize("hgvsc,hgvs_offset,expected", CASES)
def test_offset_and_side(hgvsc, hgvs_offset, expected):
    df = pd.DataFrame({"HGVSc": [hgvsc], "HGVS_OFFSET": [hgvs_offset]})
    out = sp.build_signed_intron_offset(df)
    assert (int(out["intron_offset_signed"].iloc[0]), out["splice_side"].iloc[0]) == expected


def test_cdna_name_used_when_hgvsc_missing():
    out = sp.build_signed_intron_offset(pd.DataFrame({"HGVSc": [None], "CDNA_NAME": ["NM_000256.3:c.1624-1G>A"]}))
    assert (out["intron_offset_signed"].iloc[0], out["splice_side"].iloc[0]) == (-1, "acceptor")
