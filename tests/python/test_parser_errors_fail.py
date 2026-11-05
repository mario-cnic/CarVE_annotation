"""A parsing error fails TSV_TO_TABLE instead of leaving the table without the predictor's columns (MISC-36)."""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python/parsing")))
from modules import splicing as sp  # noqa: E402


def boom(*args, **kwargs):
    raise RuntimeError("unexpected value")


@pytest.mark.parametrize("parse,frame", [
    (sp.parse_splicevault, {"SpliceVault_top_events": ["Top1:CD:-4:2:inframe"]}),
    (sp.parse_maxentscan, {"MaxEntScan_ref": [8.0], "MaxEntScan_alt": [2.0]}),
])
def test_parse_error_is_raised(monkeypatch, parse, frame):
    monkeypatch.setattr(sp, "apply_rows", boom)
    with pytest.raises(RuntimeError):
        parse(pd.DataFrame(frame))
