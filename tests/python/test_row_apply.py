"""io_qc.apply_rows / row_records behave like DataFrame.apply(axis=1) returning a pd.Series per row."""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python/parsing")))

from modules.io_qc import apply_rows, row_records  # noqa: E402


def mixed_frame():
    return pd.DataFrame({
        "i": [1, 2, 3, 4],
        "f": [0.5, np.nan, 2.0, -1.0],
        "s": ["a", None, "c+1", "."],
        "b": [True, False, True, False],
        "o": [1, "x", np.nan, 2.5],
        "unused": ["u"] * 4,
    }, index=[10, 11, 12, 13])


def describe(row):
    """Value, type and membership seen by the row function, plus mixed int/float/NaN/str outputs."""
    out = []
    for col in ("i", "f", "s", "b", "o", "missing"):
        v = row.get(col)
        out.append(f"{col}:{type(v).__name__}:{v}:{col in row}")
    i = row.get("i")
    num = np.nan if pd.isna(row.get("f")) else row.get("f") * 2
    return (";".join(out), i * 10, num, int(i) if i % 2 else np.nan, "yes" if i > 2 else "no")


COLUMNS = ["i", "f", "s", "b", "o", "missing"]
NAMES = ["desc", "int_col", "float_col", "int_or_nan", "flag"]


def reference(data):
    res = data.apply(lambda r: pd.Series(list(describe(r))), axis=1)
    res.columns = NAMES
    return res


def test_values_types_and_dtypes_match_dataframe_apply():
    data = mixed_frame()
    pd.testing.assert_frame_equal(apply_rows(data, describe, COLUMNS, NAMES), reference(data))


def test_only_present_columns_are_in_the_row():
    rows = list(row_records(mixed_frame(), ["s", "missing", "s"]))
    assert [list(r) for r in rows] == [["s"]] * 4


def test_no_columns_gives_empty_rows():
    data = mixed_frame()
    res = apply_rows(data, lambda r: (len(r), "x"), ["nope"], ["n", "x"])
    assert list(res["n"]) == [0] * 4 and list(res.index) == list(data.index)


@pytest.mark.parametrize("n", [0, 1])
def test_small_frames(n):
    data = mixed_frame().head(n)
    res = apply_rows(data, describe, COLUMNS, NAMES)
    assert len(res) == n and list(res.columns) == NAMES
    if n:
        pd.testing.assert_frame_equal(res, reference(data))
