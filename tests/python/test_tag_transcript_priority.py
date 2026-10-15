"""TRANSCRIPT_PRIORITY_TIER: the streaming Parquet writer matches the pandas path and the per-row rule."""
import json
import os
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python")))

import tag_transcript_priority as tp  # noqa: E402
from provenance_stamp import STAMP_KEY, with_stamp  # noqa: E402

CURATED = {"ENST00000545968", "ENST00000601016"}
FEATURES = ["ENST00000545968.6", "ENST00000545968", " ENST00000545968.6 ", "ENST00000000001.1", "", None,
            "nan", "ENST00000601016.5", "ENST00000601016.5.2"]
MANES = ["NM_1.1", "", ".", " . ", None, "NM_2", "nan", "  ", "\t"]


def row_rule(row):
    """The original per-row rule of tag_transcript_priority.py."""
    feature = str(row.get("Feature", "")).strip()
    if feature and feature.split(".")[0] in CURATED:
        return 1
    mane = row.get("MANE_SELECT", "")
    if pd.notna(mane) and str(mane).strip() not in ("", "."):
        return 2
    return 3


TIERS = ["Tier 1", "Tier 2", "Tier 3", "Tier 4"]


def frame(n=200):
    rng = np.random.default_rng(1)
    tiers = ["Tier 4"] * min(n, 128) + [TIERS[i] for i in rng.integers(0, 2, max(0, n - 128))]
    return pd.DataFrame({"Locus": [f"chr1:{i}" for i in range(n)],
                         "Feature": [FEATURES[i] for i in rng.integers(0, len(FEATURES), n)],
                         "MANE_SELECT": [MANES[i] for i in rng.integers(0, len(MANES), n)],
                         "score": np.where(np.arange(n) % 5 == 0, np.nan, rng.random(n)),
                         "empty": [None] * n,
                         "PRIORITY_TIER": pd.Categorical(tiers, categories=TIERS, ordered=True)})


def test_vectorized_tiers_follow_the_row_rule():
    df = frame()
    expected = [row_rule(r) for _, r in df.iterrows()]
    assert list(tp.tag_priority(df, CURATED)[tp.TIER_COLUMN]) == expected


def test_missing_columns_and_float_mane():
    df = pd.DataFrame({"MANE_SELECT": [np.nan, 1.5, np.nan]})
    assert list(tp.tag_priority(df, CURATED)[tp.TIER_COLUMN]) == [3, 2, 3]
    assert list(tp.tag_priority(pd.DataFrame({"x": [1, 2]}), CURATED)[tp.TIER_COLUMN]) == [3, 3]


def test_streaming_parquet_equals_pandas_path(tmp_path, monkeypatch):
    monkeypatch.setattr(tp, "ROW_GROUP_ROWS", 200)
    src = tmp_path / "in.pq"
    frame(500).to_parquet(src, index=False, row_group_size=128)
    stamp = json.dumps({"run_id": "r1"})
    streamed, reference = str(tmp_path / "streamed.pq"), str(tmp_path / "reference.pq")
    counts = tp.tag_parquet(str(src), streamed, CURATED, stamp)
    tagged = tp.tag_priority(pd.read_parquet(src), CURATED)
    table = pa.Table.from_pandas(tagged, preserve_index=False)
    pq.write_table(table.replace_schema_metadata(with_stamp(table.schema, stamp).metadata), reference,
                   row_group_size=200)
    a, b = pq.ParquetFile(streamed), pq.ParquetFile(reference)
    assert a.metadata.metadata == b.metadata.metadata
    assert a.schema_arrow.equals(b.schema_arrow, check_metadata=True)
    assert a.read().equals(b.read())
    assert [a.metadata.row_group(i).num_rows for i in range(a.metadata.num_row_groups)] == [200, 200, 100]
    restored = pd.read_parquet(streamed)["PRIORITY_TIER"]
    assert list(restored.cat.categories) == TIERS and restored.cat.ordered
    assert json.loads(a.metadata.metadata[STAMP_KEY.encode()]) == {"run_id": "r1"}
    tiers = a.read().column(tp.TIER_COLUMN).to_pylist()
    assert counts == {t: tiers.count(t) for t in (1, 2, 3)}


def test_streaming_handles_an_empty_table(tmp_path):
    src = tmp_path / "in.pq"
    frame(0).astype({"Feature": object, "MANE_SELECT": object}).to_parquet(src, index=False)
    out = str(tmp_path / "out.pq")
    assert tp.tag_parquet(str(src), out, CURATED) == {1: 0, 2: 0, 3: 0}
    table = pq.read_table(out)
    assert table.num_rows == 0 and table.schema.names[-1] == tp.TIER_COLUMN
