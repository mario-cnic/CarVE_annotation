"""Chunked TSV -> Parquet (modules/chunked.py) gives what the whole-file read and to_parquet give."""
import json
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python/parsing")))

from modules import chunked  # noqa: E402
from modules.io_qc import check_coordinate_order  # noqa: E402

READ = dict(na_values=".", sep="\t", header=0, index_col=False, low_memory=False)


def write_tsv(path, rows):
    header = ["CHROM", "POS", "REF", "ALT", "i", "i_nan_late", "f", "nan_then_text", "num_then_text",
              "flag", "flag_nan_late"]
    lines = ["\t".join(header)] + ["\t".join(map(str, r)) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


ROWS = [
    ("chr1", 10, "A", "C", 1, 5, 0.5, ".", 7, "True", "True"),
    ("chr1", 10, "A", "C", 2, 6, ".", ".", 8, "False", "False"),
    ("chr1", 20, "G", "T", 3, 7, 1.5, ".", 9, "True", "True"),
    ("chr1", 30, "T", "A", 4, ".", 2.0, "x", "y", "False", "."),
    ("chr1", 30, "T", "A", 5, 9, ".", "z", 10, "True", "False"),
    ("chr2", 5, "C", "G", 6, 10, 3.0, ".", 11, "False", "True"),
]


def chunked_frame(path, chunk_rows):
    kinds = chunked.whole_file_dtypes(path, READ, chunk_rows)
    dtype, objbool = chunked.read_dtype_args(kinds)
    return list(chunked.locus_chunks(path, READ, chunk_rows, dtype, objbool))


@pytest.mark.parametrize("chunk_rows", [1, 2, 3, 100])
def test_chunks_read_like_the_whole_file(tmp_path, chunk_rows):
    path = write_tsv(tmp_path / "t.tsv", ROWS)
    whole = pd.read_csv(path, **READ)
    parts = chunked_frame(path, chunk_rows)
    pd.testing.assert_frame_equal(pd.concat(parts, ignore_index=True), whole)


def test_chunks_never_split_a_variant(tmp_path):
    path = write_tsv(tmp_path / "t.tsv", ROWS)
    parts = chunked_frame(path, 1)
    keys = [tuple(p[chunked.LOCUS_KEY].drop_duplicates().itertuples(index=False)) for p in parts]
    flat = [k for ks in keys for k in ks]
    assert len(flat) == len(set(flat)) == 4
    assert [len(p) for p in parts] == [2, 1, 2, 1]


def test_whole_file_dtype_rules():
    assert chunked._combine({"int"}) == "int"
    assert chunked._combine({"int", "nan"}) == "float"
    assert chunked._combine({"int", "float"}) == "float"
    assert chunked._combine({"nan"}) == "float"
    assert chunked._combine({"float", "str"}) == "str"
    assert chunked._combine({"bool"}) == "bool"
    assert chunked._combine({"bool", "nan"}) == "objbool"
    assert chunked._combine({"bool", "int"}) == "str"


def to_parquet_whole(frames, path):
    whole = pd.concat(frames, ignore_index=True)
    whole.to_parquet(path, index=False)
    return path


def test_part_writer_equals_whole_frame_to_parquet(tmp_path):
    part1 = pd.DataFrame({"a": [1, 2], "s": [None, None], "t": [np.nan, np.nan]})
    part2 = pd.DataFrame({"a": [1.5, np.nan], "s": ["x", None], "t": ["y", None], "m": [0.1, np.nan]})
    out = str(tmp_path / "out.pq")
    writer = chunked.PartWriter(out)
    try:
        writer.add(part1)
        writer.add(part2)
        assert writer.finish() == 4
    finally:
        writer.cleanup()
    assert not os.path.exists(out + ".parts.tmp")
    ref = to_parquet_whole([part1.assign(m=np.nan)[["a", "s", "t", "m"]], part2], str(tmp_path / "ref.pq"))
    got, exp = pq.read_table(out), pq.read_table(ref)
    assert got.schema.remove_metadata() == exp.schema.remove_metadata()
    assert got.equals(exp)
    meta = lambda t: [(c["name"], c["pandas_type"], c["numpy_type"]) for c in json.loads(t.schema.metadata[b"pandas"])["columns"]]
    assert meta(got) == meta(exp)


@pytest.mark.parametrize("filled,other", [
    (pd.DataFrame({"a": [1, 2], "n": [21, 22]}), pd.DataFrame({"a": [3, 4], "n": [np.nan, np.nan]})),
    (pd.DataFrame({"a": [1, 2], "n": [21, 22]}), pd.DataFrame({"a": [3, 4]})),
    (pd.DataFrame({"a": [1, 2], "n": [True, False]}), pd.DataFrame({"a": [3, 4], "n": [np.nan, np.nan]})),
    (pd.DataFrame({"a": [1, 2], "n": [True, False]}), pd.DataFrame({"a": [3, 4]})),
])
def test_part_writer_column_filled_in_one_part_only(tmp_path, filled, other):
    """A column a step fills on every row of one part and leaves empty or out of another: the whole
    frame holds NaN there, so ints become float64 and bools become object, as in pandas."""
    for parts in ([filled, other], [other, filled]):
        out = str(tmp_path / "out.pq")
        writer = chunked.PartWriter(out)
        try:
            for part in parts:
                writer.add(part)
            writer.finish()
        finally:
            writer.cleanup()
        # The step run on the whole frame builds each column from the rows' values.
        whole = pd.DataFrame({c: [v for part in parts for v in (part[c].tolist() if c in part else [np.nan] * len(part))]
                              for c in ["a", "n"]})
        ref = str(tmp_path / "ref.pq")
        whole.to_parquet(ref, index=False)
        got, exp = pq.read_table(out), pq.read_table(ref)
        assert got.schema.remove_metadata() == exp.schema.remove_metadata()
        assert got.equals(exp)
        meta = lambda t: [(c["name"], c["pandas_type"], c["numpy_type"]) for c in json.loads(t.schema.metadata[b"pandas"])["columns"]]
        assert meta(got) == meta(exp)
        os.remove(out)


def test_part_writer_keeps_nullable_ints(tmp_path):
    parts = [pd.DataFrame({"n": pd.array([1, 2], dtype="Int64")}), pd.DataFrame({"n": pd.array([None, None], dtype="Int64")})]
    out = str(tmp_path / "out.pq")
    writer = chunked.PartWriter(out)
    try:
        for part in parts:
            writer.add(part)
        writer.finish()
    finally:
        writer.cleanup()
    exp = pq.read_table(to_parquet_whole(parts, str(tmp_path / "ref.pq")))
    got = pq.read_table(out)
    assert got.schema.remove_metadata() == exp.schema.remove_metadata() and got.equals(exp)
    assert pd.read_parquet(out)["n"].dtype == "Int64"


def test_part_writer_row_groups(tmp_path, monkeypatch):
    monkeypatch.setattr(chunked, "ROW_GROUP_ROWS", 3)
    out = str(tmp_path / "out.pq")
    writer = chunked.PartWriter(out)
    try:
        for start in (0, 2, 4):
            writer.add(pd.DataFrame({"a": [start, start + 1]}))
        writer.finish()
    finally:
        writer.cleanup()
    meta = pq.ParquetFile(out).metadata
    assert [meta.row_group(i).num_rows for i in range(meta.num_row_groups)] == [3, 3]
    assert pq.read_table(out).column("a").to_pylist() == list(range(6))


def test_part_writer_refuses_reordered_columns(tmp_path):
    writer = chunked.PartWriter(str(tmp_path / "out.pq"))
    try:
        writer.add(pd.DataFrame({"a": [1], "b": [2]}))
        writer.add(pd.DataFrame({"b": [3], "a": [4]}))
        with pytest.raises(ValueError):
            writer.finish()
    finally:
        writer.cleanup()


@pytest.mark.parametrize("chrom,pos", [
    (["chr1"] * 4 + ["chr2"] * 3, [1, 2, 2, 5, 1, 3, 4]),
    (["chr1", "chr1", "chr2", "chr1", "chr2"], [1, 2, 1, 3, 2]),
    (["chr1"] * 5, [5, 4, 6, 7, 1]),
])
def test_coordinate_order_in_chunks_equals_whole(chrom, pos):
    data = pd.DataFrame({"CHROM": chrom, "POS": pos})
    for size in (1, 2, 3, len(data)):
        order = chunked.CoordinateOrder()
        for start in range(0, len(data), size):
            order.add(data.iloc[start:start + size])
        assert order.result() == check_coordinate_order(data)
