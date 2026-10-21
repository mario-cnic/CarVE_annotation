"""Chunked TSV -> Parquet helpers for filter_variants.py: whole-file column dtypes from a first pass,
Locus-aligned chunks, and one Parquet file with a single schema written from per-chunk parts.

The output is meant to equal reading the whole TSV with pandas (`low_memory=False`) and writing it
with `DataFrame.to_parquet(index=False)`: same column dtypes on input, same rows and order, same
Arrow schema, row groups of `ROW_GROUP_ROWS` rows and the pandas footer entry.
"""
import json
import os
import shutil

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .io_qc import logger

LOCUS_KEY = ["CHROM", "POS", "REF", "ALT"]
ROW_GROUP_ROWS = 1024 * 1024  # pyarrow's default for write_table, used by DataFrame.to_parquet
BOOL_TEXT = {"True": True, "TRUE": True, "true": True, "False": False, "FALSE": False, "false": False}


def _chunk_kind(s: pd.Series) -> str:
    """dtype pandas inferred for one column of one chunk, as one of the kinds combined below."""
    if s.dtype == np.int64:
        return "int"
    if s.dtype == np.float64:
        return "nan" if s.isna().all() else "float"
    if s.dtype == np.bool_:
        return "bool"
    if s.dtype == object:
        return "objbool" if pd.api.types.infer_dtype(s, skipna=True) == "boolean" else "str"
    raise ValueError(f"column {s.name!r}: unsupported inferred dtype {s.dtype}")


def _combine(kinds: set) -> str:
    """dtype the whole-file read gives a column, from the dtypes inferred for its chunks."""
    if "str" in kinds:
        return "str"
    if kinds & {"bool", "objbool"}:
        if kinds - {"bool", "objbool", "nan"}:
            return "str"    # True/False mixed with numbers: the whole column stays text
        return "bool" if kinds == {"bool"} else "objbool"
    if kinds == {"int"}:
        return "int"
    return "float"          # int/float/all-NaN mixes and all-NaN columns read as float64


def whole_file_dtypes(path: str, read_kwargs: dict, chunk_rows: int) -> dict:
    """{column: combined kind} from a first pass over the TSV in chunks."""
    kinds = {}
    for chunk in pd.read_csv(path, chunksize=chunk_rows, **read_kwargs):
        for col in chunk.columns:
            kinds.setdefault(col, set()).add(_chunk_kind(chunk[col]))
    return {col: _combine(k) for col, k in kinds.items()}


def read_dtype_args(kinds: dict) -> tuple[dict, list]:
    """read_csv `dtype` map for the second pass, plus the columns read as text and mapped to bools."""
    dtype = {"str": str, "objbool": str, "int": np.int64, "float": np.float64, "bool": np.bool_}
    return {c: dtype[k] for c, k in kinds.items()}, [c for c, k in kinds.items() if k == "objbool"]


def locus_chunks(path: str, read_kwargs: dict, chunk_rows: int, dtype: dict, objbool_cols: list):
    """DataFrames of about `chunk_rows` rows that never split the rows of one CHROM/POS/REF/ALT;
    each has a fresh RangeIndex. Rows of one record must be consecutive in the TSV."""
    carry = None
    for chunk in pd.read_csv(path, chunksize=chunk_rows, dtype=dtype, **read_kwargs):
        for col in objbool_cols:
            chunk[col] = chunk[col].map(BOOL_TEXT, na_action="ignore").astype(object)
        if carry is not None:
            chunk = pd.concat([carry, chunk], ignore_index=True)
        last = chunk[LOCUS_KEY].iloc[-1]
        same = np.logical_and.reduce([(chunk[c] == last[c]).to_numpy() for c in LOCUS_KEY])
        start = len(same) - int(np.argmin(same[::-1])) if not same.all() else 0
        carry = chunk.iloc[start:].reset_index(drop=True)
        if start:
            yield chunk.iloc[:start].reset_index(drop=True)
    if carry is not None and len(carry):
        yield carry


class CoordinateOrder:
    """io_qc.check_coordinate_order over rows that arrive in chunks."""

    def __init__(self, chrom_col: str = "CHROM", pos_col: str = "POS"):
        self.chrom_col, self.pos_col = chrom_col, pos_col
        self.contigs, self.blocks, self.decreases, self.rows = set(), 0, 0, 0
        self.prev_chrom, self.prev_pos = None, None
        self.has_cols = True

    def add(self, data: pd.DataFrame):
        if self.chrom_col not in data.columns or self.pos_col not in data.columns:
            self.has_cols = False
            return
        if data.empty:
            return
        chrom = data[self.chrom_col].astype(str).to_numpy()
        pos = pd.to_numeric(data[self.pos_col], errors="coerce").to_numpy()
        if self.prev_chrom is not None:
            chrom = np.concatenate([[self.prev_chrom], chrom])
            pos = np.concatenate([[self.prev_pos], pos])
        new_block = np.empty(len(chrom), dtype=bool)
        new_block[0] = self.prev_chrom is None
        new_block[1:] = chrom[1:] != chrom[:-1]
        self.blocks += int(new_block.sum())
        self.decreases += int(((pos[1:] < pos[:-1]) & ~new_block[1:]).sum())
        self.contigs.update(np.unique(chrom).tolist())
        self.prev_chrom, self.prev_pos = chrom[-1], pos[-1]
        self.rows += len(data)

    def result(self) -> dict:
        if not self.rows or not self.has_cols:
            return {"checked": False, "ordered": None, "contigs": 0, "contig_blocks": 0, "position_decreases": 0}
        return {"checked": True, "ordered": self.blocks == len(self.contigs) and self.decreases == 0,
                "contigs": len(self.contigs), "contig_blocks": self.blocks, "position_decreases": self.decreases}


def _unify_type(name: str, types: list):
    """Arrow type of a column over the parts where it has values: int64 + double -> double; null when
    no part has a value."""
    found = {t for t in types if not pa.types.is_null(t)}
    if not found:
        return pa.null()
    if len(found) == 1:
        return found.pop()
    if found == {pa.int64(), pa.float64()}:
        return pa.float64()
    raise ValueError(f"column {name!r} has incompatible types across chunks: {sorted(map(str, found))}")


def _null_counts(path: str) -> tuple[dict, int]:
    """({column: null count}, rows) of a Parquet file, from the column-chunk statistics (read when missing)."""
    meta = pq.ParquetFile(path).metadata
    names = meta.schema.to_arrow_schema().names
    out = {}
    for j, name in enumerate(names):
        nulls = 0
        for i in range(meta.num_row_groups):
            stats = meta.row_group(i).column(j).statistics
            if stats is None or not stats.has_null_count:
                nulls = None
                break
            nulls += stats.null_count
        if nulls is None:
            nulls = pq.read_table(path, columns=[name]).column(0).null_count
        out[name] = nulls
    return out, meta.num_rows


def _column_order(schemas: list) -> list:
    """Columns of the widest part. A part may lack columns that a step only creates when one of its rows
    needs them (the whole frame has them, null on those rows); its other columns must keep the same order."""
    names = max((s.names for s in schemas), key=len)
    position = {n: i for i, n in enumerate(names)}
    for s in schemas:
        idx = [position.get(n) for n in s.names]
        if None in idx or idx != sorted(idx):
            raise ValueError(f"chunk columns are not a same-order subset of {len(names)} columns: "
                             f"{[n for n in s.names if n not in position][:5]}")
    return names


def _pandas_entry(metas: list, schemas: list, name: str) -> list:
    """pandas footer entries of one column, one per part that has it (None where it is absent)."""
    return [next(c for c in m["columns"] if c["name"] == name) if name in s.names else None
            for m, s in zip(metas, schemas)]


def _pandas_metadata(schemas: list, unified: pa.Schema, with_nulls: set) -> bytes:
    """pandas footer entry of the whole frame: per column, the entry of a part whose Arrow type is the
    unified one, preferring numpy dtype `object` (a column that holds NaN anywhere is object as a whole).
    An int column promoted to float64 gets a float64 entry; a bool column with nulls is object."""
    metas = [json.loads(s.metadata[b"pandas"]) for s in schemas]
    out = next(m for m, s in zip(metas, schemas) if s.names == unified.names)
    columns = []
    for field in unified:
        entries = _pandas_entry(metas, schemas, field.name)
        candidates = [e for e, s in zip(entries, schemas) if e and s.field(field.name).type == field.type]
        if candidates:
            entry = dict(next((c for c in candidates if c.get("numpy_type") == "object"), candidates[0]))
        else:
            entry = dict(next(e for e in entries if e), pandas_type="float64", numpy_type="float64")
        if field.name in with_nulls and entry.get("numpy_type") == "bool":
            entry["numpy_type"] = "object"
        columns.append(entry)
    out["columns"] = columns
    return json.dumps(out).encode()


def _conform(table: pa.Table, schema: pa.Schema) -> pa.Table:
    """The part with the unified schema: missing or all-null columns as typed nulls, others cast."""
    arrays = [table.column(f.name).cast(f.type)
              if f.name in table.column_names and table.column(f.name).null_count < table.num_rows
              else pa.nulls(table.num_rows, f.type) for f in schema]
    return pa.Table.from_arrays(arrays, schema=schema)


class PartWriter:
    """Writes each chunk's frame as a Parquet part in a temporary directory owned by this execution,
    then one output file with the unified schema; the directory is removed afterwards."""

    def __init__(self, output: str):
        self.output = output
        self.tmp = f"{output}.parts.tmp"
        if os.path.exists(self.tmp):
            raise FileExistsError(f"{self.tmp} exists: left by another execution, refusing to reuse it")
        os.makedirs(self.tmp)
        self.parts = []

    def add(self, data: pd.DataFrame):
        path = os.path.join(self.tmp, f"part_{len(self.parts):05d}.pq")
        data.to_parquet(path, index=False)
        self.parts.append(path)

    def finish(self) -> int:
        if not self.parts:
            raise ValueError("no rows to write")
        schemas = [pq.read_schema(p) for p in self.parts]
        counts = [_null_counts(p) for p in self.parts]
        empty = [{n for n, c in nulls.items() if c == rows} for nulls, rows in counts]
        names = _column_order(schemas)
        metas = [json.loads(s.metadata[b"pandas"]) for s in schemas]
        # A column with no value in a part says nothing about its type (pandas gives it float64 or object);
        # a column with no value anywhere keeps the type every part gives it.
        fields, with_nulls = [], set()
        for n in names:
            valued = [s.field(n).type for s, e in zip(schemas, empty) if n in s.names and n not in e]
            declared = {s.field(n).type for s in schemas if n in s.names}
            if not valued and len(declared) > 1:
                raise ValueError(f"column {n!r} has no value in any chunk and different types: {declared}")
            if any(n not in s.names or nulls[n] for s, (nulls, _) in zip(schemas, counts)):
                with_nulls.add(n)
            kind = _unify_type(n, valued) if valued else declared.pop()
            # numpy ints cannot hold NaN: on the whole frame a column with nulls is float64.
            numpy_ints = all(e["numpy_type"].startswith("int")
                             for e, e_empty in zip(_pandas_entry(metas, schemas, n), empty) if e and n not in e_empty)
            if pa.types.is_integer(kind) and n in with_nulls and numpy_ints:
                kind = pa.float64()
            fields.append(pa.field(n, kind))
        unified = pa.schema(fields)
        unified = unified.with_metadata({b"pandas": _pandas_metadata(schemas, unified, with_nulls)})
        bare = unified.remove_metadata()
        rows = 0
        with pq.ParquetWriter(self.output, unified) as writer:
            buffer = []
            buffered = 0
            for p in self.parts:
                table = _conform(pq.read_table(p), bare)
                buffer.append(table)
                buffered += table.num_rows
                while buffered >= ROW_GROUP_ROWS:
                    whole = pa.concat_tables(buffer)
                    writer.write_table(whole.slice(0, ROW_GROUP_ROWS).replace_schema_metadata(unified.metadata),
                                       row_group_size=ROW_GROUP_ROWS)
                    rest = whole.slice(ROW_GROUP_ROWS)
                    buffer, buffered = [rest], rest.num_rows
                    rows += ROW_GROUP_ROWS
            if buffered:
                writer.write_table(pa.concat_tables(buffer).replace_schema_metadata(unified.metadata),
                                   row_group_size=ROW_GROUP_ROWS)
                rows += buffered
        logger.info(f"Wrote {rows} rows from {len(self.parts)} chunk(s) to {self.output}")
        return rows

    def cleanup(self):
        shutil.rmtree(self.tmp, ignore_errors=True)
