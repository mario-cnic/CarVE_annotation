#!/usr/bin/env python3
"""
Adds TRANSCRIPT_PRIORITY_TIER to a VCF_TO_TABLE output table (whole-VCF entry point only —
the legacy per-gene path already filters to one transcript via --gene_set, so has no use for
this tagging step).

Tier 1: row's Feature (transcript) matches its gene's curated transcript in
        resources/gene_transcript_mapping.txt (the 217-panel's clinically-vetted choice).
Tier 2: row's own MANE_SELECT column (a standard VEP/CSQ field, already present when
        vcf_parser_pysam.py is run with --add_vep) is non-empty.
Tier 3: everything else — kept, not dropped, mirroring VEP's own --flag_pick semantics.

Transcript comparison ignores version suffixes (ENST00000545968.6 vs ENST00000545968), same
convention as src/hpc/gene_coords.sh's "drop version number from transcript".

Parquet input is processed one row group at a time with PyArrow, so memory stays at about one row
group (5.5-7.7 GB for a 1 Mi-row group of the F2005 table); the output keeps the input schema and row order, with the tier column appended, in row groups
of the same size pq.write_table uses.
"""

import argparse
import json
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from provenance_stamp import load_stamp, with_stamp

TIER_COLUMN = "TRANSCRIPT_PRIORITY_TIER"
ROW_GROUP_ROWS = 1024 * 1024  # pq.write_table's default maximum row-group length


def load_curated_mapping(path: str) -> set[str]:
    """Returns the set of curated transcript IDs (version stripped) from the 217-panel's Gen,NM,ENST file.

    Keyed by ENST on purpose: the file's gene symbols are not reliable keys (TAZ and ASNA1 are renamed to
    TAFAZZIN / GET3 in the VEP annotation, `NKX25` is a typo duplicate of NKX2-5, TAZ is listed twice), so a
    symbol-keyed lookup silently never matches those genes. An ENST identifies its gene unambiguously.
    """
    curated = set()
    with open(path) as f:
        next(f)  # header
        for line in f:
            fields = line.strip().split(",")
            if len(fields) >= 3 and fields[2].strip():
                curated.add(_strip_version(fields[2]))
    return curated


def read_table(path: str) -> pd.DataFrame:
    if path.endswith(".pq") or path.endswith(".parquet"):
        return pd.read_parquet(path)
    if path.endswith(".tsv"):
        return pd.read_csv(path, sep="\t")
    raise ValueError(f"Unsupported table format: {path}")


def write_table(df: pd.DataFrame, path: str, stamp=None) -> None:
    if path.endswith(".pq") or path.endswith(".parquet"):
        table = pa.Table.from_pandas(df, preserve_index=False)
        pq.write_table(table.replace_schema_metadata(with_stamp(table.schema, stamp).metadata), path)
    elif path.endswith(".tsv"):
        df.to_csv(path, sep="\t", index=False)
    else:
        raise ValueError(f"Unsupported table format: {path}")


def _strip_version(transcript_id: str) -> str:
    return transcript_id.split(".")[0]


def _as_text(values: pd.Series) -> pd.Series:
    """str() of every value, as the per-row rule applies it."""
    return values.map(str)


def tag_priority(df: pd.DataFrame, curated: set[str]) -> pd.DataFrame:
    """Copy of df with TRANSCRIPT_PRIORITY_TIER (1/2/3) appended."""
    feature = _as_text(df["Feature"]).str.strip() if "Feature" in df else pd.Series("", index=df.index)
    tier1 = feature.ne("") & feature.str.split(".").str[0].isin(curated)
    if "MANE_SELECT" in df:
        mane = df["MANE_SELECT"]
        tier2 = mane.notna() & ~_as_text(mane).str.strip().isin(["", "."])
    else:
        tier2 = pd.Series(False, index=df.index)
    df = df.copy()
    df[TIER_COLUMN] = np.where(tier1, 1, np.where(tier2, 2, 3)).astype("int64")
    return df


def _present(values: pa.ChunkedArray) -> pa.ChunkedArray:
    """Not null and, for floats, not NaN (pandas notna)."""
    present = pc.is_valid(values)
    if pa.types.is_floating(values.type):
        present = pc.and_(present, pc.invert(pc.fill_null(pc.is_nan(values), True)))
    return present


def _trimmed_text(values: pa.ChunkedArray) -> pa.ChunkedArray:
    return pc.utf8_trim_whitespace(pc.cast(values, pa.string()))


def tier_array(table: pa.Table, curated: pa.Array) -> pa.ChunkedArray:
    """Tier per row of an Arrow table, same rule as tag_priority."""
    names = table.schema.names
    n = table.num_rows
    if "Feature" in names:
        feature = table.column("Feature")
        text = _trimmed_text(feature)
        base = pc.list_element(pc.split_pattern(text, "."), 0)
        tier1 = pc.and_(_present(feature), pc.and_(pc.not_equal(text, ""), pc.is_in(base, curated)))
        tier1 = pc.fill_null(tier1, False)
    else:
        tier1 = pa.array(np.zeros(n, dtype=bool))
    if "MANE_SELECT" in names:
        mane = table.column("MANE_SELECT")
        empty = pc.is_in(_trimmed_text(mane), pa.array(["", "."]))
        tier2 = pc.fill_null(pc.and_(_present(mane), pc.invert(empty)), False)
    else:
        tier2 = pa.array(np.zeros(n, dtype=bool))
    return pc.if_else(tier1, 1, pc.if_else(tier2, 2, 3)).cast(pa.int64())


def _output_metadata(schema: pa.Schema, stamp) -> dict:
    """Input file metadata with the tier column added to the pandas entry and the stamp set."""
    metadata = dict(schema.metadata or {})
    metadata.pop(b"ARROW:schema", None)
    if b"pandas" in metadata:
        # Same entry pa.Table.from_pandas(df, preserve_index=False) writes, plus the tier column.
        source = json.loads(metadata[b"pandas"])
        columns = [c for c in source["columns"] if c.get("field_name") is not None
                   and not str(c["field_name"]).startswith("__index_level_")]
        columns.append({"name": TIER_COLUMN, "field_name": TIER_COLUMN, "pandas_type": "int64",
                        "numpy_type": "int64", "metadata": None})
        pandas_meta = {"index_columns": [], "column_indexes": [], "columns": columns,
                       "attributes": source.get("attributes", {}),
                       "creator": {"library": "pyarrow", "version": pa.__version__},
                       "pandas_version": pd.__version__}
        metadata[b"pandas"] = json.dumps(pandas_meta).encode()
    return with_stamp(pa.schema([], metadata=metadata), stamp).metadata


def tag_parquet(input_path: str, output_path: str, curated: set[str], stamp=None) -> dict[int, int]:
    """Streams input_path to output_path with the tier column appended; returns rows per tier."""
    source = pq.ParquetFile(input_path)
    schema = source.schema_arrow
    out_schema = pa.schema(list(schema.remove_metadata()) + [pa.field(TIER_COLUMN, pa.int64())],
                           metadata=_output_metadata(schema, stamp))
    curated_array = pa.array(sorted(curated), pa.string())
    counts = {1: 0, 2: 0, 3: 0}
    pending = []
    with pq.ParquetWriter(output_path, out_schema) as writer:
        # Output row groups of ROW_GROUP_ROWS rows, as pq.write_table cuts a whole table.
        for i in range(source.num_row_groups):
            group = source.read_row_group(i)
            tier = tier_array(group, curated_array)
            for value, count in zip(*np.unique(tier.to_numpy(), return_counts=True)):
                counts[int(value)] += int(count)
            pending.append(pa.Table.from_arrays(group.columns + [tier], schema=out_schema))
            del group, tier
            buffered = pa.concat_tables(pending)
            pending = []
            full = (buffered.num_rows // ROW_GROUP_ROWS) * ROW_GROUP_ROWS
            if full:
                writer.write_table(buffered.slice(0, full), row_group_size=ROW_GROUP_ROWS)
            # Only a partial group is carried over; an aligned input keeps one row group in memory.
            if buffered.num_rows > full:
                pending = [buffered.slice(full)]
            del buffered
        if pending:
            writer.write_table(pa.concat_tables(pending), row_group_size=ROW_GROUP_ROWS)
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--gene-transcript-mapping", default="resources/gene_transcript_mapping.txt")
    parser.add_argument("--stamp", help="Provenance stamp JSON to store in the Parquet file metadata")
    args = parser.parse_args()
    stamp = load_stamp(args.stamp)

    curated = load_curated_mapping(args.gene_transcript_mapping)
    print(f"Loaded {len(curated)} curated transcripts from {args.gene_transcript_mapping}")

    is_parquet = args.input.endswith((".pq", ".parquet")) and args.output.endswith((".pq", ".parquet"))
    if is_parquet:
        counts = tag_parquet(args.input, args.output, curated, stamp)
        tier_counts = {k: v for k, v in counts.items() if v}
        rows = sum(counts.values())
    else:
        df = tag_priority(read_table(args.input), curated)
        tier_counts = dict(df[TIER_COLUMN].value_counts().sort_index())
        rows = len(df)
        write_table(df, args.output, stamp)
    print(f"Tier distribution: {tier_counts}")
    print(f"Wrote {rows} rows to {args.output}")


if __name__ == "__main__":
    sys.exit(main())
