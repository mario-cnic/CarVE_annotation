"""Provenance stamp in the file metadata (footer) of a published Parquet table.

The stamp is the JSON written at launch to pipeline_info/TABLE_STAMP.json (run id, code commit,
parameters, database releases, tool versions). It is stored under one key-value entry, so it
travels with the .pq file; existing entries (ARROW:schema, pandas) are kept.
"""
import json
import sys
from pathlib import Path

STAMP_KEY = "carve.provenance"


def load_stamp(path):
    """Compact JSON text of the stamp file, or None when no path is given. Exits when a given file is unreadable."""
    if not path:
        return None
    try:
        return json.dumps(json.loads(Path(path).read_text()), separators=(",", ":"), sort_keys=True)
    except (OSError, ValueError) as e:
        sys.exit(f"ERROR: cannot read the provenance stamp {path}: {e}")


def with_stamp(schema, stamp):
    """The Arrow schema with the stamp added to its file metadata; unchanged when there is no stamp."""
    if stamp is None:
        return schema
    metadata = dict(schema.metadata or {})
    metadata[STAMP_KEY.encode()] = stamp.encode()
    return schema.with_metadata(metadata)
