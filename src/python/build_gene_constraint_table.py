#!/usr/bin/env python3
"""Build resources/gnomad_v4.1.2_gene_constraint.grch38.tsv + a build manifest.

One row per Ensembl gene (ENSG) from the gnomAD v4.1.2 constraint metrics TSV (GENCODE v39
transcripts, GRCh38): the MANE Select Ensembl transcript, or the canonical one when the gene has
no MANE Select. RefSeq rows are dropped (their gene_id is an Entrez ID).

Columns: pLI = lof.pLI, LOEUF = lof.oe_ci.upper, LOEUF_decile = lof.oe_ci.upper_bin_decile and
lof_z = lof.z_score (high-confidence pLoF), mis_z = mis.z_score, plus gnomAD's constraint_flags
and gene_flags. Values of flagged transcripts are kept; a missing metric stays empty.
"""
import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone

SOURCE_COLUMNS = {
    "pLI": "lof.pLI",
    "LOEUF": "lof.oe_ci.upper",
    "LOEUF_decile": "lof.oe_ci.upper_bin_decile",
    "lof_z": "lof.z_score",
    "mis_z": "mis.z_score",
}
OUT_COLUMNS = ["gene_id", "gene", "transcript", "transcript_basis", *SOURCE_COLUMNS,
               "constraint_flags", "gene_flags"]
MISSING = {"", "NA", "nan", "NaN", "None"}


def file_hashes(path):
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            sha.update(chunk)
            md5.update(chunk)
    return sha.hexdigest(), md5.hexdigest()


def flag_text(value):
    """gnomAD writes flags as a JSON list; keep them as a ';'-joined string, empty when none."""
    value = (value or "").strip()
    if value in MISSING or value == "[]":
        return ""
    return ";".join(json.loads(value))


def select_rows(reader):
    """Pick one Ensembl transcript per gene: MANE Select, else canonical."""
    chosen, counts = {}, Counter()
    for row in reader:
        if not row["transcript"].startswith("ENST"):
            counts["refseq_rows_skipped"] += 1
            continue
        counts["ensembl_rows"] += 1
        basis = "mane_select" if row["mane_select"] == "true" else "canonical" if row["canonical"] == "true" else None
        if basis is None:
            continue
        gene_id = row["gene_id"]
        if not gene_id.startswith("ENSG"):
            raise ValueError(f"Ensembl transcript {row['transcript']} has a non-ENSG gene_id {gene_id!r}")
        previous = chosen.get(gene_id)
        if previous is not None:
            if previous[0] == basis:
                raise ValueError(f"{gene_id}: two {basis} Ensembl transcripts ({previous[1]['transcript']}, {row['transcript']})")
            if previous[0] == "mane_select":
                continue
        chosen[gene_id] = (basis, row)
    return chosen, counts


def out_row(basis, row):
    out = {"gene_id": row["gene_id"], "gene": row["gene"], "transcript": row["transcript_version"] or row["transcript"],
           "transcript_basis": basis,
           "constraint_flags": flag_text(row["constraint_flags"]), "gene_flags": flag_text(row["gene_flags"])}
    for name, column in SOURCE_COLUMNS.items():
        value = (row[column] or "").strip()
        out[name] = "" if value in MISSING else value
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--constraint-tsv", required=True, help="gnomad.v4.1.2.constraint_metrics.tsv.bgz")
    p.add_argument("--source-url", required=True)
    p.add_argument("--release", required=True, help="e.g. gnomAD v4.1.2")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--note", action="append", default=[], help="free-text provenance note for the manifest")
    a = p.parse_args(argv)

    with gzip.open(a.constraint_tsv, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in ["gene", "gene_id", "transcript", "transcript_version", "canonical", "mane_select",
                               "constraint_flags", "gene_flags", *SOURCE_COLUMNS.values()] if c not in reader.fieldnames]
        if missing:
            raise SystemExit(f"constraint TSV lacks columns: {missing}")
        chosen, counts = select_rows(reader)

    rows = sorted((out_row(basis, row) for basis, row in chosen.values()), key=lambda r: r["gene_id"])
    with open(a.out_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    sha, md5 = file_hashes(a.constraint_tsv)
    manifest = {
        "built": datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds"),
        "source_url": a.source_url, "release": a.release,
        "source": a.constraint_tsv, "source_sha256": sha, "source_md5": md5,
        "assembly": "GRCh38", "source_transcripts": "GENCODE v39 (gnomAD constraint)",
        "fields": SOURCE_COLUMNS,
        "genes": len(rows),
        "genes_by_basis": dict(Counter(r["transcript_basis"] for r in rows)),
        "genes_missing": {name: sum(1 for r in rows if r[name] == "") for name in SOURCE_COLUMNS},
        "genes_with_constraint_flags": sum(1 for r in rows if r["constraint_flags"]),
        "rows_read": dict(counts),
        "notes": a.note,
        "output": a.out_tsv, "output_sha256": file_hashes(a.out_tsv)[0],
        "argv": sys.argv,
    }
    with open(a.out_tsv + ".manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
