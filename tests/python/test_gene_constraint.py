"""Synthetic checks for the gene constraint table build and its per-row join by Ensembl gene ID."""
import csv
import gzip
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src/python'))
from parsing.modules.predictors import annotate_gene_constraint, load_gene_constraint
import build_gene_constraint_table as bgc

TABLE_COLUMNS = ["gene_id", "gene", "transcript", "transcript_basis", "pLI", "LOEUF", "LOEUF_decile", "lof_z",
                 "mis_z", "constraint_flags", "gene_flags"]


def write_table(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TABLE_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in TABLE_COLUMNS})
    return str(path)


@pytest.fixture
def table(tmp_path):
    return load_gene_constraint(write_table(tmp_path / "c.tsv", [
        {"gene_id": "ENSG01", "gene": "AAA", "transcript": "ENST01.3", "pLI": "9.9e-01", "LOEUF": "2.1e-01",
         "mis_z": "3.4", "constraint_flags": "outlier_syn"},
        {"gene_id": "ENSG02", "gene": "BBB", "transcript": "ENST02.1", "mis_z": "-0.5", "constraint_flags": "no_exp_lof"},
    ]))


def test_join_by_gene_keeps_missing_as_nan(table):
    df = pd.DataFrame({
        "Gene": ["ENSG01", "ENSG02", "ENSG03", None, "-"],
        "pLI_gene_value": [0.1, 0.1, 0.1, np.nan, np.nan],
        "SYMBOL": ["AAA", "BBB", "CCC", None, None],
    }, index=[9, 4, 7, 1, 2])
    out = annotate_gene_constraint(df.copy(), table)
    assert out.index.equals(df.index)
    assert out.gene_constraint_status.tolist() == ["matched", "matched", "not_in_table", "not_applicable", "not_applicable"]
    np.testing.assert_allclose(out.pLI_gene_value.to_numpy(), [0.99, np.nan, np.nan, np.nan, np.nan])
    np.testing.assert_allclose(out.LOEUF_gene_value.to_numpy(), [0.21, np.nan, np.nan, np.nan, np.nan])
    np.testing.assert_allclose(out.mis_z_gene_value.to_numpy(), [3.4, -0.5, np.nan, np.nan, np.nan])
    assert out.gene_constraint_flags.tolist() == ["outlier_syn", "no_exp_lof", None, None, None]
    assert out.gene_constraint_transcript.tolist() == ["ENST01.3", "ENST02.1", None, None, None]
    for col in ("pLI_gene_value", "LOEUF_gene_value", "mis_z_gene_value"):
        assert out[col].dtype == np.float64


def test_without_table_rows_with_a_gene_are_not_annotated():
    out = annotate_gene_constraint(pd.DataFrame({"Gene": ["ENSG01", None]}), None)
    assert out.gene_constraint_status.tolist() == ["not_annotated", "not_applicable"]
    assert out.pLI_gene_value.isna().all() and out.LOEUF_gene_value.isna().all()


def test_table_with_repeated_gene_is_refused(tmp_path):
    path = write_table(tmp_path / "dup.tsv", [{"gene_id": "ENSG01"}, {"gene_id": "ENSG01"}])
    with pytest.raises(ValueError, match="repeats gene_id"):
        load_gene_constraint(path)


SOURCE_HEADER = ["gene", "gene_id", "transcript", "canonical", "mane_select", "transcript_version", "gene_flags",
                 "constraint_flags", "lof.pLI", "lof.oe_ci.upper", "lof.oe_ci.upper_bin_decile", "lof.z_score",
                 "mis.z_score", "lof_hc_lc.oe_ci.upper"]


def source_row(**kw):
    base = dict.fromkeys(SOURCE_HEADER, "NA")
    base.update({"canonical": "false", "mane_select": "false", "gene_flags": "[]", "constraint_flags": "[]"})
    base.update(kw)
    return base


def test_build_picks_mane_then_canonical_ensembl_only(tmp_path):
    src = tmp_path / "src.tsv.bgz"
    rows = [
        source_row(gene="AAA", gene_id="ENSG01", transcript="ENST10", transcript_version="ENST10.2", canonical="true",
                   **{"lof.pLI": "0.2", "lof.oe_ci.upper": "0.9", "mis.z_score": "1.0"}),
        source_row(gene="AAA", gene_id="ENSG01", transcript="ENST11", transcript_version="ENST11.4", mane_select="true",
                   constraint_flags='["outlier_mis","outlier_syn"]',
                   **{"lof.pLI": "0.95", "lof.oe_ci.upper": "0.30", "lof.oe_ci.upper_bin_decile": "1", "mis.z_score": "3.1"}),
        source_row(gene="AAA", gene_id="1234", transcript="NM_000001", mane_select="true", **{"lof.pLI": "0.5"}),
        source_row(gene="BBB", gene_id="ENSG02", transcript="ENST20", transcript_version="ENST20.1", canonical="true",
                   gene_flags='["low_exome_coverage"]', **{"mis.z_score": "0.2"}),
        source_row(gene="CCC", gene_id="ENSG03", transcript="ENST30", **{"lof.pLI": "0.99"}),
    ]
    with gzip.open(src, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SOURCE_HEADER, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    out = tmp_path / "out.tsv"
    bgc.main(["--constraint-tsv", str(src), "--source-url", "https://example.invalid/x", "--release", "test",
              "--out-tsv", str(out)])
    got = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
    assert got.gene_id.tolist() == ["ENSG01", "ENSG02"]
    assert got.transcript.tolist() == ["ENST11.4", "ENST20.1"]
    assert got.transcript_basis.tolist() == ["mane_select", "canonical"]
    assert got.pLI.tolist() == ["0.95", ""] and got.LOEUF.tolist() == ["0.30", ""]
    assert got.constraint_flags.tolist() == ["outlier_mis;outlier_syn", ""]
    assert got.gene_flags.tolist() == ["", "low_exome_coverage"]
    manifest = json.loads(Path(str(out) + ".manifest.json").read_text())
    assert manifest["genes"] == 2 and manifest["genes_missing"]["pLI"] == 1
    assert manifest["rows_read"] == {"ensembl_rows": 4, "refseq_rows_skipped": 1}
    loaded = load_gene_constraint(str(out))
    assert loaded.loc["ENSG01", "pLI"] == pytest.approx(0.95) and np.isnan(loaded.loc["ENSG02", "LOEUF"])


def test_build_refuses_two_mane_transcripts(tmp_path):
    src = tmp_path / "src.tsv.bgz"
    with gzip.open(src, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SOURCE_HEADER, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for t in ("ENST1", "ENST2"):
            w.writerow(source_row(gene="AAA", gene_id="ENSG01", transcript=t, transcript_version=t, mane_select="true"))
    with pytest.raises(ValueError, match="two mane_select"):
        bgc.main(["--constraint-tsv", str(src), "--source-url", "u", "--release", "r", "--out-tsv", str(tmp_path / "o.tsv")])
