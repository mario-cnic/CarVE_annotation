"""filter_variants.py CLI: --freq 0 is a filter, not 'no filter'; unused genotype options are gone."""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python/parsing")))
import filter_variants as fv  # noqa: E402

GENES = os.path.join(os.path.dirname(__file__), "../fixtures/mini_disease_gene_list.csv")


def run_transform(monkeypatch, extra):
    calls = []
    monkeypatch.setattr(fv, "filter_freq", lambda data, freq, freq_column: calls.append(freq) or data)
    args = fv.argparsing().parse_args(["--input", "x.tsv", "--output", "y.pq", "--disease-gene-list", GENES,
                                       "--skip_quality_filter", *extra])
    df = pd.DataFrame({"CHROM": ["chr1"], "POS": [100], "REF": ["A"], "ALT": ["G"], "SYMBOL": ["MYH7"],
                       "HGVSc": ["ENST00000355349.4:c.100G>A"], "Feature": ["ENST00000355349"], "MANE_SELECT": ["NM_000257.4"]})
    fv.transform(df, args, {}, None, [])
    return calls


def test_freq_zero_filters(monkeypatch):
    assert run_transform(monkeypatch, ["--freq", "0"]) == [0.0]
    assert run_transform(monkeypatch, []) == []


@pytest.mark.parametrize("option", [["--genotype_filter"], ["--genotype_group", "x"]])
def test_unused_genotype_options_removed(option):
    with pytest.raises(SystemExit):
        fv.argparsing().parse_args(["--disease-gene-list", GENES, *option])
