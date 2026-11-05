import os
import sys

import numpy as np
import pandas as pd

# Repo parsing package only; config.py would put shared/utils/src (its own filter_variants/modules) first.
PIPELINE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(PIPELINE_ROOT, "src", "python", "parsing"))

from filter_variants import parse_branchpointer, NEW_BRANCHPOINT_COLUMNS, build_newImpact  # noqa: E402


def test_labranchor_positional_columns_only():
    df = pd.DataFrame({"LaBranchoR_score": [0.95, np.nan, "0.40"], "LaBranchoR_acc_dist": [24, np.nan, "31"]})
    res = parse_branchpointer(df.copy())
    assert list(NEW_BRANCHPOINT_COLUMNS.values()) == ["LaBranchoR_score", "LaBranchoR_acc_dist", "LaBranchoR_status"]
    assert res["LaBranchoR_status"].tolist() == ["scored", "not_covered", "scored"]
    assert res.loc[2, "LaBranchoR_score"] == 0.40 and res.loc[2, "LaBranchoR_acc_dist"] == 31
    assert not any(c.startswith(("Branchpointer_", "Branchpoint_")) for c in res.columns)


def test_parse_branchpointer_missing_columns():
    df = pd.DataFrame({"Locus": ["chr1:100-A-G"]})
    res = parse_branchpointer(df.copy())
    assert list(res.columns) == ["Locus"]


def test_labranchor_score_alone_does_not_raise_new_impact():
    # PRED-18: a positional LaBranchoR hit (site score >= 0.5) made NEW_IMPACT HIGH whatever the alleles.
    df = pd.DataFrame({"IMPACT": ["MODIFIER", "LOW"], "LaBranchoR_score": [0.95, 0.80], "Locus": ["chr1:1-A-G", "chr1:2-C-T"]})
    out = build_newImpact(parse_branchpointer(df))
    assert out["NEW_IMPACT"].astype(str).tolist() == ["MODIFIER", "LOW"]


def test_labranchor_ism_columns_and_status():
    from modules.splicing import LABRANCHOR_ISM_COLUMNS
    df = pd.DataFrame({"LaBranchoR_score": [np.nan] * 3, "LaBranchoR_acc_dist": [np.nan] * 3,
                       "LaBranchoR_ISM_d_best": ["-0.8", np.nan, np.nan], "LaBranchoR_ISM_d_max": ["-0.3", np.nan, np.nan],
                       "LaBranchoR_ISM_dist": ["30", np.nan, np.nan], "LaBranchoR_ISM_n_3pss": ["2", np.nan, np.nan],
                       "Locus": ["chr1:1001-A-G", "chr1:5-C-T", "chr1:9-TA-T"]})
    res = parse_branchpointer(df)
    assert res["LaBranchoR_ISM_status"].tolist() == ["scored", "not_covered", "not_applicable"]
    assert res.loc[0, "LaBranchoR_ISM_d_best"] == -0.8 and res.loc[0, "LaBranchoR_ISM_dist"] == 30
    assert all(c in res.columns for c in LABRANCHOR_ISM_COLUMNS)
    out = build_newImpact(res.assign(IMPACT="MODIFIER"))
    assert (out["NEW_IMPACT"].astype(str) == "MODIFIER").all()           # informational only
