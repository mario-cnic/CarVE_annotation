"""PVS1 veto genes come from one list, read by both build_acmg_criteria and build_priority_tier (ACMG-4)."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src/python"))
from parsing.modules import acmg, scoring  # noqa: E402

DISEASE_GENE_LIST = str(Path(__file__).resolve().parents[1] / "fixtures/mini_disease_gene_list.csv")


def lof_rows():
    return pd.DataFrame({"SYMBOL": ["MYH7", "MYBPC3"], "Consequence": ["stop_gained"] * 2,
                         "pLI_gene_value": [1.0, 1.0], "NEW_IMPACT": ["HIGH"] * 2, "QC_STATUS": ["PASS"] * 2})


def pvs1(df):
    return ["PVS1" in c for c in acmg.build_acmg_criteria(df.copy())["ACMG_CRITERIA"]]


def lof_points(df, monkeypatch):
    monkeypatch.setattr(scoring, "apply_disease_phenotype_curation", lambda data, gene_list, phenotype: data)
    out = scoring.build_priority_tier(df.copy(), DISEASE_GENE_LIST)
    return out["VARIANT_PRIORITY_SCORE"].tolist()


def test_myh7_lof_is_vetoed_in_both_places(monkeypatch):
    assert acmg.PVS1_VETO_GENES == frozenset({"MYH7"})
    assert pvs1(lof_rows()) == [False, True]
    myh7, mybpc3 = lof_points(lof_rows(), monkeypatch)
    assert mybpc3 - myh7 == 35.0          # only the non-vetoed LoF gets the LoF points


def test_both_places_follow_the_one_list(monkeypatch):
    monkeypatch.setattr(acmg, "PVS1_VETO_GENES", frozenset({"MYBPC3"}))
    assert pvs1(lof_rows()) == [True, False]
    myh7, mybpc3 = lof_points(lof_rows(), monkeypatch)
    assert myh7 - mybpc3 == 35.0


def test_explanation_and_rationale_follow_the_veto(monkeypatch):
    # ACMG-7: a vetoed LoF got "LoF (+35.0)" and a LoF rationale although no LoF points were added.
    monkeypatch.setattr(scoring, "apply_disease_phenotype_curation", lambda data, gene_list, phenotype: data)
    df = lof_rows().assign(gnomADv4_AF_grpmax_joint=0.0)
    out = scoring.build_priority_tier(df, DISEASE_GENE_LIST)
    myh7, mybpc3 = out.iloc[0], out.iloc[1]
    assert "LoF (+35.0)" not in myh7["SCORE_EXPLANATION"] and "LoF (+35.0)" in mybpc3["SCORE_EXPLANATION"]
    assert "Loss-of-Function" not in myh7["TIER_RATIONALE"] and "Loss-of-Function" in mybpc3["TIER_RATIONALE"]
