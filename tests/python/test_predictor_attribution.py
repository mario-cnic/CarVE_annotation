"""Predictor scores must land only on the row of their own gene / transcript.

Synthetic data only. Multi-entry values use '&' as in the real table (',' in the raw VCF INFO).
Gene IDs/symbols are public GENCODE/HGNC identifiers used as labels, not patient data.
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python/parsing")))

from modules import gene_identity as gi  # noqa: E402
from modules import splicing as sp  # noqa: E402

PMS2, AIMP2 = "ENSG00000122512", "ENSG00000106305"
MADD, MYBPC3 = "ENSG00000110514", "ENSG00000134571"
PRKN = "ENSG00000185345"
OTHER = "ENSG00000999999"


@pytest.fixture
def identity(tmp_path):
    hgnc = tmp_path / "hgnc.tsv"
    hgnc.write_text(
        "hgnc_id\tsymbol\tprev_symbol\talias_symbol\tensembl_gene_id\n"
        f"HGNC:1\tPMS2\t\t\t{PMS2}\n"
        f"HGNC:2\tAIMP2\t\t\t{AIMP2}\n"
        f"HGNC:3\tPRKN\tPARK2\t\t{PRKN}\n"
        f"HGNC:4\tMADD\t\t\t{MADD}\n"
        f"HGNC:5\tMYBPC3\t\t\t{MYBPC3}\n"
        f"HGNC:6\tDISAGREE\t\t\tENSG00000000006\n"
    )
    sym = tmp_path / "map.tsv"
    sym.write_text(
        "symbol\tensg\toverlap_bp\trunner_up_ensg\trunner_up_bp\n"
        f"PMS2\t{PMS2}\t1\t\t0\nAIMP2\t{AIMP2}\t1\t\t0\nPARK2\t{PRKN}\t1\t\t0\n"
        f"MADD\t{MADD}\t1\t\t0\nMYBPC3\t{MYBPC3}\t1\t\t0\n"
        "DISAGREE\tENSG00000000007\t1\t\t0\nCOORDONLY\tENSG00000000008\t1\t\t0\n"
    )
    cur = tmp_path / "curated.txt"
    cur.write_text("Gen,NM,ENST\nGENEX,NM_000001.2,ENST00000000011\n")
    return gi.GeneIdentity(str(sym), str(hgnc), str(cur))


def sai(symbol, ds, allele="A", dp=1):
    return f"{allele}|{symbol}|{ds}|0.00|0.00|0.00|{dp}|{dp}|{dp}|{dp}"


def one(df, col):
    return df[col].iloc[0]


# ------------------------------------------------------------------ SpliceAI
def test_spliceai_overlapping_genes_each_row_gets_its_own_entry(identity):
    value = sai("AIMP2", "0.90") + "&" + sai("PMS2", "0.05")
    df = pd.DataFrame({
        "Locus": ["chr7:1-C-A"] * 4, "SpliceAI": [value] * 4,
        "Gene": [PMS2, AIMP2, OTHER, None], "SYMBOL": ["PMS2", "AIMP2", "OTHERG", None],
    })
    out = sp.parse_spliceai_custom(df, identity=identity)
    assert out["spliceai_custom_MAX"].iloc[0] == pytest.approx(0.05)   # PMS2 row: PMS2's score, not AIMP2's
    assert out["spliceai_custom_MAX"].iloc[1] == pytest.approx(0.90)
    assert np.isnan(out["spliceai_custom_MAX"].iloc[2]) and out["spliceai_custom_match"].iloc[2] == gi.NO_ENTRY
    assert np.isnan(out["spliceai_custom_MAX"].iloc[3]) and out["spliceai_custom_match"].iloc[3] == gi.NO_ROW_GENE
    assert np.allclose(out["spliceai_custom_anygene_MAX"], 0.90)  # variant-level max kept separately
    assert list(out["spliceai_custom_match"].iloc[:2]) == [gi.MATCHED, gi.MATCHED]


def test_spliceai_multi_gene_value_is_no_longer_dropped(identity):
    """The old parser split on ',' only and returned NaN for every '&'-joined multi-gene value."""
    df = pd.DataFrame({"Locus": ["chr7:1-C-A"], "SpliceAI": [sai("PMS2", "0.30") + "&" + sai("AIMP2", "0.10")],
                       "Gene": [PMS2], "SYMBOL": ["PMS2"]})
    assert one(sp.parse_spliceai_custom(df, identity=identity), "spliceai_custom_MAX") == pytest.approx(0.30)


def test_spliceai_alias_symbol_resolved_by_two_routes(identity):
    df = pd.DataFrame({"Locus": ["chr6:1-C-A"], "SpliceAI": [sai("PARK2", "0.40")], "Gene": [PRKN], "SYMBOL": ["PRKN"]})
    out = sp.parse_spliceai_custom(df, identity=identity)
    assert one(out, "spliceai_custom_match") == gi.ALIAS_RESOLVED
    assert one(out, "spliceai_custom_MAX") == pytest.approx(0.40)


def test_spliceai_routes_disagree_or_unknown_are_not_scored(identity):
    for symbol, status in (("DISAGREE", gi.AMBIGUOUS), ("NOSUCHGENE", gi.UNRESOLVED)):
        df = pd.DataFrame({"Locus": ["chr1:1-C-A"], "SpliceAI": [sai(symbol, "0.80")],
                           "Gene": ["ENSG00000000007"], "SYMBOL": ["OTHERNAME"]})
        out = sp.parse_spliceai_custom(df, identity=identity)
        assert np.isnan(one(out, "spliceai_custom_MAX")) and one(out, "spliceai_custom_match") == status
        assert one(out, "spliceai_custom_anygene_MAX") == pytest.approx(0.80)


def test_spliceai_allele_must_match_row_alt(identity):
    df = pd.DataFrame({"Locus": ["chr1:1-C-T"], "SpliceAI": [sai("PMS2", "0.70", allele="A")],
                       "Gene": [PMS2], "SYMBOL": ["PMS2"]})
    out = sp.parse_spliceai_custom(df, identity=identity)
    assert np.isnan(one(out, "spliceai_custom_MAX")) and one(out, "spliceai_custom_match") == gi.NO_ENTRY


def test_spliceai_single_gene_matches_previous_behaviour(identity):
    df = pd.DataFrame({"Locus": ["chr1:1-C-A"], "SpliceAI": ["A|PMS2|0.01|0.02|0.30|0.04|-5|6|7|8"],
                       "Gene": [PMS2], "SYMBOL": ["PMS2"]})
    out = sp.parse_spliceai_custom(df, identity=identity)
    assert [one(out, c) for c in ("spliceai_custom_DS_AG", "spliceai_custom_DS_AL", "spliceai_custom_DS_DG",
                                  "spliceai_custom_DS_DL")] == [0.01, 0.02, 0.30, 0.04]
    assert [one(out, c) for c in ("spliceai_custom_DP_AG", "spliceai_custom_DP_AL", "spliceai_custom_DP_DG",
                                  "spliceai_custom_DP_DL")] == [-5, 6, 7, 8]
    assert one(out, "spliceai_custom_MAX") == pytest.approx(0.30) and one(out, "spliceai_custom_SYMBOL") == "PMS2"


def test_spliceai_requires_identity():
    with pytest.raises(ValueError):
        sp.parse_spliceai_custom(pd.DataFrame({"SpliceAI": ["x"]}), identity=None)


def test_spliceai_fallback_to_vep_only_when_custom_absent():
    df = pd.DataFrame({"spliceai_custom_MAX": [0.0, np.nan, 0.5], "spliceAI_MAX": [0.3, 0.3, 0.1]})
    scores, methods = sp._extract_spliceai_details(df)
    assert list(scores) == [0.0, 0.3, 0.5]          # a genuine 0.00 custom score is not replaced
    assert methods.iloc[0].startswith("Custom") and methods.iloc[1].startswith("VEP")


# ------------------------------------------------------------------ Pangolin
def pang(gene, score):
    return f"{gene}|-10:{score}|20:-0.0|Warnings:"


def test_pangolin_two_genes_each_row_gets_its_own_score():
    value = pang(MADD + ".20", "0.62") + "&" + pang(MYBPC3 + ".12", "0.03")
    df = pd.DataFrame({"Pangolin": [value] * 3, "Gene": [MYBPC3, MADD, OTHER]})
    out = sp.parse_pangolin(df)
    assert out["Pangolin_max_score"].iloc[0] == pytest.approx(0.03)     # MYBPC3 row must not inherit MADD's 0.62
    assert out["Pangolin_max_score"].iloc[1] == pytest.approx(0.62)
    assert np.isnan(out["Pangolin_max_score"].iloc[2])
    assert list(out["Pangolin_status"]) == ["scored", "scored", "not_covered"]
    assert out["Pangolin_match"].iloc[2] == gi.NO_ENTRY
    assert np.allclose(out["Pangolin_anygene_max"], 0.62)


def test_pangolin_missing_and_no_row_gene():
    df = pd.DataFrame({"Pangolin": [None, pang(MADD, "0.2")], "Gene": [MADD, None]})
    out = sp.parse_pangolin(df)
    assert list(out["Pangolin_match"]) == [gi.NO_PREDICTION, gi.NO_ROW_GENE]
    assert out["Pangolin_max_score"].isna().all() and (out["Pangolin_status"] == "not_covered").all()


# ------------------------------------------------------------------ SPiP
def spip(nm, symbol, pred, allele="A"):
    f = [""] * 30
    f[0], f[2], f[3], f[4], f[11], f[12] = allele, "NTR", "00 %", f"{pred}", nm, symbol
    f[20], f[21], f[24], f[25] = "No", "0", "0", "None"
    return "|".join(f)


def test_spip_row_takes_its_own_transcript_not_the_first_entry():
    value = spip("NM_000009.1", "GENEZ", "0.90") + "&" + spip("NM_000008.3", "GENEZ", "0.10")
    df = pd.DataFrame({"Locus": ["chr1:1-C-A"] * 2, "SPiP": [value] * 2, "Feature": ["ENST00000000001"] * 2,
                       "MANE_SELECT": ["NM_000008.3", None], "MANE_PLUS_CLINICAL": [None, None], "Gene": [OTHER] * 2})
    ident = _bare_identity()
    out = sp.parse_spip(df, identity=ident)
    assert out["SPiP_prediction"].iloc[0] == pytest.approx(0.10)         # MANE transcript, not the first entry (0.90)
    assert out["SPiP_match"].iloc[0] == gi.MATCHED
    assert np.isnan(out["SPiP_prediction"].iloc[1]) and out["SPiP_match"].iloc[1] == gi.NOT_APPLICABLE_TRANSCRIPT
    assert out["SPiP_anygene_max_prediction"].iloc[1] == pytest.approx(0.90)


def _bare_identity(curated_nm=None):
    ident = gi.GeneIdentity.__new__(gi.GeneIdentity)
    ident._curated_nm = curated_nm or {}
    ident._hgnc_current, ident._hgnc_other = {}, {}
    return ident


def test_spip_curated_transcript_and_missing_nm_entry():
    value = spip("NM_000001.2", "GENEX", "0.50")
    ident = _bare_identity({"ENST00000000011": {"NM_000001"}})
    df = pd.DataFrame({"Locus": ["chr1:1-C-A"] * 2, "SPiP": [value] * 2,
                       "Feature": ["ENST00000000011", "ENST00000000012"],
                       "MANE_SELECT": [None, "NM_999999.1"], "MANE_PLUS_CLINICAL": [None, None]})
    out = sp.parse_spip(df, identity=ident)
    assert out["SPiP_prediction"].iloc[0] == pytest.approx(0.50) and out["SPiP_mechanism"].iloc[0] == "complex_splicing"
    assert np.isnan(out["SPiP_prediction"].iloc[1]) and out["SPiP_match"].iloc[1] == gi.NO_ENTRY


def test_spip_never_borrows_another_genes_entry_and_keeps_error_status():
    df = pd.DataFrame({"Locus": ["chr1:1-C-A"] * 2,
                       "SPiP": [spip("NM_000009.1", "NEIGHBOUR", "0.90"), "chr1:1 caused an error in SPiP execution"],
                       "Feature": ["ENST00000000001"] * 2, "MANE_SELECT": ["NM_000008.3"] * 2,
                       "MANE_PLUS_CLINICAL": [None, None]})
    out = sp.parse_spip(df, identity=_bare_identity())
    assert np.isnan(out["SPiP_prediction"].iloc[0]) and out["SPiP_match"].iloc[0] == gi.NO_ENTRY
    assert out["SPiP_status"].iloc[1] == "error" and out["SPiP_match"].iloc[1] == gi.PREDICTOR_ERROR


def test_spip_requires_identity():
    with pytest.raises(ValueError):
        sp.parse_spip(pd.DataFrame({"SPiP": ["x"]}), identity=None)


# ------------------------------------------------------------------ SPiP labelled same-gene columns
def test_spip_samegene_columns_label_other_transcript_result(identity):
    value = spip("NM_000009.1", "PMS2", "0.90") + "&" + spip("NM_000008.3", "PMS2", "0.10") + "&" + spip("NM_000050.1", "AIMP2", "0.60")
    df = pd.DataFrame({"Locus": ["chr7:1-C-A"] * 3, "SPiP": [value] * 3, "Gene": [PMS2] * 3, "SYMBOL": ["PMS2"] * 3,
                       "Feature": ["ENST00000000001", "ENST00000000002", "ENST00000000003"],
                       "MANE_SELECT": ["NM_000008.3", None, "NM_777777.1"], "MANE_PLUS_CLINICAL": [None] * 3})
    out = sp.parse_spip(df, identity=identity)
    # row 0: own transcript present -> strict columns filled; same-gene columns describe the GENE (max over PMS2's 2 NMs)
    assert out["SPiP_prediction"].iloc[0] == pytest.approx(0.10) and out["SPiP_match"].iloc[0] == gi.MATCHED
    assert (out["SPiP_samegene_prediction"] == 0.90).all() and (out["SPiP_samegene_transcript"] == "NM_000009").all()
    assert (out["SPiP_samegene_n_transcripts"] == 2).all()                       # AIMP2's entry is not counted
    # row 1: no RefSeq equivalent -> strict empty, but the labelled gene-level value is there
    assert np.isnan(out["SPiP_prediction"].iloc[1]) and out["SPiP_match"].iloc[1] == gi.NOT_APPLICABLE_TRANSCRIPT
    # row 2: has an NM that SPiP lacks -> strict empty (no_entry), gene-level still labelled
    assert out["SPiP_match"].iloc[2] == gi.NO_ENTRY and out["SPiP_samegene_prediction"].iloc[2] == 0.90


def test_spip_samegene_resolves_old_symbol_through_hgnc(identity):
    df = pd.DataFrame({"Locus": ["chr6:1-C-A"], "SPiP": [spip("NM_000100.1", "PARK2", "0.40")], "Gene": [PRKN],
                       "SYMBOL": ["PRKN"], "Feature": ["ENST00000000009"], "MANE_SELECT": [None], "MANE_PLUS_CLINICAL": [None]})
    out = sp.parse_spip(df, identity=identity)
    assert one(out, "SPiP_samegene_prediction") == pytest.approx(0.40) and one(out, "SPiP_match") == gi.NOT_APPLICABLE_TRANSCRIPT


# ------------------------------------------------------------------ dbNSFP per-transcript selection
from modules import predictors as pr  # noqa: E402


def dbn(feature, **cols):
    base = {"Ensembl_transcriptid": "ENST00000000001&ENST00000000002&ENST00000000003",
            "SIFT_score": "0.01&.&0.30", "AlphaMissense_score": "0.9&0.5&0.1", "REVEL_score": "0.7",
            "MetaRNN_score": "0.8&0.2&0.4"}
    base.update(cols)
    return {"Feature": feature, **base}


ALIGNED = ["Ensembl_transcriptid", "SIFT_score", "AlphaMissense_score", "MetaRNN_score"]


def test_dbnsfp_row_keeps_its_own_transcript_value():
    df = pd.DataFrame([dbn("ENST00000000003.5"), dbn("ENST00000000001"), dbn("ENST00000000002"), dbn("ENST00000999999")])
    out = pr.parse_dbnsfp_by_row_transcript(df, ALIGNED)
    assert list(out["AlphaMissense_score"].iloc[:3]) == ["0.1", "0.9", "0.5"]
    assert out["SIFT_score"].iloc[0] == "0.30" and np.isnan(out["SIFT_score"].iloc[2])      # '.' = missing for that transcript
    assert list(out["dbNSFP_match"]) == [pr.DBNSFP_MATCHED] * 3 + [pr.DBNSFP_NO_ENTRY]
    assert out["AlphaMissense_score"].iloc[3] is np.nan or pd.isna(out["AlphaMissense_score"].iloc[3])
    assert (out["REVEL_score"] == "0.7").all()                                                # single-valued column untouched
    assert out["Ensembl_transcriptid"].iloc[1] == "ENST00000000001"
    assert out["dbNSFP_transcripts_all"].iloc[0].count("&") == 2                              # original list kept for audit


def test_dbnsfp_anytranscript_max_keeps_the_old_variant_level_view():
    df = pd.DataFrame([dbn("ENST00000000002")])
    out = pr.parse_dbnsfp_by_row_transcript(df, ALIGNED)
    assert out["MetaRNN_score"].iloc[0] == "0.2" and out["MetaRNN_score_anytranscript_max"].iloc[0] == pytest.approx(0.8)


def test_dbnsfp_rows_without_dbnsfp_and_missing_resource():
    df = pd.DataFrame([dbn("ENST00000000001"), {"Feature": "ENST00000000009", "Ensembl_transcriptid": None,
                                                "SIFT_score": None, "AlphaMissense_score": None, "MetaRNN_score": None}])
    out = pr.parse_dbnsfp_by_row_transcript(df, ALIGNED)
    assert out["dbNSFP_match"].iloc[1] == pr.DBNSFP_NO_PREDICTION and pd.isna(out["SIFT_score"].iloc[1])
    with pytest.raises(ValueError):
        pr.parse_dbnsfp_by_row_transcript(df, None)


def test_dbnsfp_misaligned_list_is_not_used():
    df = pd.DataFrame([dbn("ENST00000000001", AlphaMissense_score="0.9&0.5")])           # 2 items for 3 transcripts
    out = pr.parse_dbnsfp_by_row_transcript(df, ALIGNED)
    assert pd.isna(out["AlphaMissense_score"].iloc[0]) and out["dbNSFP_match"].iloc[0] == pr.DBNSFP_MISALIGNED


def test_dbnsfp_single_transcript_list_not_matching_the_row_is_not_copied_to_it():
    """A one-transcript dbNSFP list must not leak onto rows of other transcripts (the COL6A5-style case)."""
    one_tx = {"Ensembl_transcriptid": "ENST00000265379", "SIFT_score": "0.0", "AlphaMissense_score": "0.8363",
              "MetaRNN_score": "0.9"}
    df = pd.DataFrame([{"Feature": "ENST00000312481", **one_tx}, {"Feature": "ENST00000265379", **one_tx}])
    out = pr.parse_dbnsfp_by_row_transcript(df, ALIGNED)
    assert out["dbNSFP_match"].iloc[0] == pr.DBNSFP_NO_ENTRY
    assert out[["Ensembl_transcriptid", "SIFT_score", "AlphaMissense_score", "MetaRNN_score"]].iloc[0].isna().all()
    assert out["AlphaMissense_score"].iloc[1] == "0.8363" and out["Ensembl_transcriptid"].iloc[1] == "ENST00000265379"


def test_dbnsfp_single_value_for_multi_transcript_list_is_kept_for_matched_rows_only():
    df = pd.DataFrame([dbn("ENST00000000002", SIFT_score="0.05"), dbn("ENST00000999999", SIFT_score="0.05")])
    out = pr.parse_dbnsfp_by_row_transcript(df, ALIGNED)
    assert out["SIFT_score"].iloc[0] == "0.05" and pd.isna(out["SIFT_score"].iloc[1])
