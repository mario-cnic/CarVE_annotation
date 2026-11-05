"""Disease core panels hold real gene symbols (ACMG-9: `AGGL` was listed instead of `AGL`)."""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src/python"))
from parsing.modules import disease_hpo as dh  # noqa: E402

PANELS = [dh.DCM_CORE, dh.HCM_CORE, dh.ARVC_CORE, dh.CHANNELOPATHY_CORE, dh.AORTIC_CORE, dh.STORAGE_INFILTRATIVE_CORE]


def test_every_core_gene_is_a_gencode_symbol():
    # resources/spliceai_symbol_to_ensg.grch38.gencode45.tsv: GENCODE 45 symbols of SpliceAI's gene set
    table = REPO / "resources/spliceai_symbol_to_ensg.grch38.gencode45.tsv"
    symbols = {l.split("\t")[0] for l in table.read_text().splitlines()[1:]}
    unknown = sorted(g for panel in PANELS for g in panel if g not in symbols)
    assert unknown == []


def test_agl_is_a_storage_core_gene(tmp_path):
    gl = tmp_path / "genes.csv"
    gl.write_text("GENE,NAME,SOURCE,EVIDENCE\nAGL,amylo-alpha-1 6-glucosidase,Cardiomyopathies,Secondary\n")
    assert dh.load_cardiac_disease_gene_curation(str(gl))["AGL"]["evidence_tier"] == "Definitive / Primary Evidence"
    out = dh.apply_disease_phenotype_curation(pd.DataFrame({"SYMBOL": ["AGL"]}), str(gl), "storage")
    assert out["DISEASE_GENE_EVIDENCE"].iloc[0] == "Definitive / Primary Evidence"
