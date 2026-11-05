"""build_cDNA_MANE keeps HGVSc's missing values missing (MISC-33)."""
import io
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src/python/parsing")))
from modules.io_qc import build_cDNA_MANE  # noqa: E402

TSV = ("CHROM\tPOS\tREF\tALT\tHGVSc\tMANE_SELECT\tFeature\n"
       "chr1\t100\tA\tG\tENST1.2:c.100A>G\tNM_1.1\tENST1\n"
       "chr1\t200\tC\tT\t.\t.\t.\n")


def test_missing_hgvsc_stays_missing():
    data = pd.read_csv(io.StringIO(TSV), na_values=".", sep="\t", header=0, index_col=False, low_memory=False)
    out = build_cDNA_MANE(data)
    assert out["HGVSc"].iloc[0] == "ENST1.2:c.100A>G" and pd.isna(out["HGVSc"].iloc[1])
    assert out["CDNA_NAME"].iloc[0] == "NM_1.1:c.100A>G" and pd.isna(out["CDNA_NAME"].iloc[1])
