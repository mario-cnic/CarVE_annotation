import os
import sys
import pytest
import pandas as pd
import numpy as np

# Add shared/utils/src to path
sys.path.insert(0, "/home/mruizp/data_lab_PGP/shared/utils/src")
from filter_variants import parse_branchpointer, NEW_BRANCHPOINT_COLUMNS

def test_parse_branchpointer_status():
    df = pd.DataFrame({
        "Branchpointer_prob": [0.95, np.nan, 0.40],
        "Branchpointer_U2_energy": [-6.2, np.nan, -4.5],
        "Branchpoint_disrupted": ["YES", "NO", "YES"],
        "LaBranchoR_score": [0.95, np.nan, np.nan],
        "LaBranchoR_acc_dist": [24, np.nan, np.nan],
    })
    
    res = parse_branchpointer(df.copy())
    
    # Verify all expected columns exist
    for col in NEW_BRANCHPOINT_COLUMNS.values():
        assert col in res.columns
        
    # Check status contract semantics
    assert res.loc[0, "Branchpoint_status"] == "scored"
    assert res.loc[0, "LaBranchoR_status"] == "scored"
    assert res.loc[0, "Branchpoint_disrupted"] == "YES"
    
    assert res.loc[1, "Branchpoint_status"] == "not_covered"
    assert res.loc[1, "LaBranchoR_status"] == "not_covered"
    assert res.loc[1, "Branchpoint_disrupted"] == "NO"
    
    assert res.loc[2, "Branchpoint_status"] == "scored"
    assert res.loc[2, "LaBranchoR_status"] == "not_covered"

def test_parse_branchpointer_missing_columns():
    df = pd.DataFrame({"Locus": ["chr1:100-A-G"]})
    res = parse_branchpointer(df.copy())
    assert "Locus" in res.columns
