"""
ClinVar significance masks from the source-selected classification.
"""
import pandas as pd
from .clinvar import annotate_clinvar, clinvar_masks


def _extract_clinvar_status(df: pd.DataFrame):
    """Extract source-selected ClinVar masks and expose selection provenance."""
    classification, _, _, _ = annotate_clinvar(df)
    return clinvar_masks(classification)
