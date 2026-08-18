import sys
import os
import pytest
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python')))

import generate_clinical_prioritization_report as gcpr

def test_generate_clinical_prioritization_report(tmp_path):
    pq_path = tmp_path / "BAG3.parsed.clean.pq"
    df = pd.DataFrame({
        'Locus': ['chr10:119660000-119660001', 'chr10:119670000-119670001', 'chr10:119680000-119680001'],
        'SYMBOL': ['BAG3', 'BAG3', 'BAG3'],
        'HGVSc': ['ENST00000001:c.100G>A', 'ENST00000001:c.200+1G>A', 'ENST00000001:c.-10A>G'],
        'HGVSp': ['ENSP00000001:p.Arg34Gln', '-', '-'],
        'Consequence': ['missense_variant', 'splice_donor_variant', '5_prime_UTR_variant'],
        'PRIORITY_TIER': [
            'Tier 2 (Likely Deleterious / Strong Candidate)',
            'Tier 1 (Critical Pathogenic Candidate)',
            'Tier 3 (VUS / Moderate Potential)'
        ],
        'VARIANT_PRIORITY_SCORE': [62.5, 95.0, 35.0],
        'NEW_IMPACT': ['MODERATE', 'HIGH', 'MODERATE'],
        'IMPACT': ['MODERATE', 'HIGH', 'MODIFIER'],
        'gnomADv4_AF_grpmax_joint': [0.00001, 0.0, 0.0001],
        'spliceai_custom_MAX': [0.05, 0.95, 0.0],
        'spliceai_custom_DS_AG': [0.0, 0.0, 0.0],
        'spliceai_custom_DS_AL': [0.0, 0.0, 0.0],
        'spliceai_custom_DS_DG': [0.05, 0.95, 0.0],
        'spliceai_custom_DS_DL': [0.0, 0.0, 0.0],
        'am_pathogenicity': [0.85, 0.0, 0.0],
        'am_class': ['likely_pathogenic', 'likely_benign', 'likely_benign'],
        'REVEL_score': [0.78, 0.0, 0.0],
        'CADD_PHRED': [28.0, 35.0, 12.0],
        'CLNSIG': ['Uncertain_significance', 'Pathogenic', 'not_provided'],
        '5UTR_consequence': [None, None, 'uAUG_gained'],
        'pLI_gene_value': [0.98, 0.98, 0.98],
        'gene_priority': ['Tier_1', 'Tier_1', 'Tier_1'],
        'MANE_SELECT': ['ENST00000001.8', 'ENST00000001.8', 'ENST00000001.8']
    })
    df.to_parquet(pq_path)

    out_html = tmp_path / "BAG3_clinical_report.html"
    success = gcpr.generate_gene_report(str(pq_path), str(out_html))

    assert success is True
    assert out_html.exists()
    content = out_html.read_text()
    
    # Verify 4 tabs are present
    assert "General / Overview" in content
    assert "Splicing Alterations" in content
    assert "5'/3' UTR & Translation" in content
    assert "Missense Pathogenicity" in content
    
    # Verify Provenance & Disclaimers
    assert "Custom SpliceAI Active" in content
    assert "BAG3" in content
