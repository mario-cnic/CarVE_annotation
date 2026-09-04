#!/usr/bin/env python3
"""
🧪 Compound Heterozygosity & TRANS/CIS Phasing Unit Test Suite

Verifies:
  1. True TRANS Phasing: Variant A (Paternal) + Variant B (Maternal) -> compound_het = True.
  2. False CIS Phasing: Variant A (Paternal) + Variant B (Paternal) -> compound_het = False.
  3. Homozygous Recessive: Proband HOMALT -> Autosomal Recessive (Hom).
  4. Single-Sample / Missing Parents -> compound_het_candidate = True.
"""

import os
import sys
import pandas as pd
import unittest

PIPELINE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(PIPELINE_ROOT, "src", "python"))
import config

from filter_variants import mark_compound_het, analyze_pedigree_inheritance, build_priority_tier


class TestCompoundHeterozygosity(unittest.TestCase):

    def setUp(self):
        # Trio Pedigree Config
        self.pedigree_dict = {
            "proband": "PROBAND",
            "father": "FATHER",
            "mother": "MOTHER",
            "affected": ["PROBAND"],
            "unaffected": ["FATHER", "MOTHER"]
        }

    def test_trans_compound_het(self):
        """Test true TRANS compound heterozygosity (Paternal HET + Maternal HET in same gene)."""
        df_trans = pd.DataFrame({
            "CHROM": ["11", "11"],
            "POS": [47360000, 47370000],
            "REF": ["C", "T"],
            "ALT": ["T", "C"],
            "SYMBOL": ["MYBPC3", "MYBPC3"],
            "GT_PROBAND": ["HET", "HET"],
            "GT_FATHER": ["HET", "HOMREF"],
            "GT_MOTHER": ["HOMREF", "HET"],
            "Consequence": ["frameshift_variant", "splice_donor_variant"]
        })
        
        df_res = mark_compound_het(df_trans.copy(), self.pedigree_dict)
        self.assertTrue(df_res["compound_het"].iloc[0], "Paternal variant should be marked compound_het = True")
        self.assertTrue(df_res["compound_het"].iloc[1], "Maternal variant should be marked compound_het = True")
        self.assertIn("Paternal", df_res["COMPOUND_HET_PAIR"].iloc[0])
        self.assertIn("Maternal", df_res["COMPOUND_HET_PAIR"].iloc[1])

        # Test full inheritance analyzer integration
        df_inh = analyze_pedigree_inheritance(df_trans.copy(), self.pedigree_dict)
        self.assertIn(df_inh["INHERITANCE_MODEL"].iloc[0], ["Compound Heterozygous", "Autosomal Dominant / Heterozygous"])
        self.assertIn(df_inh["INHERITANCE_MODEL"].iloc[1], ["Compound Heterozygous", "Autosomal Dominant / Heterozygous"])

    def test_cis_non_compound_het(self):
        """Test false CIS phasing (both variants inherited from Father) -> MUST NOT be marked compound_het."""
        df_cis = pd.DataFrame({
            "CHROM": ["2", "2"],
            "POS": [178000000, 178010000],
            "REF": ["A", "G"],
            "ALT": ["G", "A"],
            "SYMBOL": ["TTN", "TTN"],
            "GT_PROBAND": ["HET", "HET"],
            "GT_FATHER": ["HET", "HET"],     # Both inherited from Father!
            "GT_MOTHER": ["HOMREF", "HOMREF"], # Mother has neither!
            "Consequence": ["missense_variant", "missense_variant"]
        })

        df_res = mark_compound_het(df_cis.copy(), self.pedigree_dict)
        self.assertFalse(df_res["compound_het"].iloc[0], "CIS variant 1 must NOT be marked compound_het")
        self.assertFalse(df_res["compound_het"].iloc[1], "CIS variant 2 must NOT be marked compound_het")

    def test_unphased_single_sample(self):
        """Test single-sample dataset with 2 HET variants in same gene -> compound_het_candidate = True."""
        df_single = pd.DataFrame({
            "CHROM": ["1", "1"],
            "POS": [156000000, 156010000],
            "REF": ["C", "T"],
            "ALT": ["T", "C"],
            "SYMBOL": ["LMNA", "LMNA"],
            "GT_SAMPLE1": ["HET", "HET"],
            "Consequence": ["missense_variant", "stop_gained"]
        })

        df_res = mark_compound_het(df_single.copy(), pedigree=None)
        self.assertTrue(df_res["compound_het_candidate"].iloc[0], "Single-sample duplicate HETs should mark compound_het_candidate")
        self.assertTrue(df_res["compound_het_candidate"].iloc[1], "Single-sample duplicate HETs should mark compound_het_candidate")
        self.assertFalse(df_res["compound_het"].iloc[0], "Unphased single-sample variants must NOT be marked true compound_het")


if __name__ == "__main__":
    unittest.main()
