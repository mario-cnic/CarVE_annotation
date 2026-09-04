#!/usr/bin/env python3
"""
Unit Test Suite for Multi-Gene Cohort Master Dashboard Generator
"""

import os
import sys
import unittest
import pandas as pd
import tempfile

PIPELINE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(PIPELINE_ROOT, "src", "python"))
import config
from generate_cohort_master_dashboard import generate_dashboard_html

class TestCohortMasterDashboard(unittest.TestCase):
    def setUp(self):
        self.sample_df = pd.DataFrame({
            "Locus": ["11:47353297:G:A", "2:178000000:A:G"],
            "SYMBOL": ["MYBPC3", "TTN"],
            "Consequence": ["stop_gained", "missense_variant"],
            "HGVSc": ["c.526C>T", "c.1504C>T"],
            "HGVSp": ["p.Arg176Ter", "p.Pro502Ser"],
            "PRIORITY_TIER": ["Class 1 (High Priority)", "Class 2 (Moderate Priority)"],
            "VARIANT_PRIORITY_SCORE": [95.0, 75.0],
            "gnomADv4_AF_grpmax_joint": [0.00001, 0.0005],
            "REVEL_score": [0.92, 0.65],
            "am_pathogenicity": [0.98, 0.45],
            "SPiP_prediction": ["Alteration", "No_alteration"],
            "spliceai_custom_MAX": [0.85, 0.05]
        })

    def test_generate_dashboard_html(self):
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as tmp:
            tmp_path = tmp.name

        try:
            generate_dashboard_html(self.sample_df, "test_cohort_run", tmp_path)
            self.assertTrue(os.path.exists(tmp_path))
            self.assertGreater(os.path.getsize(tmp_path), 1000)

            with open(tmp_path, "r", encoding="utf-8") as f:
                html_txt = f.read()
                self.assertIn("MYBPC3", html_txt)
                self.assertIn("TTN", html_txt)
                self.assertIn("Multi-Gene Cohort Master Variant Dashboard", html_txt)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

if __name__ == "__main__":
    unittest.main()
