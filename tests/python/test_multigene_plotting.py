#!/usr/bin/env python3
"""
🧪 Multi-Gene / WES / WGS Plotting & Visualization Unit Test Suite

Verifies:
  1. Multi-Gene Dataset: create_transcript_visualization_figure() auto-switches to Chromosomal Manhattan Plot.
  2. Single-Gene Dataset: create_transcript_visualization_figure() renders single-gene transcript track.
  3. HTML Report Generation: generate_gene_report() runs cleanly on multi-gene datasets.
"""

import os
import sys
import tempfile
import pandas as pd
import unittest

PIPELINE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(PIPELINE_ROOT, "src", "python"))

from generate_clinical_prioritization_report import create_transcript_visualization_figure, generate_gene_report


class TestMultiGenePlotting(unittest.TestCase):

    def setUp(self):
        # Create Multi-Gene Synthetic Dataset (10 genes across 5 chromosomes)
        self.df_multigene = pd.DataFrame({
            "CHROM": ["chr1", "chr1", "chr2", "chr11", "chr11", "chr17", "chrX"],
            "POS": [100000, 200000, 500000, 47328406, 47350000, 7577000, 1500000],
            "Locus": ["1:100000-A-G", "1:200000-C-T", "2:500000-T-C", "11:47328406-C-A", "11:47350000-G-A", "17:7577000-C-T", "X:1500000-A-T"],
            "SYMBOL": ["LMNA", "LMNA", "TTN", "MYBPC3", "MYBPC3", "TP53", "DMD"],
            "Consequence": ["missense_variant", "synonymous_variant", "stop_gained", "missense_variant", "stop_gained", "missense_variant", "frameshift_variant"],
            "HGVSc": ["c.100A>G", "c.200C>T", "c.500T>C", "c.526C>A", "c.1200G>A", "c.750C>T", "c.1500delA"],
            "HGVSp": ["p.Lys34Glu", "p.Pro67Pro", "p.Trp167*", "p.Arg176Ser", "p.Trp400*", "p.Arg251Trp", "p.Lys501fs"],
            "VARIANT_PRIORITY_SCORE": [45.0, 10.0, 85.0, 20.0, 50.0, 70.0, 90.0],
            "PRIORITY_TIER": [
                "Tier 2 (Likely Deleterious / Strong Candidate)",
                "Tier 4 (Benign / Tolerated)",
                "Tier 1 (Critical Pathogenic Candidate)",
                "Tier 4 (Benign / Tolerated)",
                "Tier 2 (Likely Deleterious / Strong Candidate)",
                "Tier 1 (Critical Pathogenic Candidate)",
                "Tier 1 (Critical Pathogenic Candidate)"
            ],
            "am_pathogenicity": [0.8, 0.1, 0.95, 0.3, 0.85, 0.9, 0.99],
            "SPLICE_MAX_UNIFIED": [0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.9]
        })

    def test_multigene_figure_switch(self):
        """Multi-gene dataset should render Chromosomal Manhattan Plot."""
        fig = create_transcript_visualization_figure(self.df_multigene)
        self.assertIsNotNone(fig, "Multi-gene figure should not be None")
        self.assertIn("Chromosomal Manhattan Distribution", fig.layout.title.text, "Figure title should indicate Manhattan distribution")

    def test_singlegene_figure_rendering(self):
        """Single-gene dataset should render single-gene transcript map."""
        df_single = self.df_multigene[self.df_multigene["SYMBOL"] == "MYBPC3"].copy()
        fig = create_transcript_visualization_figure(df_single)
        self.assertIsNotNone(fig, "Single-gene figure should not be None")

    def test_multigene_report_generation(self):
        """Multi-gene dataset should generate full 5-tab HTML report cleanly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_pq = os.path.join(tmpdir, "multigene.pq")
            tmp_html = os.path.join(tmpdir, "multigene_report.html")
            
            self.df_multigene.to_parquet(tmp_pq)
            success = generate_gene_report(tmp_pq, tmp_html)
            self.assertTrue(success, "Report generation should return True")
            self.assertTrue(os.path.exists(tmp_html), "Multi-gene HTML report should be created")


if __name__ == "__main__":
    unittest.main()
