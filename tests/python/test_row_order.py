"""Row order of the parsed table: kept in the input's coordinate order, checked, not re-sorted."""
import importlib.util
import os

import pandas as pd

# Loaded by path under its own name: other tests put a different `modules` package (shared/utils) on sys.path.
IO_QC = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/python/parsing/modules/io_qc.py'))
spec = importlib.util.spec_from_file_location("vendored_io_qc", IO_QC)
io_qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(io_qc)
check_coordinate_order = io_qc.check_coordinate_order


def frame(rows):
    return pd.DataFrame(rows, columns=["CHROM", "POS"])


def test_sorted_blocks_are_ordered_whatever_the_contig_order():
    df = frame([("chr2", 5), ("chr2", 5), ("chr2", 9), ("chr1", 3), ("chrX", 1), ("chrX", 7)])
    out = check_coordinate_order(df)
    assert out["ordered"] and out["contigs"] == 3 and out["contig_blocks"] == 3


def test_a_contig_that_reappears_is_not_ordered():
    out = check_coordinate_order(frame([("chr1", 1), ("chr2", 1), ("chr1", 9)]))
    assert out["ordered"] is False and out["contig_blocks"] == 3 and out["contigs"] == 2


def test_a_position_decrease_inside_a_contig_is_not_ordered():
    out = check_coordinate_order(frame([("chr1", 10), ("chr1", 5), ("chr2", 1)]))
    assert out["ordered"] is False and out["position_decreases"] == 1


def test_position_restart_at_a_new_contig_is_fine_and_ties_are_allowed():
    assert check_coordinate_order(frame([("1", 900), ("1", 900), ("2", 1)]))["ordered"]


def test_numeric_not_lexicographic_positions():
    assert check_coordinate_order(frame([("chr1", 99), ("chr1", 100), ("chr1", 1000)]))["ordered"]


def test_missing_columns_or_no_rows_are_not_checked():
    assert check_coordinate_order(pd.DataFrame({"Locus": ["chr1:1-A-C"]}))["checked"] is False
    assert check_coordinate_order(frame([]))["checked"] is False
