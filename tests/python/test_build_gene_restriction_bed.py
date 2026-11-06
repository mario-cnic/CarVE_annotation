"""Gene-restriction BED: symbol / case / Ensembl-ID matching and primary-contig fallback (MISC-4)."""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "src/python/build_gene_restriction_bed.py"
sys.path.insert(0, str(REPO / "src/python"))
import build_gene_restriction_bed as gb  # noqa: E402


def gtf(path, genes):
    """genes: (chrom, start, end, strand, gene_id, gene_name)."""
    path.write_text("".join(
        f'{c}\tsrc\tgene\t{s}\t{e}\t.\t{st}\t.\tgene_id "{gid}"; gene_name "{n}";\n' for c, s, e, st, gid, n in genes))
    return path


MANE = [("chr1", 1000, 2000, "+", "ENSG00000000001.5", "GENEA"),
        ("chr10", 5000, 6000, "-", "ENSG00000000002.3", "C10orf71"),
        ("chrX", 9000, 9500, "+", "ENSG00000000003.1", "TAFAZZIN"),
        ("chr11_KN196481v1_fix", 35351, 52488, "-", "ENSG00000000004.2", "SLC37A4")]
PRIMARY = [("chr11", 119024000, 119031000, "-", "ENSG00000000004.9", "SLC37A4"),
           ("chrM", 3307, 4262, "+", "ENSG00000198888.2", "MT-ND1"),
           ("chr11_KN196481v1_fix", 1, 9, "-", "ENSG00000000005.1", "ONLYALT")]


def test_resolution(tmp_path):
    mane = gb.read_gtf_genes(str(gtf(tmp_path / "mane.gtf", MANE)))
    primary = gb.read_gtf_genes(str(gtf(tmp_path / "primary.gtf", PRIMARY)))
    listed = ["GENEA", "C10ORF71", "ENSG00000000003", "SLC37A4", "MT-ND1", "ONLYALT", "NOSUCH"]
    found, missing = gb.resolve_genes(listed, mane, primary)
    how = {g: [r[4] for r in recs] for g, recs in found.items()}
    assert how["GENEA"] == ["mane_symbol"]
    assert how["C10ORF71"] == ["mane_symbol_case"]
    assert how["ENSG00000000003"] == ["mane_ensembl_id"]
    assert how["SLC37A4"] == ["mane_symbol", "primary_gtf_symbol"]          # alt record kept, primary added
    assert "MT-ND1" not in how                                               # no introns: left out
    assert missing == {"MT-ND1": "mitochondrial", "ONLYALT": "no_match",     # primary GTF record on alt: not used
                       "NOSUCH": "no_match"}
    found, missing = gb.resolve_genes(["SLC37A4"], mane)                      # without the primary GTF
    assert missing == {"SLC37A4": "mane_alt_contig_only"} and len(found["SLC37A4"]) == 1


def run(tmp_path, *extra):
    out = tmp_path / "out.bed"
    subprocess.run([sys.executable, str(SCRIPT), "--gene-list", str(tmp_path / "genes.csv"), "--mane-gtf",
                    str(tmp_path / "mane.gtf"), "--output", str(out), "--margin", "100", *extra],
                   check=True, capture_output=True, text=True)
    return [l.split("\t") for l in out.read_text().splitlines()]


def test_exact_matches_give_the_previous_bed(tmp_path):
    (tmp_path / "genes.csv").write_text("GENE,NAME\nGENEA,a\nC10orf71,b\nTAFAZZIN,c\n")
    gtf(tmp_path / "mane.gtf", MANE)
    rows = run(tmp_path)
    # what the previous script wrote for the same input: both contig spellings, 0-based start - margin
    assert rows == [["1", "899", "2100", "GENEA", ".", "+"], ["10", "4899", "6100", "C10orf71", ".", "-"],
                    ["X", "8899", "9600", "TAFAZZIN", ".", "+"], ["chr1", "899", "2100", "GENEA", ".", "+"],
                    ["chr10", "4899", "6100", "C10orf71", ".", "-"], ["chrX", "8899", "9600", "TAFAZZIN", ".", "+"]]


def test_primary_gtf_adds_the_primary_region(tmp_path):
    (tmp_path / "genes.csv").write_text("GENE,NAME\nSLC37A4,a\n")
    gtf(tmp_path / "mane.gtf", MANE)
    gtf(tmp_path / "primary.gtf", PRIMARY)
    unmatched = tmp_path / "unmatched.tsv"
    rows = run(tmp_path, "--primary-gtf", str(tmp_path / "primary.gtf"), "--unmatched-output", str(unmatched))
    assert {r[0] for r in rows} == {"11", "chr11", "11_KN196481v1_fix", "chr11_KN196481v1_fix"}
    assert unmatched.read_text() == ""


HGNC = [  # symbol, status, alias_symbol, prev_symbol, ensembl_gene_id
    ("DRC4", "Approved", "", "GAS8", "ENSG00000000010"),
    ("TAFAZZIN", "Approved", "BTHS", "TAZ", "ENSG00000000003"),
    ("WWTR1", "Approved", "TAZ|SHARED", "", "ENSG00000000020"),
    ("OTHER", "Approved", "SHARED|FRA16E", "", "ENSG00000000030"),
    ("FRA16E", "Approved", "", "", ""),
    ("OLDGENE", "Entry Withdrawn", "", "", "ENSG00000000040"),
]


def test_renamed_genes_through_hgnc(tmp_path):
    table = tmp_path / "hgnc.txt"
    table.write_text("symbol\tstatus\talias_symbol\tprev_symbol\tensembl_gene_id\n" + "".join(
        f'{s}\t{st}\t"{a}"\t"{p}"\t{e}\n' for s, st, a, p, e in HGNC))
    hgnc = gb.read_hgnc(str(table))
    # approved > previous > alias; shared names and an approved symbol without an Ensembl ID map nowhere
    assert hgnc["DRC4"] == hgnc["GAS8"] == "ENSG00000000010" and hgnc["TAZ"] == "ENSG00000000003"
    assert "SHARED" not in hgnc and "FRA16E" not in hgnc and "OLDGENE" not in hgnc
    mane = gb.read_gtf_genes(str(gtf(tmp_path / "mane.gtf", MANE)))
    primary = gb.read_gtf_genes(str(gtf(tmp_path / "primary.gtf",
                                        PRIMARY + [("chr17", 5000, 6000, "+", "ENSG00000000010.4", "GAS8")])))
    found, missing = gb.resolve_genes(["DRC4", "TAZ", "SHARED", "FRA16E"], mane, primary, hgnc)
    assert {g: [r[4] for r in recs] for g, recs in found.items()} == {"DRC4": ["primary_gtf_hgnc"], "TAZ": ["mane_hgnc"]}
    assert found["DRC4"][0][:2] == ("chr17", 5000)
    assert missing == {"SHARED": "no_match", "FRA16E": "no_match"}
