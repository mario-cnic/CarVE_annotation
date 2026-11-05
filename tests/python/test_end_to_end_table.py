"""VCF -> VCF_TO_TSV (vcf_parser_pysam.py) -> TSV_TO_TABLE (filter_variants.py, whole and chunked) ->
TAG_TRANSCRIPT_PRIORITY, with the repo's identity resources. Synthetic records; public gene IDs as labels.

Checks that every record reaches the table, chunked output equals whole-file output, and a few per-row
contracts (gene/transcript attribution, missing values, intron offset)."""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("vcf_parser")
REPO = Path(__file__).resolve().parents[2]
RES = REPO / "resources"
MYH7, TTN = "ENSG00000092054", "ENSG00000155657"
COLS = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "HGVSc", "HGVSp",
        "CDS_position", "STRAND", "VARIANT_CLASS", "MANE_SELECT", "MANE_PLUS_CLINICAL", "NMD", "am_pathogenicity",
        "am_class", "MaxEntScan_ref", "MaxEntScan_alt", "ClinVarVCF_CLNSIG", "ClinVarVCF_CLNREVSTAT", "gnomADv4_AF_joint",
        "gnomADv4_AF_grpmax_joint", "gnomADv4_AN_joint", "Ensembl_transcriptid", "REVEL_score", "MetaRNN_score",
        "SpliceVault_top_events", "SpliceVault_site_pos", "ReMM", "AVI_PHRED", "AVI_RAW", "5UTR_annotation", "5UTR_consequence"]


def csq(**kw):
    return "|".join(str(kw.get(c, "")) for c in COLS)


RECORDS = [  # chrom, pos, ref, alt, CSQ entries, other INFO
    ("chr1", 1000000, "G", "A", [], ""),                                                    # no CSQ
    ("chr1", 2000000, "C", "G", [csq(Allele="G", Consequence="intergenic_variant", IMPACT="MODIFIER", VARIANT_CLASS="SNV")], ""),
    ("chr2", 178527000, "AT", "A", [csq(Allele="-", Consequence="frameshift_variant", IMPACT="HIGH", SYMBOL="TTN", Gene=TTN,
                                        Feature_type="Transcript", Feature="ENST00000589042", HGVSc="ENST00000589042.5:c.5000del",
                                        CDS_position="5000", STRAND="-1", VARIANT_CLASS="deletion", MANE_SELECT="NM_001267550.2",
                                        gnomADv4_AF_joint="0.00001", gnomADv4_AF_grpmax_joint="0.00002", gnomADv4_AN_joint="150000")], ""),
    ("chr14", 23412800, "C", "T", [
        csq(Allele="T", Consequence="missense_variant", IMPACT="MODERATE", SYMBOL="MYH7", Gene=MYH7, Feature_type="Transcript",
            Feature="ENST00000355349", HGVSc="ENST00000355349.4:c.1000G>A", CDS_position="1000", STRAND="-1", VARIANT_CLASS="SNV",
            MANE_SELECT="NM_000257.4", am_pathogenicity="0.95", am_class="likely_pathogenic",
            Ensembl_transcriptid="ENST00000355349&ENST00000999999", REVEL_score="0.91&0.20", MetaRNN_score="0.8&0.1",
            ClinVarVCF_CLNSIG="Pathogenic", ClinVarVCF_CLNREVSTAT="reviewed_by_expert_panel", AVI_PHRED="12.5", AVI_RAW="0.3", ReMM="0.4"),
        csq(Allele="T", Consequence="intron_variant", IMPACT="MODIFIER", SYMBOL="MYH7", Gene=MYH7, Feature_type="Transcript",
            Feature="ENST00000999999", HGVSc="ENST00000999999.1:c.-45+3G>A", STRAND="-1", VARIANT_CLASS="SNV")],
     f"SpliceAI=T|MYH7|0.01|0.00|0.02|0.00|1|2|3|4;Pangolin={MYH7}.1|5:0.3|-3:-0.5|Warnings:"),
]


@pytest.fixture(scope="module")
def tables(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("e2e")
    header = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={c},length=250000000>" for c in ("chr1", "chr2", "chr14")]
    header += ['##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: '
               + "|".join(COLS) + '">', '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="x">',
               '##INFO=<ID=Pangolin,Number=.,Type=String,Description="x">', "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    lines = []
    for chrom, pos, ref, alt, entries, info in RECORDS:
        fields = ([f"CSQ={','.join(entries)}"] if entries else []) + ([info] if info else [])
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t{';'.join(fields) or '.'}")
    (tmp / "in.vcf").write_text("\n".join(header + lines) + "\n")
    (tmp / "cols.txt").write_text("\n".join(COLS) + "\n")
    (tmp / "hgnc.tsv").write_text(f"symbol\tensembl_gene_id\tprev_symbol\talias_symbol\nMYH7\t{MYH7}\t\t\nTTN\t{TTN}\t\t\n")
    py = sys.executable

    def run(*args):
        res = subprocess.run([py, *map(str, args)], cwd=tmp, capture_output=True, text=True)
        assert res.returncode == 0, res.stderr[-3000:]

    run(REPO / "src/python/parsing/vcf_parser_pysam.py", "--input", "in.vcf", "--output", "scratch.tsv",
        "--vep_columns", "cols.txt", "--add_info", "--add_vep", "--overwrite")
    common = ["--input", "scratch.tsv", "--hgnc-table", "hgnc.tsv",
              "--disease-gene-list", REPO / "tests/fixtures/mini_disease_gene_list.csv",
              "--spliceai-symbol-map", RES / "spliceai_symbol_to_ensg.grch38.gencode45.tsv",
              "--enst-spip-map", RES / "enst_to_spip_nm.grch38.gencode45.tsv",
              "--dbnsfp-aligned-columns", RES / "dbnsfp_transcript_aligned_columns.txt",
              "--gene-constraint-table", RES / "gnomad_v4.1.2_gene_constraint.grch38.tsv", "--skip_quality_filter"]
    out = {}
    for name, extra in (("whole", []), ("chunk1", ["--chunk-rows", "1"]), ("chunk2", ["--chunk-rows", "2"])):
        run(REPO / "src/python/parsing/filter_variants.py", *common, "--output", f"{name}.pq", *extra)
        out[name] = pd.read_parquet(tmp / f"{name}.pq")
    run(REPO / "src/python/tag_transcript_priority.py", "--input", "chunk2.pq", "--output", "final.pq",
        "--gene-transcript-mapping", RES / "gene_transcript_mapping.txt")
    out["final"] = pd.read_parquet(tmp / "final.pq")
    return out


def test_every_record_reaches_the_table(tables):
    expected = {f"{c}:{p}-{r}-{a}" for c, p, r, a, _, _ in RECORDS}
    assert set(tables["whole"]["Locus"]) == expected and len(tables["final"]) == 5


@pytest.mark.parametrize("name", ["chunk1", "chunk2"])
def test_chunked_equals_whole(tables, name):
    pd.testing.assert_frame_equal(tables[name], tables["whole"])


def test_row_contracts(tables):
    t = tables["final"].set_index("Feature", drop=False)
    mane, other = t.loc["ENST00000355349"], t.loc["ENST00000999999"]
    assert mane["REVEL_score"] == 0.91 and pd.isna(other["REVEL_score"])           # own transcript's dbNSFP value
    assert mane["Pangolin_max_score"] == 0.5 and mane["Pangolin_match"] == "matched"
    assert (other["intron_offset_signed"], other["splice_side"]) == (3, "donor")    # c.-45+3
    assert mane["TRANSCRIPT_PRIORITY_TIER"] in (1, 2)
    no_csq = tables["final"][tables["final"]["Locus"] == "chr1:1000000-G-A"].iloc[0]
    assert pd.isna(no_csq["HGVSc"]) and pd.isna(no_csq["Feature"])
    assert tables["whole"]["NEW_IMPACT"].cat.categories.tolist() == ["HIGH", "MODERATE", "LOW", "MODIFIER"]
