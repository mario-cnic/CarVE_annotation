"""Seeded random records through VCF_TO_TSV and TSV_TO_TABLE: chunked output equals whole-file output,
one row per CSQ entry (one for a record without CSQ), dbNSFP values only from the row's own transcript."""
import random
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("vcf_parser")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_end_to_end_table import COLS, MYH7, REPO, RES, TTN, csq  # noqa: E402


def random_records(seed):
    rng = random.Random(seed)
    genes = [("MYH7", MYH7, "ENST00000355349", "NM_000257.4"), ("TTN", TTN, "ENST00000589042", "NM_001267550.2")]
    def maybe(v, p=0.5): return v if rng.random() < p else ""
    recs = []
    for i in range(60):
        chrom = rng.choice(["chr1", "chr2", "chr14"]); pos = 1000 + i * 37
        ref = rng.choice("ACGT"); alt = rng.choice([b for b in "ACGT" if b != ref])
        if rng.random() < 0.15: ref, alt = ref + "T", ref
        entries = []
        for _ in range(rng.choice([0, 1, 1, 2, 3])):
            sym, gid, tx, nm = rng.choice(genes)
            entries.append(csq(Allele=alt if len(alt) == 1 and len(ref) == 1 else "-", Consequence=rng.choice(["missense_variant", "intron_variant", "synonymous_variant", "stop_gained", "splice_region_variant&intron_variant", "5_prime_UTR_variant"]),
                IMPACT=rng.choice(["HIGH", "MODERATE", "LOW", "MODIFIER"]), SYMBOL=sym, Gene=gid, Feature_type="Transcript", Feature=tx if rng.random() < 0.6 else f"ENST0000099{i:04d}",
                HGVSc=maybe(f"{tx}.1:c.{rng.randint(1,900)}{rng.choice(['', '+5', '-3'])}{ref}>{alt}", 0.7), CDS_position=maybe(str(rng.randint(1, 900))), STRAND="-1",
                VARIANT_CLASS="SNV" if len(ref) == len(alt) else "deletion", MANE_SELECT=maybe(nm), am_pathogenicity=maybe(f"{rng.random():.3f}"), am_class=maybe(rng.choice(["likely_benign", "ambiguous", "likely_pathogenic"])),
                MaxEntScan_ref=maybe(f"{rng.uniform(-2,12):.2f}"), MaxEntScan_alt=maybe(f"{rng.uniform(-2,12):.2f}"),
                ClinVarVCF_CLNSIG=maybe(rng.choice(["Pathogenic", "Benign", "Uncertain_significance", "Conflicting_classifications_of_pathogenicity"]), 0.3),
                ClinVarVCF_CLNREVSTAT=maybe("criteria_provided&_single_submitter", 0.3), gnomADv4_AF_joint=maybe(f"{rng.random()/100:.6f}"), gnomADv4_AF_grpmax_joint=maybe(f"{rng.random()/50:.6f}"),
                gnomADv4_AN_joint=maybe("150000"), Ensembl_transcriptid=maybe(f"{tx}&ENST00000999999"), REVEL_score=maybe(f"{rng.random():.3f}&{rng.random():.3f}"), MetaRNN_score=maybe(f"{rng.random():.3f}"),
                SpliceVault_top_events=maybe(f"Top1:CD:{rng.randint(-50,50)}:{rng.randint(1,9)}:inframe", 0.2), SpliceVault_site_pos=maybe(str(pos + 10), 0.2),
                ReMM=maybe(f"{rng.random():.3f}"), AVI_PHRED=maybe(f"{rng.uniform(0,30):.2f}"), AVI_RAW=maybe(f"{rng.random():.3f}"),
                **{"5UTR_consequence": maybe("5_prime_UTR_premature_start_codon_gain_variant", 0.1)}))
        info = []
        if rng.random() < 0.5: info.append(f"SpliceAI={alt}|MYH7|{rng.random():.2f}|0.00|{rng.random():.2f}|0.00|1|2|3|4")
        if rng.random() < 0.4: info.append(f"Pangolin={MYH7}.1|5:{rng.random():.2f}|-3:-{rng.random():.2f}|Warnings:")
        recs.append((chrom, pos, ref, alt, entries, ";".join(info)))
    recs.sort(key=lambda r: (["chr1", "chr2", "chr14"].index(r[0]), r[1]))
    return recs


def write_vcf(path, recs):
    header = ["##fileformat=VCFv4.2"] + [f"##contig=<ID={c},length=250000000>" for c in ("chr1", "chr2", "chr14")]
    header += ['##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: '
               + "|".join(COLS) + '">', '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="x">',
               '##INFO=<ID=Pangolin,Number=.,Type=String,Description="x">', "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    lines = [f"{c}\t{p}\t.\t{r}\t{a}\t60\tPASS\t{';'.join(([f'CSQ={chr(44).join(e)}'] if e else []) + ([i] if i else [])) or '.'}"
             for c, p, r, a, e, i in recs]
    path.write_text("\n".join(header + lines) + "\n")


@pytest.mark.parametrize("seed", [7])
def test_random_records(tmp_path, seed):
    recs = random_records(seed)
    write_vcf(tmp_path / "in.vcf", recs)
    (tmp_path / "cols.txt").write_text("\n".join(COLS) + "\n")
    (tmp_path / "hgnc.tsv").write_text(f"symbol\tensembl_gene_id\tprev_symbol\talias_symbol\nMYH7\t{MYH7}\t\t\nTTN\t{TTN}\t\t\n")

    def run(*args):
        res = subprocess.run([sys.executable, *map(str, args)], cwd=tmp_path, capture_output=True, text=True)
        assert res.returncode == 0, res.stderr[-3000:]

    run(REPO / "src/python/parsing/vcf_parser_pysam.py", "--input", "in.vcf", "--output", "scratch.tsv",
        "--vep_columns", "cols.txt", "--add_info", "--add_vep", "--overwrite")
    common = ["--input", "scratch.tsv", "--hgnc-table", "hgnc.tsv",
              "--disease-gene-list", REPO / "tests/fixtures/mini_disease_gene_list.csv",
              "--spliceai-symbol-map", RES / "spliceai_symbol_to_ensg.grch38.gencode45.tsv",
              "--enst-spip-map", RES / "enst_to_spip_nm.grch38.gencode45.tsv",
              "--dbnsfp-aligned-columns", RES / "dbnsfp_transcript_aligned_columns.txt",
              "--gene-constraint-table", RES / "gnomad_v4.1.2_gene_constraint.grch38.tsv", "--skip_quality_filter"]
    tables = {}
    for name, extra in (("whole", []), ("c1", ["--chunk-rows", "1"]), ("c7", ["--chunk-rows", "7"])):
        run(REPO / "src/python/parsing/filter_variants.py", *common, "--output", f"{name}.pq", *extra)
        tables[name] = pd.read_parquet(tmp_path / f"{name}.pq")
    whole = tables["whole"]
    for name in ("c1", "c7"):
        pd.testing.assert_frame_equal(tables[name], whole)
    assert len(whole) == sum(len(e) or 1 for *_, e, _ in recs)
    matched = whole[whole["dbNSFP_match"] == "matched_transcript"]
    own = matched.apply(lambda r: str(r["Feature"]).split(".")[0] in str(r["dbNSFP_transcripts_all"]).split("&"), axis=1)
    assert own.all()
