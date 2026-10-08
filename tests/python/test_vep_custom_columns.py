"""Static checks tying the VEP --custom annotations to the parser column whitelist."""
import re
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
VEP_NF = (REPO / "modules/local/vep.nf").read_text()
COLS = (REPO / "resources/all_but_old_gnomad_vep_cols.txt").read_text().split()


def custom_columns(short_name):
    """CSQ column names VEP emits for one --custom source: the short name plus short_name_<field>."""
    m = re.search(rf"short_name={short_name},.*?fields=(\S+)", VEP_NF)
    assert m, f"no --custom with short_name={short_name} in vep.nf"
    fields = [f for f in m.group(1).replace("\\", "").split("%") if f]
    return [short_name] + [f"{short_name}_{f}" for f in fields]


def test_clinvar_columns_are_whitelisted():
    missing = [c for c in custom_columns("ClinVarVCF") if c not in COLS]
    assert not missing, f"ClinVar columns not in the parser whitelist: {missing}"


def test_clinvar_vcf_param_is_used_and_declared():
    assert "file=${params.clinvar_vcf}" in VEP_NF
    assert re.search(r"^\s*clinvar_vcf\s*=", (REPO / "nextflow.config").read_text(), re.M)
    assert "clinvar_vcf" in (REPO / "run_annotate_vcf.sh").read_text()


def test_clinvar_columns_do_not_collide_with_other_columns_by_case():
    counts = Counter(x.lower() for x in COLS)
    clash = [c for c in custom_columns("ClinVarVCF") if counts[c.lower()] > 1]
    assert not clash, f"ClinVar columns that differ from another column only by case: {clash}"


def test_gnomad_requests_are_unique_and_keep_homozygote_counts():
    columns = custom_columns("gnomADv4")
    assert len(columns) == len(set(columns))
    assert all(c in COLS for c in columns)
    assert not any("age_hist_hom" in c for c in columns + COLS)
    for field in ("joint", "joint_nfe", "joint_afr", "grpmax_joint", "exomes", "genomes"):
        assert f"gnomADv4_nhomalt_{field}" in columns
    bash = (REPO / "src/hpc/annotate_vep_vars.sh").read_text()
    fields = re.search(r"short_name=gnomADv4,.*?fields=(\S+)", bash).group(1).split("%")
    assert columns[1:] == [f"gnomADv4_{f}" for f in fields]


def test_cached_gnomad_is_disabled_and_selector_is_unique():
    assert len(COLS) == len(set(COLS))
    assert not any(c.startswith(("gnomADe_", "gnomADg_")) or c in ("MAX_AF", "MAX_AF_POPS") for c in COLS)
    for source in (VEP_NF, (REPO / "src/hpc/annotate_vep_vars.sh").read_text()):
        flags = set(re.findall(r"--([a-z_0-9]+)", source.split("--no_stats")[0]))
        assert not flags & {"everything", "af_gnomad", "af_gnomade", "af_gnomadg", "max_af"}
        assert {"check_existing", "regulatory", "mane", "af", "af_1kg"} <= flags
        assert "--plugin Mastermind," in source
        assert "--plugin LoFtool," in source


def test_dbnsfp_deduplication_preserves_scores_and_population_fields():
    fields = [x for x in (REPO / "resources/vep_dbnsfp_fields.txt").read_text().splitlines() if x and not x.startswith("#")]
    assert len(fields) == len(set(fields)) == 454
    assert not {"TSL", "APPRIS"} & set(fields)
    assert {"Ensembl_transcriptid", "REVEL_score", "MetaRNN_score", "gnomAD_exomes_AF", "gnomAD_genomes_AF"} <= set(fields)
    aligned = (REPO / "resources/dbnsfp_transcript_aligned_columns.txt").read_text().splitlines()
    assert not {"TSL", "APPRIS"} & set(aligned)
    assert "${dbnsfp_fields}" in VEP_NF
    assert "vep_dbnsfp_fields" in (REPO / "run_annotate_vcf.sh").read_text()
