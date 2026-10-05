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
