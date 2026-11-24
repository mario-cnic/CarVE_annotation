"""SPiP fields are read by the names of the SPiP v2.1 header (PRED-5): the config's expected header equals the
parser's field list, SPIP_ANNOTATE's real script block fails on a reordered header, and the parser fails on an
entry with another field count."""
import gzip
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src/python"))
import annotation_cache as ac  # noqa: E402
from parsing.modules import splicing as sp  # noqa: E402

FORMAT = "|".join(sp.SPIP_FIELDS)


def test_config_header_equals_parser_fields():
    m = re.search(r"^\s*spip_info_format\s*=\s*'([^']*)'", (REPO / "nextflow.config").read_text(), re.M)
    assert m and m.group(1) == FORMAT
    assert len(sp.SPIP_FIELDS) == 37 and len(set(sp.SPIP_FIELDS)) == 37


def test_entry_with_other_field_count_fails():
    good = "|".join(["A", "id", "NTR", "00 %", "0.1"] + [""] * 32)
    assert len(sp._spip_entries(good + "," + good)) == 2
    with pytest.raises(ValueError, match="36 fields, expected 37"):
        sp._spip_entries(good + "," + good.rsplit("|", 1)[0])


def test_fields_unpacked_by_name():
    f = dict.fromkeys(sp.SPIP_FIELDS, "")
    f.update(Interpretation="Alter BP", InterConfident="40 %", SPiPscore="0.62", mutInPBarea="Yes",
             deltaESRscore="-0.2", probaCryptMut="0.1", classProbaCryptMut="No")
    out = sp._spip_fields(list(f.values()))
    assert out == ["Alter BP", "40 %", 0.62, "40 %", 0.62, "bp_disruption", "Yes", -0.2, 0.1, "scored"]


TOOLS = {t: shutil.which(t) for t in ("bgzip", "tabix")}


def render(script: str, **values) -> str:
    for k, v in values.items():
        script = script.replace("${" + k + "}", v)
    return script.replace("\\$", "$").replace("\\\\\n", "\\\n")


@pytest.mark.skipif(not all(TOOLS.values()), reason="bgzip/tabix not available")
@pytest.mark.parametrize("fields, ok", [(sp.SPIP_FIELDS, True),
                                        (sp.SPIP_FIELDS[:3] + (sp.SPIP_FIELDS[4], sp.SPIP_FIELDS[3]) + sp.SPIP_FIELDS[5:], False)])
def test_spip_annotate_checks_header(tmp_path, fields, ok):
    header = ("##fileformat=VCFv4.2\n##SPiP output v2.1\n"
              f"##INFO=<ID=SPiP,Number=.,Type=String,Description=\"SPiPs predictions of splicing defects. "
              f"Format: '{'|'.join(fields)}'\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    fake = tmp_path / "fake_rscript"
    fake.write_text('#!/bin/bash\nwhile [ $# -gt 0 ]; do [ "$1" = --output ] && out="$2"; shift; done\n'
                    f"cat > \"$out\" <<'EOF'\n{header}EOF\n")
    fake.chmod(0o755)
    with gzip.open(tmp_path / "in.vcf.gz", "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    block = ac.script_block((REPO / "modules/local/spip.nf").read_text(), "SPIP_ANNOTATE")
    body = block.split('"""')[1]
    script = render(body, **{"raw_out": "p.annSPiP.vcf", "vcf": "in.vcf.gz", "params.spip_rscript": str(fake),
                             "params.spip_script": "spip.r", "task.cpus": "1", "params.spip_max_lines": "10",
                             "params.spip_info_format": FORMAT, "params.bgzip": TOOLS["bgzip"],
                             "params.tabix": TOOLS["tabix"], "meta.partition_id": "p"})
    res = subprocess.run(["bash", "-ue", "-c", script], cwd=tmp_path, capture_output=True, text=True)
    if ok:
        assert res.returncode == 0, res.stderr
        assert (tmp_path / "p.annSPiP.vcf.gz.tbi").exists()
    else:
        assert res.returncode == 1 and "SPiP INFO header fields differ" in res.stderr
        assert not (tmp_path / "p.annSPiP.vcf.gz").exists()
