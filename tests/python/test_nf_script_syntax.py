"""Every process `script:` block in modules/local/*.nf, rendered with dummy values, is valid bash (MISC-15 lint).

A -stub-run never executes script blocks, so a shell error there only shows on the cluster. Groovy
interpolations `${...}` are replaced by a dummy word, Groovy escapes are undone, and `bash -n` checks the
result. `bash -n` accepts a line starting with `//` (it parses as a command), which is what broke
PANGOLIN_ANNOTATE (MISC-15), so such lines are checked separately.
"""
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src/python"))
import annotation_cache as ac  # noqa: E402

MODULES = sorted((REPO / "modules/local").glob("*.nf"))
PROCESSES = [(nf, name) for nf in MODULES for name in re.findall(r"^process\s+(\w+)\s*\{", nf.read_text(), re.M)]


def strip_interpolations(text: str) -> str:
    """Replaces each ${...} (nested braces allowed) with a dummy word."""
    out, i = [], 0
    while i < len(text):
        if text.startswith("${", i):
            depth, j = 0, i + 1
            while j < len(text):
                depth += {"{": 1, "}": -1}.get(text[j], 0)
                if depth == 0:
                    break
                j += 1
            out.append("DUMMY")
            i = j + 1
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


def render(block: str) -> str:
    body = block.split('"""')[1]
    body = strip_interpolations(body)
    return body.replace("\\\\", "\0").replace("\\$", "$").replace("\0", "\\")


def test_every_module_has_processes():
    assert len(PROCESSES) >= 17


@pytest.mark.parametrize("nf,process", PROCESSES, ids=[p for _, p in PROCESSES])
def test_script_block_is_valid_bash(nf, process):
    script = render(ac.script_block(nf.read_text(), process))
    res = subprocess.run(["bash", "-n"], input=script, capture_output=True, text=True)
    assert res.returncode == 0, f"{nf.name} {process}: {res.stderr}\n{script}"
    groovy_comments = [l for l in script.splitlines() if l.strip().startswith("//")]
    assert not groovy_comments, f"{nf.name} {process}: Groovy comment inside the bash script: {groovy_comments}"


CODE_KEYS = set(re.findall(r"^\s+(\w+)\s*:\s*codeMd5\(", (REPO / "annotate_vcf.nf").read_text(), re.M))


@pytest.mark.parametrize("nf,process", PROCESSES, ids=[p for _, p in PROCESSES])
def test_process_running_repo_code_names_its_code_md5(nf, process):
    # TODO ordered item 6: -resume must rerun a task whose code changed; Nextflow only sees it through the
    # `# code ${params.code_md5.<key>}` line (the md5 is computed in annotate_vcf.nf).
    block = ac.script_block(nf.read_text(), process)            # without the code line
    full = re.search(r"^process\s+" + process + r"\s*\{(.*?)^}", nf.read_text(), re.S | re.M).group(1)
    runs_code = re.search(r"params\.\w+_(script|python)\b|params\.pangolin_repo|params\.vep_dbnsfp_fields", block)
    keys = re.findall(r"# code \$\{params\.code_md5\.(\w+)\}", full)
    if runs_code:
        assert len(keys) == 1 and keys[0] in CODE_KEYS, f"{process}: {keys} not in {sorted(CODE_KEYS)}"
    else:
        assert not keys


def test_cache_identity_ignores_the_code_line():
    text = 'process P {\n    script:\n    """\n    # code ${params.code_md5.spliceai}\n    run it\n    """\n}\n'
    assert ac.script_block(text, "P") == ac.script_block(text.replace("    # code ${params.code_md5.spliceai}\n", ""), "P")
