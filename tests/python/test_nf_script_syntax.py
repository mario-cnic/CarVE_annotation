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
