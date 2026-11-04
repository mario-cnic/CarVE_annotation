"""Every script the pipeline launches imports as in production (C3).

The other tests put src/python (and, through config.py, shared/utils/src) on sys.path before importing,
which hid C1: filter_variants.py failed to import in production while every test passed. Here each
script runs as its own process, `python <script> --help`, from an unrelated directory, without
PYTHONPATH and without importing config, as Nextflow and run_annotate_vcf.sh call it.
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

# Scripts referenced by nextflow.config (*_script params) and run_annotate_vcf.sh.
SCRIPTS = [
    "src/python/parsing/filter_variants.py",
    "src/python/parsing/vcf_parser_pysam.py",
    "src/python/annotate_branchpointer.py",
    "src/python/annotation_cache.py",
    "src/python/check_vcf_assembly.py",
    "src/python/extract_genotypes.py",
    "src/python/genotypes_to_wide.py",
    "src/python/split_vcf_chunks.py",
    "src/python/tag_transcript_priority.py",
    "src/python/write_run_manifest.py",
    "src/python/build_run_qc_report.py",
]
# Needs the spliceai/tensorflow env, not installed in the test env; covered by test_annotate_spliceai.py.
NOT_HERE = ["src/python/annotate_spliceai.py"]
# Third-party packages a script needs from its own conda env (the script is skipped where they are missing).
EXTERNAL = {"src/python/parsing/vcf_parser_pysam.py": "vcf_parser"}


def test_list_matches_the_pipeline():
    config = (REPO / "nextflow.config").read_text()
    launcher = (REPO / "run_annotate_vcf.sh").read_text()
    referenced = {s for s in SCRIPTS + NOT_HERE if s in config.replace("${projectDir}/", "")
                  or s in launcher.replace("$REPO_ROOT/", "")}
    assert referenced == set(SCRIPTS + NOT_HERE)


@pytest.mark.parametrize("script", SCRIPTS)
def test_script_imports_as_in_production(script, tmp_path):
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME")}
    if script in EXTERNAL:
        probe = subprocess.run([sys.executable, "-c", f"import {EXTERNAL[script]}"], cwd=tmp_path, env=env,
                               capture_output=True)
        if probe.returncode:
            pytest.skip(f"{EXTERNAL[script]} not installed in this env")
    res = subprocess.run([sys.executable, str(REPO / script), "--help"], cwd=tmp_path, env=env,
                         capture_output=True, text=True, timeout=120)
    assert res.returncode == 0, res.stderr[-2000:]
    assert "usage" in res.stdout.lower()
