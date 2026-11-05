"""run_annotate_vcf.sh with stub `nextflow` and `python3` on PATH: the score-cache segment list a run reads."""
import os
import stat
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

NEXTFLOW = """#!/bin/bash
if [[ "$1" == config ]]; then echo "annotation_cache_dir = '$FAKE_CACHE'"; exit 0; fi
printf '%s\\n' "$@" > "$FAKE_LOG"
"""
# python3 stand-in: version-id prints $FAKE_VID; list-segments writes one segment named after $FAKE_TAG.
PYTHON3 = """#!/bin/bash
case "$1" in
  *write_run_manifest.py)
    while [[ $# -gt 0 ]]; do [[ "$1" == --outdir ]] && mkdir -p "$2/pipeline_info"; shift; done; exit 0 ;;
  *annotation_cache.py)
    if [[ "$2" == version-id ]]; then echo "$FAKE_VID"; exit 0; fi
    if [[ "$2" == list-segments ]]; then
      while [[ $# -gt 0 ]]; do [[ "$1" == --output ]] && echo "/cache/$FAKE_TAG.vcf.gz" > "$2"; shift; done; echo 1; exit 0
    fi ;;
esac
exit 0
"""


def launch(tmp_path, vid, tag, *extra):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    for name, text in (("nextflow", NEXTFLOW), ("python3", PYTHON3)):
        (bin_dir / name).write_text(text)
        (bin_dir / name).chmod(0o755)
    (tmp_path / "in.vcf.gz").touch()
    env = dict(os.environ, PATH=f"{bin_dir}:{os.environ['PATH']}", FAKE_VID=vid, FAKE_TAG=tag,
               FAKE_CACHE=str(tmp_path / "cache"), FAKE_LOG=str(tmp_path / "nf_args.txt"))
    subprocess.run(["bash", str(REPO / "run_annotate_vcf.sh"), "--input_vcf", "in.vcf.gz", "--outdir", "out", *extra],
                   cwd=tmp_path, env=env, check=True, capture_output=True, text=True)
    args = (tmp_path / "nf_args.txt").read_text().splitlines()
    seg_list = args[args.index("--annotation_cache_segments_spliceai") + 1]
    return args[args.index("--annotation_cache_version_spliceai") + 1], Path(seg_list).read_text().strip()


def test_resume_keeps_the_segment_list_of_the_same_version(tmp_path):
    assert launch(tmp_path, "v1", "first") == ("v1", "/cache/first.vcf.gz")
    # same scoring code: the frozen list of the first launch is reused, not re-listed
    assert launch(tmp_path, "v1", "second", "-resume") == ("v1", "/cache/first.vcf.gz")


def test_resume_after_a_code_change_lists_the_new_versions_segments(tmp_path):
    launch(tmp_path, "v1", "first")
    # scoring code changed between launches: new version id, so its own segments, not v1's
    assert launch(tmp_path, "v2", "second", "-resume") == ("v2", "/cache/second.vcf.gz")
