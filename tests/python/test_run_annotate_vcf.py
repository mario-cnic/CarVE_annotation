"""run_annotate_vcf.sh with stub `nextflow` and `python3` on PATH: the score-cache segment list a run reads."""
import os
import stat
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# nextflow stand-in: the --dump_params run writes the params the wrapper resolves (a --fasta value is
# echoed back, as Nextflow would resolve it); the real run records its arguments.
NEXTFLOW = """#!/bin/bash
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
  [[ "${args[$i]}" == --dump_params ]] && dump="${args[$((i + 1))]}"
  [[ "${args[$i]}" == --fasta ]] && fasta="${args[$((i + 1))]}"
done
if [[ -n "$dump" ]]; then
  printf '{"annotation_cache_dir": "%s", "fasta": "%s"}' "$FAKE_CACHE" "${fasta:-/config/ref.fa}" > "$dump"; exit 0
fi
printf '%s\\n' "$@" > "$FAKE_LOG"
"""
# python3 stand-in: version-id prints $FAKE_VID; list-segments writes one segment named after $FAKE_TAG.
PYTHON3 = """#!/bin/bash
[[ "$1" == -c ]] && exec "$REAL_PYTHON" "$@"
case "$1" in
  *write_run_manifest.py)
    printf '%s\\n' "$@" > "$FAKE_LOG.manifest_$2"
    while [[ $# -gt 0 ]]; do [[ "$1" == --outdir ]] && mkdir -p "$2/pipeline_info"; shift; done; exit 0 ;;
  *annotation_cache.py)
    if [[ "$2" == version-id ]]; then echo "$FAKE_VID"; exit 0; fi
    if [[ "$2" == list-segments ]]; then
      while [[ $# -gt 0 ]]; do [[ "$1" == --output ]] && echo "/cache/$FAKE_TAG.vcf.gz" > "$2"; shift; done; echo 1; exit 0
    fi ;;
esac
exit 0
"""


def launch(tmp_path, vid, tag, *extra, cwd=None):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    for name, text in (("nextflow", NEXTFLOW), ("python3", PYTHON3)):
        (bin_dir / name).write_text(text)
        (bin_dir / name).chmod(0o755)
    (tmp_path / "in.vcf.gz").touch()
    env = dict(os.environ, PATH=f"{bin_dir}:{os.environ['PATH']}", FAKE_VID=vid, FAKE_TAG=tag, REAL_PYTHON=sys.executable,
               FAKE_CACHE=str(tmp_path / "cache"), FAKE_LOG=str(tmp_path / "nf_args.txt"))
    subprocess.run(["bash", str(REPO / "run_annotate_vcf.sh"), "--input_vcf", str(tmp_path / "in.vcf.gz"), "--outdir",
                    str(tmp_path / "out"), *extra], cwd=cwd or tmp_path, env=env, check=True, capture_output=True, text=True)
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


def test_command_line_params_reach_the_manifest_and_the_run(tmp_path):
    # MISC-42: the manifest took values from `nextflow config`, which ignores --param values; a relative
    # path was resolved by nextflow from the repo root.
    (tmp_path / "ref.fa").write_text(">chr1\nACGT\n")
    launch(tmp_path, "v1", "first", "--fasta", "ref.fa")
    manifest_args = (tmp_path / "nf_args.txt.manifest_launch").read_text().splitlines()
    run_args = (tmp_path / "nf_args.txt").read_text().splitlines()
    assert f"fasta={tmp_path}/ref.fa" in manifest_args
    assert run_args[run_args.index("--fasta") + 1] == f"{tmp_path}/ref.fa"
    assert run_args[run_args.index("-profile") + 1] == "standard"   # the default profile is passed, as recorded
