#!/usr/bin/env bash
#$ -P BIGN
#$ -N download_spliceai2
#$ -A PGP
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -o _log/download_spliceai2.stdout
#$ -e _log/download_spliceai2.stderr
#
# Download SpliceAI2 precomputed scores (GRCh38, v2.0) from Hugging Face at a
# pinned revision, verify sha256 against the LFS metadata, tabix-index each
# file and write a manifest. Re-running resumes: finished files are kept,
# existing indexes are not rebuilt.
#
# Usage (from the repo root):
#   qsub src/tools/download_spliceai2.sh --test     # indel chr21 only (0.16 GB)
#   qsub src/tools/download_spliceai2.sh            # all GRCh38 files (~348 GB)
# Options:
#   --out-dir DIR   destination root (default below)
#   --test          download only README, gene table and indel/chr21
#
# Output: $OUT_DIR/precomputed_scores_v2.0/GRCh38/{snv,indel}/chr*.tsv.gz(.tbi),
#         $OUT_DIR/SHA256SUMS[.test], $OUT_DIR/MANIFEST[.test].txt

set -euo pipefail

REPO="illumina-ai/SpliceAI2-data"
REVISION="cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9"
OUT_DIR="/data_lab_PGP/resources/annotation/spliceai2"
ENV_HF="/data_lab_PGP/shared/utils/conda_envs/whisper/bin"
ENV_TABIX="/data_lab_PGP/shared/utils/conda_envs/spliceai_env/bin"
THREADS="${NSLOTS:-4}"
TEST=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --test) TEST=1; shift ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

PREFIX="precomputed_scores_v2.0"
if [[ $TEST -eq 1 ]]; then
    INCLUDES=("$PREFIX/README.md" "$PREFIX/GRCh38/protein_coding_genes.tsv" "$PREFIX/GRCh38/indel/chr21.tsv.gz")
    SUFFIX=".test"
else
    INCLUDES=("$PREFIX/README.md" "$PREFIX/GRCh38/*")
    SUFFIX=""
fi

echo "Start $(date -Is) host=$(hostname) out=$OUT_DIR test=$TEST revision=$REVISION"
mkdir -p "$OUT_DIR"

# Fail early without network or credentials.
TOKEN_FILE="${HF_TOKEN_PATH:-${HF_HOME:-$HOME/.cache/huggingface}/token}"
if [[ -z "${HF_TOKEN:-}" && ! -s "$TOKEN_FILE" ]]; then
    echo "ERROR: no Hugging Face token at $TOKEN_FILE; run 'hf auth login' on the cluster once" >&2; exit 1
fi
"$ENV_HF/hf" auth whoami || { echo "ERROR: Hugging Face unreachable from $(hostname) (token present)" >&2; exit 1; }

# Expected checksums from the pinned revision.
SUMS="$OUT_DIR/SHA256SUMS$SUFFIX"
"$ENV_HF/python" -I - "$REPO" "$REVISION" "$TEST" > "$SUMS" <<'PY'
import sys
from huggingface_hub import HfApi
repo, rev, test = sys.argv[1], sys.argv[2], sys.argv[3] == "1"
keep = {"precomputed_scores_v2.0/GRCh38/indel/chr21.tsv.gz"}
for f in HfApi().list_repo_tree(repo, repo_type="dataset", revision=rev,
                                path_in_repo="precomputed_scores_v2.0/GRCh38", recursive=True):
    if getattr(f, "lfs", None) is None:
        continue
    if test and f.path not in keep:
        continue
    print(f"{f.lfs.sha256}  {f.path}")
PY
echo "Expected LFS files: $(wc -l < "$SUMS")"

INC_ARGS=()
for p in "${INCLUDES[@]}"; do INC_ARGS+=(--include "$p"); done
"$ENV_HF/hf" download "$REPO" --repo-type dataset --revision "$REVISION" \
    "${INC_ARGS[@]}" --local-dir "$OUT_DIR" --max-workers "$THREADS"

echo "Verifying sha256 $(date -Is)"
( cd "$OUT_DIR" && xargs -a "$SUMS" -d '\n' -P "$THREADS" -I{} sh -c 'echo "{}" | sha256sum -c --quiet -' ) \
    || { echo "ERROR: sha256 mismatch" >&2; exit 1; }
echo "sha256 OK"

echo "Indexing $(date -Is)"
awk '{print $2}' "$SUMS" | grep '\.tsv\.gz$' | while read -r f; do
    [[ -s "$OUT_DIR/$f.tbi" ]] && continue
    echo "$OUT_DIR/$f"
done | xargs -r -P "$THREADS" -I{} "$ENV_TABIX/tabix" -s 1 -b 2 -e 2 -S 1 {}

FIRST=$(awk '{print $2}' "$SUMS" | grep '\.tsv\.gz$' | head -1)
{
    echo "dataset	$REPO"
    echo "revision	$REVISION"
    echo "subset	$([[ $TEST -eq 1 ]] && echo test-indel-chr21 || echo GRCh38-all)"
    echo "assembly	GRCh38 (GENCODE v48 protein-coding genes, per dataset README)"
    echo "finished	$(date -Is)"
    echo "host	$(hostname)"
    echo "hf_cli	$("$ENV_HF/hf" version 2>&1 | head -1)"
    echo "tabix	$("$ENV_TABIX/tabix" --version | head -1)"
    echo "tabix_args	-s 1 -b 2 -e 2 -S 1"
    echo "files	$(wc -l < "$SUMS") (sha256 verified, list in $(basename "$SUMS"))"
    echo "bytes	$(awk '{print $2}' "$SUMS" | (cd "$OUT_DIR" && xargs du -cb | tail -1 | cut -f1))"
    echo "header	$( (zcat "$OUT_DIR/$FIRST" | head -1) 2>/dev/null || true)"
} > "$OUT_DIR/MANIFEST$SUFFIX.txt"
cat "$OUT_DIR/MANIFEST$SUFFIX.txt"
echo "First records of $FIRST:"
( zcat "$OUT_DIR/$FIRST" | head -3 ) 2>/dev/null || true
echo "Done $(date -Is)"
