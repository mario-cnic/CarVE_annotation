#!/bin/bash
# Compiles UTRAnnotatorPGP in the VEP container and checks its output is the same under several Perl hash seeds.
# Usage: tests/perl/run_in_vep_container.sh [path/to/vep.sif]
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SIF="${1:-/data_lab_PGP/resources/sif_images/vep.sif}"

singularity exec --cleanenv --bind /data_lab_PGP:/data_lab_PGP "$SIF" bash -s "$REPO_ROOT" <<'EOS'
set -euo pipefail
repo="$1"
export PERL5LIB="$repo/src/external/vep_plugins${PERL5LIB:+:$PERL5LIB}"
perl -c "$repo/src/external/vep_plugins/UTRAnnotatorPGP.pm"
perl -MUTRAnnotatorPGP -e 'UTRAnnotatorPGP->can("format_utr_effects") or die "format_utr_effects missing\n"'
ref=""
for seed in 0 1 2 3 4 5 6 7; do
    out=$(PERL_HASH_SEED=$seed PERL_PERTURB_KEYS=2 perl "$repo/tests/perl/test_utrannotator_pgp.pl")
    if [[ -z "$ref" ]]; then ref="$out"; elif [[ "$out" != "$ref" ]]; then echo "FAIL: output differs at seed $seed" >&2; exit 1; fi
done
echo "PASS: compiled; same output under 8 hash seeds"
echo "$ref"
EOS
