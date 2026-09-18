#!/usr/bin/env bash
# Pulls real SGE accounting data (maxvmem, cpu, wallclock, slots) for a list of
# job IDs, joining it back to the (gene, step) each one came from. Run this on
# the HPC login node.
#
# Input: a TSV of <gene>\t<step>\t<jobnumber> rows (see
#        src/tools/extract_run_job_ids.sh, or build one by hand).
#
# Usage: bash resources/pull_qacct_snapshot.sh <job_ids.tsv> <output.tsv>
# Takes a while for a large job list (~2000 jobs took several minutes) — safe
# to background:
#   nohup bash resources/pull_qacct_snapshot.sh RUNS/<RUN>/job_ids.tsv RUNS/<RUN>/qacct.tsv \
#       > RUNS/<RUN>/pull_qacct_snapshot.log 2>&1 &
#
# Output (and the input job-ID list) are run-specific data, not code — write
# them under RUNS/<RUN_NAME>/, which is gitignored, not under resources/.

set -u

in="${1:?Usage: $0 <job_ids.tsv> <output.tsv>}"
out="${2:?Usage: $0 <job_ids.tsv> <output.tsv>}"

printf "gene\tstep\tjobnumber\tjobname\tslots\tmaxvmem\tcpu\tru_wallclock\texit_status\n" > "$out"

total=$(wc -l < "$in")
i=0
while IFS=$'\t' read -r gene step jobnum; do
    i=$((i+1))
    if (( i % 200 == 0 )); then
        echo "[$i/$total] ..." >&2
    fi
    rec=$(qacct -j "$jobnum" 2>/dev/null)
    [ -z "$rec" ] && continue
    jobname=$(echo "$rec"    | awk '/^jobname/{print $2}')
    slots=$(echo "$rec"      | awk '/^slots/{print $2}')
    maxvmem=$(echo "$rec"    | awk '/^maxvmem/{print $2}')
    cpu=$(echo "$rec"        | awk '/^cpu/{print $2}')
    wallclock=$(echo "$rec"  | awk '/^ru_wallclock/{print $2}')
    exitst=$(echo "$rec"     | awk '/^exit_status/{print $2}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$gene" "$step" "$jobnum" "$jobname" "$slots" "$maxvmem" "$cpu" "$wallclock" "$exitst" >> "$out"
done < "$in"

echo "Done. $(wc -l < "$out") rows written to $out" >&2
