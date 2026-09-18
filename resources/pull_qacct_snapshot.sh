#!/usr/bin/env bash
# Pulls real SGE accounting data (maxvmem, cpu, wallclock, slots) for every job ID
# recorded in run_more_genes_20260817_1143_job_ids.tsv, joining it back to the
# (gene, step) it came from. Run this on the HPC login node.
#
# Usage: bash resources/pull_qacct_snapshot.sh
# Takes a while (~2000 individual qacct calls) — safe to background:
#   nohup bash resources/pull_qacct_snapshot.sh > resources/pull_qacct_snapshot.log 2>&1 &

set -u
cd "$(dirname "$0")/.."

in="resources/run_more_genes_20260817_1143_job_ids.tsv"
out="resources/run_more_genes_20260817_1143_qacct.tsv"

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
