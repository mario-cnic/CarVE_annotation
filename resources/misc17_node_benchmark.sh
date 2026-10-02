#!/usr/bin/env bash
# MISC-17 node benchmark: is a node slower PER CPU-SECOND than the others? Run as an SGE job pinned to one node
# (HPC only; read-only except for its own temp dir). Compare the printed lines across nodes.
#
#   1. node facts: CPU model, logical CPUs, clock MHz seen by the kernel, load, other jobs' CPU use
#   2. fixed CPU work (pure Python loop + single-thread numpy matmul), timed as CPU seconds AND wall seconds
#   3. (optional) the REAL workload: the first N variants of one S223 chunk through src/python/annotate_spliceai.py,
#      same script/model/FASTA as the pipeline; run the SAME chunk on every node.
#
# Usage: misc17_node_benchmark.sh [CHUNK.vcf.gz N_VARIANTS FASTA]   (SpliceAI -d from $SPLICEAI_D, default 10000 = pipeline value)
# Output: stdout (use qsub -o). Lines starting "BENCH" are the comparable numbers.
set -u
PY=/data_lab_PGP/shared/utils/conda_envs/spliceai_env/bin/python3
HERE=/data_lab_PGP/pipelines/annotation_pipeline_new   # absolute: SGE runs a spooled COPY of this script

echo "BENCH host=$(hostname) date=$(date -Is) slots=${NSLOTS:-?} job=${JOB_ID:-?}"
echo "--- lscpu"; lscpu | grep -E "Model name|^CPU\(s\)|Thread|Core|Socket|MHz|NUMA node\(s\)|Hypervisor|Virtualization type"
echo "--- MHz per logical CPU (value x count)"; grep MHz /proc/cpuinfo | awk '{printf "%d\n",$4}' | sort -n | uniq -c
echo "--- governor / turbo"; cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo "no cpufreq governor file"
cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null | sed 's/^/no_turbo=/'
echo "--- load / memory"; cat /proc/loadavg; free -g | head -2
echo "--- top 8 CPU users on the node right now (other jobs)"; ps -eo user,pcpu,nlwp,comm --sort=-pcpu | head -9
echo "--- cgroup cpu limit of this job"; cat /proc/self/cgroup | head -3; cat /sys/fs/cgroup/cpu.max 2>/dev/null

echo "--- fixed CPU work"
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$PY" - <<'EOF'
import time, resource
def run(name, fn, reps=3):
    for i in range(reps):
        c0, w0 = time.process_time(), time.perf_counter()
        fn()
        print(f"BENCH {name} rep{i} cpu_s={time.process_time()-c0:.2f} wall_s={time.perf_counter()-w0:.2f}")
def pyloop():
    s = 0
    for i in range(30_000_000):
        s += i * i % 7
import numpy as np
a = np.random.default_rng(0).random((1500, 1500))
def matmul():
    for _ in range(3): a @ a
def memcopy():
    b = np.ones(200_000_000 // 8 * 8, dtype=np.float64)
    for _ in range(5): b.copy()
run("python_loop", pyloop); run("numpy_matmul_1thread", matmul); run("numpy_memcopy_1.6GB_x5", memcopy)
ru = resource.getrusage(resource.RUSAGE_SELF)
print(f"BENCH rusage utime={ru.ru_utime:.1f} stime={ru.ru_stime:.1f} nvcsw={ru.ru_nvcsw} nivcsw={ru.ru_nivcsw}")
EOF

if [ $# -ge 3 ]; then
    chunk="$1"; n="$2"; fasta="$3"
    tmp="$(mktemp -d "${TMPDIR:-/tmp}/misc17_bench.XXXXXX")"
    trap 'rm -rf "$tmp"' EXIT   # only the dir created above
    echo "--- real workload: first $n variants of $chunk"
    zcat "$chunk" | awk -v n="$n" '/^#/ {print; next} c<n {print; c++}' > "$tmp/in.vcf"
    /usr/bin/time -f "BENCH spliceai_real n=$n cpu_user_s=%U cpu_sys_s=%S wall_s=%e maxrss_kb=%M" \
        "$PY" "$HERE/src/python/annotate_spliceai.py" "$tmp/in.vcf" "$tmp/out.vcf" "$fasta" -d "${SPLICEAI_D:-10000}" \
        > "$tmp/stdout.log" 2>&1
    tail -n 3 "$tmp/stdout.log"
    echo "BENCH spliceai_real_out_lines=$(grep -vc '^#' "$tmp/out.vcf" 2>/dev/null || echo 0)"
fi
