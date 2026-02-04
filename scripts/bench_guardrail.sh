#!/usr/bin/env bash
set -euo pipefail

FASTQ="${1:-data/hg002/SRR26901703_1.fastq}"
THREADS="${THREADS:-14}"
MIN_CHECK_GIB_S="${MIN_CHECK_GIB_S:-0.50}"
WARMUP_RUNS="${WARMUP_RUNS:-1}"
MEASURED_RUNS="${MEASURED_RUNS:-3}"

if [[ ! -f "$FASTQ" ]]; then
  echo "missing FASTQ: $FASTQ" >&2
  exit 2
fi

zig build -Doptimize=ReleaseFast >/dev/null

read -r SEC GIB_S < <(python3 - <<'PY' "$FASTQ" "$THREADS" "$WARMUP_RUNS" "$MEASURED_RUNS"
import os, sys, time, subprocess, statistics
fastq=sys.argv[1]
threads=int(sys.argv[2])
warmups=max(0, int(sys.argv[3]))
runs=max(1, int(sys.argv[4]))
size=os.path.getsize(fastq)/(1024**3)
cmd=f'./zig-out/bin/zdash check --quiet --threads={threads} {fastq} > /dev/null 2>&1'
for _ in range(warmups):
    subprocess.run(cmd,shell=True,check=True)
times=[]
for _ in range(runs):
    t0=time.perf_counter(); subprocess.run(cmd,shell=True,check=True); t1=time.perf_counter()
    times.append(t1-t0)
dt=statistics.median(times)
print(f"{dt:.6f} {size/dt:.6f}")
PY
)

echo "zdash-check median (${MEASURED_RUNS} runs, ${WARMUP_RUNS} warmup): ${SEC}s (${GIB_S} GiB/s)"

awk -v v="$GIB_S" -v min="$MIN_CHECK_GIB_S" 'BEGIN { exit !(v+0 >= min+0) }'
echo "guardrail passed: ${GIB_S} >= ${MIN_CHECK_GIB_S}"
