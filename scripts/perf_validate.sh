#!/usr/bin/env bash
set -euo pipefail

FASTQ="${1:-data/hg002/SRR26901703_1.fastq}"
THREADS="${THREADS:-14}"
WARMUP_RUNS="${WARMUP_RUNS:-1}"
MEASURED_RUNS="${MEASURED_RUNS:-3}"
TMP_DIR="${TMP_DIR:-/tmp/zdash-perf-validate}"

if [[ ! -f "$FASTQ" ]]; then
  echo "missing FASTQ: $FASTQ" >&2
  exit 2
fi

mkdir -p "$TMP_DIR"

export ZIG_GLOBAL_CACHE_DIR="${ZIG_GLOBAL_CACHE_DIR:-/tmp/zig-global-cache}"
mkdir -p "$ZIG_GLOBAL_CACHE_DIR"

zig build -Doptimize=ReleaseFast >/dev/null

echo "=== Throughput Bench ==="
python3 - <<'PY' "$FASTQ" "$THREADS" "$WARMUP_RUNS" "$MEASURED_RUNS"
import os, sys, time, subprocess, statistics
fastq = sys.argv[1]
threads = int(sys.argv[2])
warmups = max(0, int(sys.argv[3]))
runs = max(1, int(sys.argv[4]))
size_gib = os.path.getsize(fastq) / (1024 ** 3)
cmds = {
    "check": f"./zig-out/bin/zdash check --quiet --threads={threads} {fastq} > /dev/null 2>&1",
    "scan": f"./zig-out/bin/zdash scan --quiet --threads={threads} {fastq} > /dev/null 2>&1",
    "stats": f"./zig-out/bin/zdash stats --quiet --threads={threads} {fastq} > /dev/null 2>&1",
}
def p90(values):
    if len(values) == 1:
        return values[0]
    idx = int(round(0.9 * (len(values) - 1)))
    return sorted(values)[idx]

print("mode,p50_s,p90_s,p50_gib_s,p90_gib_s,runs_s")
for mode, cmd in cmds.items():
    for _ in range(warmups):
        subprocess.run(cmd, shell=True, check=True)
    ts = []
    for _ in range(runs):
        t0 = time.perf_counter()
        subprocess.run(cmd, shell=True, check=True)
        t1 = time.perf_counter()
        ts.append(t1 - t0)
    med = statistics.median(ts)
    p90s = p90(ts)
    print(f"{mode},{med:.6f},{p90s:.6f},{(size_gib/med):.6f},{(size_gib/p90s):.6f}," + ";".join(f"{x:.6f}" for x in ts))
PY

if [[ -f "${FASTQ}.gz" ]]; then
  echo
  echo "=== Gzip Mode Check (${FASTQ}.gz) ==="
  python3 - <<'PY' "${FASTQ}.gz" "$THREADS"
import subprocess,time,sys
gz=sys.argv[1]
threads=int(sys.argv[2])
for mode in ("temp","stream"):
    cmd=f"./zig-out/bin/zdash check --quiet --threads={threads} --gzip-mode {mode} {gz} > /dev/null 2>&1"
    t0=time.perf_counter()
    subprocess.run(cmd,shell=True,check=True)
    t1=time.perf_counter()
    print(f"{mode},{t1-t0:.6f}")
PY
fi

echo
echo "=== Accuracy Checks ==="
./zig-out/bin/zdash check --json tests/fixtures/valid_small.fastq > "$TMP_DIR/check.json"
./zig-out/bin/zdash scan --json tests/fixtures/valid_small.fastq > "$TMP_DIR/scan.json"
./zig-out/bin/zdash stats --json tests/fixtures/valid_small.fastq > "$TMP_DIR/stats.json"
gzip -c tests/fixtures/valid_small.fastq > "$TMP_DIR/valid_small.fastq.gz"
./zig-out/bin/zdash check --json --gzip-mode stream "$TMP_DIR/valid_small.fastq.gz" > "$TMP_DIR/check_gz.json"
set +e
./zig-out/bin/zdash check --json --max-errors 2 tests/fixtures/bad_two_errors.fastq > "$TMP_DIR/bad.json"
EC=$?
set -e
if [[ "$EC" -ne 4 ]]; then
  echo "accuracy check failed: expected bad fixture exit 4, got $EC" >&2
  exit 1
fi
python3 - <<'PY' "$TMP_DIR/check.json" "$TMP_DIR/scan.json" "$TMP_DIR/stats.json" "$TMP_DIR/check_gz.json" "$TMP_DIR/bad.json"
import json, sys
c = json.load(open(sys.argv[1]))
s = json.load(open(sys.argv[2]))
t = json.load(open(sys.argv[3]))
g = json.load(open(sys.argv[4]))
b = json.load(open(sys.argv[5]))
assert c["status"] == "ok"
assert s["status"] == "ok"
assert t["status"] == "ok"
assert g["status"] == "ok"
assert c["total_reads"] == s["total_reads"] == t["total_reads"]
assert c["total_bases"] == s["total_bases"] == t["total_bases"]
assert c["total_reads"] == g["total_reads"]
assert c["total_bases"] == g["total_bases"]
assert b["status"] == "failed"
print("accuracy_ok=true")
PY

echo "artifacts=$TMP_DIR"
