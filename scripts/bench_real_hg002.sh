#!/usr/bin/env bash
set -euo pipefail

ACCESSION="${1:-SRR26901703}"
OUT_DIR="${2:-data/hg002}"
RUNS="${RUNS:-5}"
THREADS="${THREADS:-$(python3 - <<'PY'
import os
print(os.cpu_count() or 1)
PY
)}"

mkdir -p "$OUT_DIR"

echo "Resolving ENA FASTQ links for accession: $ACCESSION"
ENA_TSV_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACCESSION}&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes&format=tsv"
ENA_TSV="$(curl -fsSL "$ENA_TSV_URL")"

if [[ "$(printf '%s\n' "$ENA_TSV" | wc -l | tr -d ' ')" -lt 2 ]]; then
  echo "ERROR: accession not found in ENA filereport: $ACCESSION" >&2
  exit 1
fi

LINE="$(printf '%s\n' "$ENA_TSV" | sed -n '2p')"
FASTQ_FTP="$(printf '%s\n' "$LINE" | cut -f2)"
FASTQ_BYTES="$(printf '%s\n' "$LINE" | cut -f3)"
R1_FTP="$(printf '%s\n' "$FASTQ_FTP" | cut -d';' -f1)"
R1_BYTES="$(printf '%s\n' "$FASTQ_BYTES" | cut -d';' -f1)"
R1_HTTPS="https://${R1_FTP}"
R1_GZ="${OUT_DIR}/$(basename "$R1_FTP")"
R1_FASTQ="${R1_GZ%.gz}"

echo "R1 URL: $R1_HTTPS"
echo "R1 compressed bytes (ENA): $R1_BYTES"

if [[ ! -f "$R1_GZ" ]]; then
  echo "Downloading: $R1_GZ"
  curl -fL --retry 5 --retry-delay 2 -o "$R1_GZ" "$R1_HTTPS"
else
  echo "Reusing existing download: $R1_GZ"
fi

if [[ ! -f "$R1_FASTQ" ]]; then
  echo "Decompressing: $R1_FASTQ"
  gzip -dc "$R1_GZ" > "$R1_FASTQ"
else
  echo "Reusing existing decompressed file: $R1_FASTQ"
fi

echo "Building zdash in ReleaseFast"
zig build -Doptimize=ReleaseFast

echo "Benchmarking on: $R1_FASTQ"
python3 - "$R1_FASTQ" "$THREADS" "$RUNS" <<'PY'
import os, subprocess, time, statistics, sys, shutil

fastq = sys.argv[1]
threads = int(sys.argv[2])
runs_n = int(sys.argv[3])
size = os.path.getsize(fastq)
gib = size / (1024 ** 3)

tools = {
    "zdash-full": f"zig-out/bin/zdash --quiet --threads={threads} --chunk-bytes=8388608 --profile=full {fastq} > /dev/null 2>&1",
    "zdash-validate-stats": f"zig-out/bin/zdash --quiet --threads={threads} --chunk-bytes=8388608 --profile=validate-stats {fastq} > /dev/null 2>&1",
    "zdash-stats-only": f"zig-out/bin/zdash --quiet --threads={threads} --chunk-bytes=8388608 --profile=stats-only {fastq} > /dev/null 2>&1",
    "seqtk": f"seqtk seq {fastq} > /dev/null 2>&1",
    "seqkit": f"seqkit stats {fastq} > /dev/null 2>&1",
    "fastp": f"fastp -i {fastq} -o /dev/null -A -Q -L -w {threads} -j /tmp/fastp_bench.json -h /tmp/fastp_bench.html > /dev/null 2>&1",
}

noodles_bin = "/tmp/noodles_bench/target/release/noodles_bench"
if os.path.exists(noodles_bin):
    tools["noodles-rust"] = f"{noodles_bin} {fastq} > /dev/null 2>&1"

print(f"file={fastq}")
print(f"file_bytes={size}")
print(f"threads={threads}")
print(f"runs={runs_n}")
print()

available = {}
for name, cmd in tools.items():
    exe = cmd.split()[0]
    if exe.startswith("/") or shutil.which(exe):
        available[name] = cmd
    else:
        print(f"skip {name}: missing executable '{exe}'")

print()
print("tool,avg_s,best_s,avg_gib_s,best_gib_s,runs")
for name, cmd in available.items():
    subprocess.run(cmd, shell=True, check=True)  # warmup
    runs = []
    for _ in range(runs_n):
        t0 = time.perf_counter()
        subprocess.run(cmd, shell=True, check=True)
        t1 = time.perf_counter()
        runs.append(t1 - t0)
    avg = statistics.mean(runs)
    best = min(runs)
    print(
        f"{name},{avg:.6f},{best:.6f},{(gib/avg):.3f},{(gib/best):.3f}," +
        ";".join(f"{x:.6f}" for x in runs)
    )
PY

echo
echo "Done."
