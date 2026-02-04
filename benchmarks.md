# Z-DASH Benchmark Notes

Date: February 4, 2026  
Machine: macOS (aarch64), 12 logical CPUs, Zig 0.15.2

This file tracks recent competitive benchmarks and methodology.

Quick local validation command (throughput + correctness):

```bash
WARMUP_RUNS=1 MEASURED_RUNS=3 scripts/perf_validate.sh data/hg002/SRR26901703_1.fastq
```

This emits p50/p90 timings + GiB/s for `check/scan/stats` and runs fixture-based accuracy checks.

## Dataset

- Path: `/tmp/zdash_bench.fastq`
- Size: `125,088,890` bytes (~0.116 GiB)
- Shape: `400,000` reads, length `150`
- Type: synthetic FASTQ (uncompressed)

Generation command used:

```bash
python3 - <<'PY'
import random, os
path='/tmp/zdash_bench.fastq'
random.seed(42)
reads=400000
seq_len=150
bases='ACGTN'
quals=[chr(33+i) for i in range(41)]
with open(path,'w') as f:
    for i in range(reads):
        seq=''.join(random.choice(bases) for _ in range(seq_len))
        q=''.join(random.choice(quals) for _ in range(seq_len))
        f.write(f'@r{i}\n{seq}\n+\n{q}\n')
print(path, os.path.getsize(path))
PY
```

## Real HG002 Dataset Source (recommended)

Use ENA (European Nucleotide Archive) for stable direct FASTQ links.

- ENA search API (sample alias HG002):
  - `https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=sample_alias%3D%22HG002%22%20AND%20instrument_platform%3DILLUMINA&fields=run_accession,study_accession,sample_alias,library_layout,fastq_ftp,fastq_bytes&limit=20&format=tsv`
- Example accession used in scripts:
  - `SRR26901703` (paired-end; script benchmarks R1)

Convenience script added:

```bash
scripts/bench_real_hg002.sh SRR26901703 data/hg002
```

This script:
- resolves ENA download URL
- downloads `*.fastq.gz`
- decompresses to `*.fastq`
- builds `zdash` in ReleaseFast
- runs multi-tool benchmark and prints CSV-style results

## Build Mode

`zdash` was benchmarked in ReleaseFast:

```bash
zig build -Doptimize=ReleaseFast
```

## Recent Optimization Pass (implemented before this run)

- removed chunk pre-pass recounting (`countLines`) and switched to single-pass chunk metadata build
- added `--quiet` mode to remove summary I/O from timed runs
- added profile modes for fairer comparisons:
  - `full`
  - `validate-stats`
  - `stats-only`
- switched low-quality decision from per-read float division to integer threshold compare
- changed worker scheduling to hybrid mode:
  - static ranges for smaller chunk counts
  - atomic batched fetch for larger chunk counts
- added mmap sequential access hint (`madvise` where available)

## Compared Commands

- `zdash-full`: `zig-out/bin/zdash --quiet --threads=12 --chunk-bytes=8388608 --profile=full /tmp/zdash_bench.fastq`
- `zdash-validate-stats`: `zig-out/bin/zdash --quiet --threads=12 --chunk-bytes=8388608 --profile=validate-stats /tmp/zdash_bench.fastq`
- `zdash-stats-only`: `zig-out/bin/zdash --quiet --threads=12 --chunk-bytes=8388608 --profile=stats-only /tmp/zdash_bench.fastq`
- `seqtk`: `seqtk seq /tmp/zdash_bench.fastq > /dev/null`
- `seqkit`: `seqkit stats /tmp/zdash_bench.fastq > /dev/null`
- `fastp`: `fastp -i /tmp/zdash_bench.fastq -o /dev/null -A -Q -L -w 12 -j /tmp/fastp_bench.json -h /tmp/fastp_bench.html`
- `noodles-rust`: custom Rust benchmark binary using `noodles-fastq` crate (`/tmp/noodles_bench/target/release/noodles_bench /tmp/zdash_bench.fastq`)

Note on noodles: there is no standard noodles FASTQ sanitizer CLI directly equivalent to `zdash`; this comparison uses a custom parser benchmark executable.

## Results (5 measured runs each, ReleaseFast build)

| Tool | Avg Time (s) | Best Time (s) | Avg GiB/s | Best GiB/s |
|---|---:|---:|---:|---:|
| zdash-full | 0.0702 | 0.0692 | 1.660 | 1.683 |
| zdash-validate-stats | 0.0681 | 0.0661 | 1.710 | 1.762 |
| zdash-stats-only | 0.0677 | 0.0667 | 1.721 | 1.747 |
| seqtk | 0.1023 | 0.1010 | 1.139 | 1.153 |
| seqkit | 0.0867 | 0.0844 | 1.344 | 1.381 |
| fastp | 0.3256 | 0.3172 | 0.358 | 0.367 |
| noodles-rust | 0.2123 | 0.2105 | 0.549 | 0.553 |

## Interpretation

- In this run on this dataset, `zdash` outperforms the tested baselines on throughput.
- `zdash` profile differences are small on this synthetic file; validation checks are not dominating runtime after current optimizations.
- Previous slow results were largely due to benchmarking debug builds instead of ReleaseFast.

## Caveats

- This is one synthetic dataset and one machine.
- `seqkit stats` is not a strict sanitizer equivalent; profile-based comparisons improve fairness but are still not perfect.
- A larger realistic dataset (>=10GB, mixed quality/base distributions, possibly gz workflows) is still needed before broad claims.

## Real Dataset Profiling Update (HG002 R1, Feb 4, 2026)

Command:

```bash
./zig-out/bin/zdash --quiet --profile-internal --threads 10 data/hg002/SRR26901703_1.fastq
```

Top 3 internal hotspots:

1. `validateFastqChunk`: 193.680s (98.2%)
2. `analyzeSequenceKernel`: 2.290s (1.2%)
3. `analyzeQualityKernel`: 1.244s (0.6%)

Observed end-to-end runtime (same dataset, same thread count):

- `47.31s real` (`/usr/bin/time -l ./zig-out/bin/zdash --quiet --threads 10 data/hg002/SRR26901703_1.fastq`)

Additional local tuning runs:

- `49.13s real` (`--threads 10`) after parser-path refactor iteration
- `43.67s real` (`--threads 12`) on the same file and build (best in this pass)
- `45.40s real` (`--threads 12 --chunk-bytes 1048576`)
- `48.75s real` (`--threads 12 --chunk-bytes 67108864`)
- `44.74s real` (`--threads 12`) after introducing SIMD-batched newline offset scanning
- `46.18s real` (`--threads 10`) after introducing SIMD-batched newline offset scanning
- `43.19s real` (`--threads 12`) after adding speculative 4-line fast-path on top of batched scanning
- `58.47s real` (`--threads 6`) and `52.56s real` (`--threads 8`) in thread sweep (physical-core-only was slower here)
- `42.68s real` (`--threads 14`) and `43.95s real` (`--threads 16`) in thread sweep under dynamic scheduling

Takeaway:

- SIMD math is not the bottleneck; record framing/control flow inside chunk validation dominates.
- current best local setting in this pass is `--threads 12` with default chunk size.
- batched newline scanning + speculative framing recovered and slightly improved the best local `--threads 12` timing.

## Scheduler Update (Static Contiguous Ranges)

Change:

- switched multithread chunk scheduling from hybrid/dynamic to static contiguous ranges for this workload

Why:

- with `mmap`, contiguous per-thread regions improved readahead/page locality and reduced cross-thread page churn

Observed runs after this change:

- `44.03s real` (`--threads 12`)
- `39.93s real` (`--threads 12`, repeat)
- `37.35s real` and `39.35s real` (`--threads 14`, repeats)

Current upshot:

- this is the first change that materially moved real HG002 runtime below the mid-40s in repeated local runs
- still not at competitor tier yet, but closes a meaningful part of the gap

## Hot-Path Fusion Update

Change:

- added `analyzeSequenceHot` / `analyzeQualityHot` for speculative and fallback loops
- these bypass profiler wrapper overhead in normal runs while preserving profiler behavior when `--profile-internal` is enabled
- added lightweight prefetch on upcoming record starts when buffered newline offsets are available

Observed runs after this change:

- `45.45s real` (`--threads 12`)
- `34.77s real` (`--threads 14`)
- `35.59s real` (`--threads 14`, repeat)

Internal profile snapshot (`--threads 14`, `--profile-internal`):

- `validateFastqChunk`: ~62.1%
- `readLine` (newline refill): ~36.3%
- `analyzeSequenceKernel`: ~1.0%

Current upshot:

- best repeated local runtime is now in the mid-30s on this HG002 file
- remaining gap is still primarily framing/control flow, but overall throughput moved meaningfully closer to competitor territory

## Milestone 6B Update: Producer Thread + Assume-Valid Mode

Changes:

- added parser mode flag: `--mode strict|assume-valid`
- `assume-valid` uses fused seq+qual accumulation and skips per-character validation
- enabled producer-consumer execution in `assume-valid` using a lock-free chunk queue
- producer thread prefetches upcoming chunk starts and feeds worker queue

Latest 1-run results (no warmup, HG002 R1, threads=14):

| Tool | Time (s) | GiB/s |
|---|---:|---:|
| zdash-strict-full | 42.935 | 0.399 |
| zdash-assume-full | 8.103 | 2.113 |
| zdash-assume-stats | 7.830 | 2.186 |
| seqtk | 16.891 | 1.014 |
| seqkit | 5.900 | 2.902 |
| fastp | 20.299 | 0.843 |
| noodles-rust | 21.108 | 0.811 |

Upshot:

- strict mode remains slower than top parsers on this dataset
- assume-valid mode crossed 2 GiB/s and now beats `seqtk`, `fastp`, and `noodles-rust` in this run
- `seqkit` still leads in this no-warmup run

## Strict Loop + Case-Fold Validation Pass

Changes:

- rewrote strict loop to operate on newline offsets (`fillTo(4)`) instead of repeated `nextLine` hot-path calls
- kept tail-only line parsing fallback for partial trailing records
- switched sequence validation to ASCII case-fold checks in SIMD/scalar strict kernels
- normalized strict read-index accounting to `read_index_offset + total_reads + 1` before increment

Latest 1-run results (no warmup, HG002 R1, threads=14):

| Tool | Time (s) | GiB/s |
|---|---:|---:|
| zdash-strict-full | 33.614 | 0.509 |
| zdash-assume-full | 15.026 | 1.139 |
| zdash-assume-stats | 13.040 | 1.313 |
| seqtk | 18.157 | 0.943 |
| seqkit | 5.591 | 3.062 |
| fastp | 19.514 | 0.877 |
| noodles-rust | 20.668 | 0.828 |

Upshot:

- strict mode improved materially in this pass and now beats `fastp`/`noodles` in this single run
- assume-valid remains faster than strict but slower than the prior best warm-cache spike
- `seqkit` remains the leader

## Streaming Chunk Producer + Ring Newline Buffer Pass

Changes:

- multithreaded chunk boundaries are now produced on the fly by the producer thread (no upfront `buildChunks()` prepass in worker path)
- `BatchedLineScanner` newline offsets now use a circular ring buffer (removed compaction memcopy)
- strict loop switched to larger refill targets (`fillTo(1024)`) and inner batched record processing

Latest 1-run results (no warmup, HG002 R1, threads=14):

| Tool | Time (s) | GiB/s |
|---|---:|---:|
| zdash-strict-full | 26.210 | 0.653 |
| zdash-assume-full | 13.465 | 1.272 |
| zdash-assume-stats | 11.810 | 1.450 |
| seqtk | 17.698 | 0.967 |
| seqkit | 6.081 | 2.816 |
| fastp | 18.638 | 0.919 |
| noodles-rust | 20.756 | 0.825 |

Upshot:

- strict mode moved up substantially and now beats `seqtk`, `fastp`, and `noodles-rust` in this run
- `seqkit` still leads on this benchmark profile

## Command UX + Retest (check/scan/stats)

Changes:

- added user-facing commands:
  - `zdash check` (strict validation + stats)
  - `zdash scan` (fast assume-valid scan)
  - `zdash stats` (fast summary-focused mode)

Latest 1-run results (no warmup, HG002 R1, threads=14):

| Tool | Time (s) | GiB/s |
|---|---:|---:|
| zdash-check | 30.673 | 0.558 |
| zdash-scan | 11.733 | 1.459 |
| zdash-stats | 11.681 | 1.466 |
| seqtk | 17.545 | 0.976 |
| seqkit | 5.729 | 2.989 |
| fastp | 18.790 | 0.911 |
| noodles-rust | 20.913 | 0.819 |

Upshot:

- strict `check` now clearly beats `fastp` and `noodles-rust` in this run
- `scan/stats` are materially faster than strict and beat `seqtk`
- `seqkit` remains the speed leader on pure stats

Quick regression retest after gzip/stdin plumbing changes (HG002 R1, 1 run, no warmup, threads=14):

| Tool | Time (s) | GiB/s |
|---|---:|---:|
| zdash-check | 28.626 | 0.598 |
| zdash-scan | 13.221 | 1.295 |
| zdash-stats | 12.441 | 1.376 |
| seqkit | 6.276 | 2.728 |

Quick spot-check after `--max-errors` / context-artifact additions (HG002 R1, 1 run, threads=14):

| Tool | Time (s) | GiB/s |
|---|---:|---:|
| zdash-check | 23.753 | 0.721 |
| seqkit | 6.178 | 2.771 |
