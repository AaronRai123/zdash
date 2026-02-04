# Z-DASH

Trust-first FASTQ tool written in Zig.

Z-DASH validates FASTQ integrity, reports actionable corruption diagnostics, and computes summary stats with fast scan modes for large data.

## One Command Start

```bash
./run.sh
```

That builds `zdash` and runs CI-safe validation on `tests/fixtures/valid_small.fastq`.

Common variants:

```bash
./run.sh your.fastq
./run.sh your.fastq --max-errors 10
./run.sh doctor --json
```

## Install (No Zig Required)

Fastest install path is the release installer script:

```bash
curl -fsSL https://raw.githubusercontent.com/AaronRai123/zdash/main/scripts/install.sh | bash
```

After install:

```bash
zdash check file.fastq
cat file.fastq | zdash stats --json -
```

Homebrew tap flow (after publishing release checksums into tap formula):

```bash
brew tap AaronRai123/tap
brew install AaronRai123/tap/zdash
```

Formula template for the tap lives in `Formula/zdash.rb`.

## Status

Active development (pre-v1.0).  
Roadmap and task tracking live in `TODO.md`.

## Goals

- Fail fast on broken inputs with precise diagnostics (read index, byte offset, line type, bad char, caret context).
- Provide practical command modes:
  - `check`: strict correctness + stats
  - `scan`: faster assume-valid scan
  - `stats`: fastest summary mode
- Keep performance competitive via mmap, multicore execution, SIMD kernels, and producer/consumer chunking.

## CLI

```bash
zdash <check|scan|stats|repair|sample|explain|compare> [options] <input>

# strict validation + stats
zdash check sample.fastq

# fast scan mode
zdash scan sample.fastq

# summary-first mode
zdash stats sample.fastq

# machine-readable output for CI/pipelines
zdash check --json sample.fastq
zdash check --report-json zdash_report.json sample.fastq
zdash check --json --json-schema-version 1.0.0 sample.fastq
zdash check --gha-annotations sample.fastq
zdash --ci sample.fastq

# preset + config driven runs
zdash --preset strict-ci sample.fastq
zdash --config .zdash.example.toml sample.fastq

# gzip/bgzip input
zdash check sample.fastq.gz
zdash check --gzip-mode stream sample.fastq.gz
# use legacy temp-file bridge only if needed for compatibility
zdash check --gzip-mode temp sample.fastq.gz

# stream from stdin
cat sample.fastq | zdash stats --json -

# collect multiple validation failures
zdash check --max-errors 10 sample.fastq

# write failing read window for debugging
zdash check --extract-error-context=bad_window.fastq sample.fastq

# repair bad inputs
zdash repair --repair-mode drop-bad-records --output repaired.fastq sample.fastq

# deterministic sampling
zdash sample --seed 42 --fraction 0.01 --output sample_1pct.fastq sample.fastq
zdash sample --seed 42 --n 100000 --output sample_100k.fastq sample.fastq

# sample JSON summary (requires file output)
zdash sample --json --seed 42 --n 100000 --output sample_100k.fastq sample.fastq

# explain a JSON report with actionable hints
zdash explain zdash_report.json

# compare two JSON reports
zdash compare --against after_report.json before_report.json

# environment diagnostics
zdash doctor
zdash doctor --json
```

Validation errors are developer-friendly:

```text
ERROR: sequence contains invalid base (allowed: A,C,G,T,N, lowercase)
Read: 12483102
Line: sequence
Byte offset: 8192441203
Bad char: 'X' (0x58)
Context: ACGTXGTT
             ^
```

## Development Setup

Prerequisites:
- Zig `0.15.x` (tested with `0.15.2`)

Build:

```bash
zig build
```

Run:

```bash
zig build run -- check <input.fastq>
```

## Roadmap

Milestones are tracked in `TODO.md`:

- v0.1 correctness-first parser and summary stats
- v0.2 mmap and chunked zero-copy processing
- v0.3 multicore execution
- v0.4 SIMD kernels
- v0.5 competitive benchmarks vs `seqtk`, `fastp`, `seqkit`
- v1.0 hardening and CI

## Benchmarks

Latest benchmark notes and result tables are tracked in `benchmarks.md`.
Current competitive focus compares `zdash-check`, `zdash-scan`, and `zdash-stats` against `seqtk`, `seqkit`, `fastp`, and a `noodles` benchmark binary.
Use `scripts/perf_validate.sh` to run a quick throughput + accuracy validation pass (reports p50/p90 timings and throughput).

## JSON Contract

- JSON outputs (`--json`, `--report-json`) follow a versioned schema.
- See `JSON_SCHEMA.md` for compatibility guarantees and field definitions.
- `--json-schema-version <v>` enforces an exact schema version at runtime.
- If `--json` is used and validation fails, `zdash_report.json` is auto-written unless `--report-json` is provided.

## Data Source Support

- Uncompressed FASTQ: `*.fastq`
- Gzip/bgzip FASTQ: `*.gz`, `*.bgz` (default in-memory stream decode; `--gzip-mode temp` for legacy temp-file flow)
- Stdin stream: `-`

## CI

GitHub Actions CI runs on Linux and macOS:

- `zig build test`
- `zig build -Doptimize=ReleaseFast`
- smoke CLI checks for `check/scan/stats/sample/repair`

## Ops Docs

- `OPERATIONS.md` contains production command patterns and failure handling.
- `scripts/bench_guardrail.sh` provides a simple throughput guardrail check.
- `pitch.ms` is a concise product pitch with current capabilities and limitations.

## License

MIT. See `LICENSE`.
