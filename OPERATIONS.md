# Z-DASH Operations Guide

Quick production usage patterns.

## Common Commands

- Gate pipeline with strict validation:
  - `zdash check --json --json-schema-version 1.0.0 --report-json zdash_report.json input.fastq`
  - `zdash check --json --gha-annotations input.fastq` (for GitHub Actions logs)
  - `zdash --ci input.fastq` (single-flag CI preset)
- Fast health scan:
  - `zdash scan --json input.fastq.gz`
  - `zdash scan --json --gzip-mode stream input.fastq.gz`
- Stats-only run:
  - `zdash stats --json input.fastq`

## Failure Handling

- Collect first 10 failures:
  - `zdash check --max-errors 10 input.fastq`
- Write failure context window:
  - `zdash check --extract-error-context=error_window.fastq input.fastq`
- Tune context window size and emit a human debug report:
  - `zdash check --extract-error-context=error_window.fastq --error-window-reads 20 --extract-error-debug error_debug.txt input.fastq`
- Exit code `4` means validation failure; stop downstream pipeline.

## Repair Workflows

- Drop broken records:
  - `zdash repair --repair-mode drop-bad-records --output repaired.fastq input.fastq`
- Truncate at first bad record:
  - `zdash repair --repair-mode truncate-to-last-good --output repaired.fastq input.fastq`
- Keep rejected records:
  - `--emit-bad-records rejected.fastq`

## Sampling Workflows

- Fractional deterministic sample:
  - `zdash sample --seed 42 --fraction 0.01 --output sample.fastq input.fastq`
- Fixed-size deterministic sample:
  - `zdash sample --seed 42 --n 100000 --output sample.fastq input.fastq`
- JSON sample report:
  - `zdash sample --json --seed 42 --n 100000 --output sample.fastq input.fastq`

## Input Types

- FASTQ file: `input.fastq`
- gzip/bgzip: `input.fastq.gz`, `input.fastq.bgz`
- stdin: `-` (example: `cat input.fastq | zdash stats --json -`)
- gzip mode:
  - `--gzip-mode stream` (default, in-memory decode)
  - `--gzip-mode temp` (legacy temp-file bridge)
  - expect `stream` to be substantially faster on large `.gz` inputs since it avoids temp-file write+read.

## Presets + Config

- Presets:
  - `--preset strict-ci` (check + strict + full + json)
  - `--preset fast-scan` (scan + assume-valid + validate-stats)
  - `--preset qc-only` (stats + assume-valid + stats-only)
- Config file:
  - `zdash --config .zdash.example.toml input.fastq`

## Doctor

- Environment health check:
  - `zdash doctor`
  - `zdash doctor --json`

## Report Tooling

- Explain failure reports:
  - `zdash explain zdash_report.json`
- Compare before/after report metrics:
  - `zdash compare --against repaired_report.json original_report.json`

## Speed + Accuracy Validation

- Quick performance + correctness pass:
  - `WARMUP_RUNS=1 MEASURED_RUNS=3 scripts/perf_validate.sh data/hg002/SRR26901703_1.fastq`
- Output includes p50/p90 seconds and GiB/s for `check/scan/stats`, plus fixture-based accuracy checks.
