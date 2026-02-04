# Z-DASH Operations Guide

Quick production usage patterns.

## Common Commands

- Gate pipeline with strict validation:
  - `zdash check --json --report-json zdash_report.json input.fastq`
- Fast health scan:
  - `zdash scan --json input.fastq.gz`
- Stats-only run:
  - `zdash stats --json input.fastq`

## Failure Handling

- Collect first 10 failures:
  - `zdash check --max-errors 10 input.fastq`
- Write failure context window:
  - `zdash check --extract-error-context=error_window.fastq input.fastq`
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
