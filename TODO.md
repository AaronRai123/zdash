# Z-DASH Project TODO

Living checklist for building a high-speed FASTQ sanitizer in Zig.

## Project Goal

Build a single-binary CLI (`zdash`) that validates FASTQ, computes summary stats, and flags low-quality reads in one pass, with performance competitive with `seqtk`, `fastp`, and `seqkit`.

## Milestone 0 - Repo Setup (Done / In Progress)

- [x] Initialize Zig project (`zig init`)
- [x] Confirm local Zig toolchain works (`zig build`)
- [x] Rename binary target from `ZDash` to `zdash`
- [x] Add `README.md` with scope, usage, and performance goals
- [x] Add `LICENSE`
- [ ] Keep `documentation.md` updated after every code change (ongoing)

## Milestone 1 - Correctness-First Sanitizer (v0.1)

### CLI + UX
- [x] Parse CLI args: `zdash <input.fastq>`
- [x] Add `--help` and `--version`
- [x] Define exit codes (`0` pass, non-zero corruption/error)
- [x] Print summary block with stable formatting

### FASTQ Parsing + Validation
- [x] Implement strict 4-line FASTQ record parser
- [x] Validate record line prefixes (`@` header, `+` separator)
- [x] Validate sequence/quality length match per read
- [x] Validate sequence characters (`A,C,G,T,N` + lowercase policy)
- [x] Validate quality ASCII range (Phred+33)
- [x] Handle last-line/truncated-file corruption detection
- [x] Report first error with byte offset and read index
### Statistics
- [x] Track total reads
- [x] Track total bases
- [x] Track GC base count
- [x] Track summed quality scores
- [x] Compute GC percent and average quality score
- [x] Define low-quality read rule (default threshold)
- [x] Track low-quality read count

### Testing (Core)
- [x] Unit tests for parser on valid/invalid mini FASTQ fixtures
- [x] Unit tests for GC and quality math
- [x] Golden test for final summary output format

## Milestone 2 - I/O and Memory Model (v0.2)

### mmap + Chunking
- [x] Add file-backed `mmap` path for large inputs
- [x] Keep fallback buffered-read mode for unsupported environments
- [x] Implement chunk partitioning aligned to FASTQ record boundaries
- [x] Ensure zero-copy slices into mapped memory (no per-read alloc)

### Robustness
- [x] Validate behavior on tiny, huge, and malformed files
- [x] Add file size and processed byte counters
- [x] Add optional progress display for long runs

### Testing
- [x] Integration tests across both read modes
- [x] Boundary tests for chunk splits (record straddling edges)

## Milestone 3 - Parallel Engine (v0.3)

### Threading Architecture
- [x] Implement worker pool sized to CPU cores (configurable)
- [x] Give each worker local counters (avoid contention in hot path)
- [x] Merge worker counters at end (single reduction step)
- [x] Add atomics only where truly required

### Correctness Under Concurrency
- [x] Ensure deterministic final totals across runs
- [x] Add tests for race-prone scenarios (small chunks, many threads)
- [x] Add thread-safety assertions in debug mode

### CLI + Ops Additions (new)
- [x] Add `--threads <N>` to configure worker count
- [x] Add `--chunk-bytes <N>` for tuning chunk size without recompiling
- [ ] Add explicit startup config line in output (threads/chunk/io mode)
- [ ] Add stable error code table to docs for pipeline integration

## Milestone 4 - SIMD Kernels (v0.4)

### Sequence Kernels
- [x] Implement SIMD counting for `G/g` and `C/c`
- [x] Implement SIMD validation mask for allowed base characters
- [x] Add scalar fallback for tail bytes and unsupported targets

### Quality Kernels
- [x] Implement SIMD accumulation of quality score bytes
- [x] Preserve exactness vs scalar reference implementation

### Verification
- [x] Bit-for-bit parity tests: SIMD vs scalar on random datasets
- [x] Benchmark kernels in isolation (microbench harness)

## Milestone 5 - Benchmarking and Competitive Goals (v0.5)

### Benchmark Harness
- [x] Add `--bench` mode (elapsed time, GB/s, reads/s)
- [x] Record CPU model, core count, storage type, Zig version
- [x] Stabilize benchmark methodology (warmup + repeated runs)
- [x] Add `--quiet` mode to suppress summary I/O during benchmarking
- [x] Remove extra chunk pre-pass (`countLines`) to avoid benchmarking overhead
- [x] Add apples-to-apples benchmark profile definitions per tool (what work is counted)
- [ ] Add reproducible benchmark script (`scripts/bench_competitors.sh` + CSV output)
- [ ] Benchmark on larger realistic datasets (>=10GB) to reduce startup noise
- [ ] Add optional mmap advisory hints (`madvise`) and test impact
- [ ] Add privileged hardware-counter profile runbook (IPC, cycles, instructions)
- [ ] Split framing vs kernel time inside `validateFastqChunk` with finer internal counters
- [x] Prototype SIMD newline scanner (`\n` index kernel) to reduce framing overhead
- [ ] Add CLI A/B switch to force parser mode (speculative fast-path on/off) for benchmarking
- [ ] Add thread/chunk autotuning helper (`--autotune`) and persist best local config
- [ ] Add scheduler A/B switch (`--schedule static|dynamic`) to preserve benchmark comparability
- [x] Add fused hot-path helpers that skip profiler wrapper overhead in non-profile runs
- [x] Add parser mode switch (`strict|assume-valid`) and fused assume-valid fast path
- [x] Add producer-consumer chunk queue path for assume-valid mode
- [x] Rewrite strict hot loop to operate on newline offsets (tail-only line fallback)
- [x] Apply ASCII case-fold sequence validation in strict SIMD/scalar kernels
- [x] Remove multithread upfront chunk prepass by generating chunk boundaries in producer
- [x] Replace newline offset compaction with circular ring buffer
- [ ] Rename internal profiler counters (`readLine` -> `newlineRefill`) to avoid ambiguity
- [ ] Add dedicated real-dataset benchmark script mode for thread sweep (`8/10/12`) with CSV output

### Competitive Targets
- [ ] Beat `seqtk` by 3-5x on multi-core runs
- [ ] Beat `fastp` by ~2x on representative workloads
- [ ] Beat `seqkit` by ~10x on representative workloads
- [ ] Publish benchmark tables and repro commands in `README.md`

## Milestone 6 - Production Hardening (v1.0)

- [ ] Add clear error messages and troubleshooting hints
- [ ] Add logging verbosity levels (`--quiet`, `--verbose`)
- [ ] Add configurable quality threshold and validation strictness
- [x] Add `--max-errors N` and multi-error reporting (developer debugging economics)
- [x] Add `--extract-error-context` to write failing reads +/- window
- [x] Add `--json` output for CI/pipeline automation ergonomics
- [x] Add `--report-json <path>` artifact output for CI
- [x] Add `repair` command (`drop-bad-records`, `truncate-to-last-good`)
- [x] Add deterministic `sample` command (`--seed`, `--fraction` or `--n`)
- [x] Add compressed input support strategy (`.gz`) or explicit non-goal
- [x] Add stdin input support (`-`) for streaming workflows
- [x] Add CI (build + tests on Linux/macOS)
- [x] Add benchmark guardrail script for check/scan/stats performance regression spotting
- [x] Add operations runbook doc (`OPERATIONS.md`)
- [ ] Freeze output format for downstream tooling compatibility
- [ ] Move tests into dedicated test files when parser grows (improves maintainability)
- [ ] Add regression FASTQ fixture corpus (`tests/fixtures/`)
- [ ] Add sanitizer/fuzz pass for parser edge cases

## Backlog (Post-v1)

- [ ] Paired-end FASTQ support
- [ ] Optional JSON output for pipelines
- [ ] Read-level filtering output file(s)
- [ ] Adapter/primer trimming modes
- [ ] ARM/x86 target-specific kernel tuning
- [ ] Optional JSON summary output (`--json`) for workflow engines
- [ ] Add provenance block (command, version, timestamp) to output/header

## Tracking Conventions

- Keep this file updated in each PR.
- Use checkboxes for completion.
- When closing a major item, add a short note in commit message.
- If scope changes, add a dated note below.

## Scope Change Log

- 2026-02-04: Initial roadmap created.
- 2026-02-04: Milestone 0 completed (binary rename, README, MIT license).
- 2026-02-04: Milestone 2 completed (mmap/buffered modes, chunking, progress, integration tests).
- 2026-02-04: Milestone 3 completed (threaded chunk workers, deterministic reductions, concurrency tests).
- 2026-02-04: Milestone 4 completed (SIMD sequence/quality kernels, parity tests, microbench mode).
- 2026-02-04: Milestone 5 harness implemented (`--bench`, metadata capture, warmup+repeat).
- 2026-02-04: Perf pass added `--quiet`, dynamic scheduling, profile modes, and new competitive benchmark results.
- 2026-02-04: Replaced line-by-line chunk parser with a single-pass parser state machine and re-profiled on real HG002 data.
- 2026-02-04: Evaluated parser path variants and thread/chunk tuning on HG002; best local run in this pass used `--threads 12`.
- 2026-02-04: Added SIMD-batched newline offset scanner (`BatchedLineScanner`) and re-ran HG002 timing/profile checks.
- 2026-02-04: Added speculative 4-line record fast-path with strict fallback; best HG002 local timing improved to ~43.19s (`--threads 12`).
- 2026-02-04: Switched to static contiguous chunk scheduling; HG002 repeated runs improved into ~37-40s range with `--threads 14`.
- 2026-02-04: Added fused hot-path kernel helpers + prefetch hint; HG002 repeated runs improved further to ~34.77-35.59s (`--threads 14`).
- 2026-02-04: Implemented Milestone 6B (`--mode assume-valid` + producer queue); assume-valid reached ~2.11-2.19 GiB/s on HG002 in 1-run benchmark.
- 2026-02-04: Applied strict-loop newline-offset rewrite + case-fold base validation; strict mode improved to ~0.51 GiB/s in latest HG002 1-run benchmark.
- 2026-02-04: Added streaming chunk producer (no upfront chunk list in multithread path) + ring newline buffer; strict mode improved to ~0.65 GiB/s in latest HG002 1-run benchmark.
- 2026-02-04: Added user-facing `check|scan|stats` commands and detailed validation diagnostics to improve developer economics.
- 2026-02-04: Added gzip/bgzip input support and stdin (`-`) ingestion; updated CLI docs and JSON workflow examples.
- 2026-02-04: Added `--max-errors` multi-error collection and `--extract-error-context` FASTQ artifact output for debugging workflows.
- 2026-02-04: Added `repair` and deterministic `sample` commands plus `--report-json` CI artifact output.
- 2026-02-04: Added JSON schema version contract (`JSON_SCHEMA.md`), integrated gzip decode path, and Linux/macOS GitHub Actions CI workflow.
- 2026-02-04: Added regression fixtures + repair golden tests, benchmark guardrail script, and operations guide.
