# Z-DASH Code Documentation

This document explains the current implementation of `zdash` in detail.

## 1. Current Scope

The code currently implements:

- CLI argument handling
- user-facing command modes (`check`, `scan`, `stats`)
- additional command modes (`repair`, `sample`)
- FASTQ input loading with selectable I/O mode (`auto`, `mmap`, `buffered`)
- gzip/bgzip input decompression support via integrated `std.compress.flate` decode (default in-memory stream mode, optional temp-file mode)
- stdin input support (`-`)
- Runtime tuning flags (`--threads`, `--chunk-bytes`)
- Work profile selection (`--profile`)
- Parser mode selection (`--mode strict|assume-valid`)
- Full-pipeline benchmark mode (`--bench`)
- FASTQ validation (strict 4-line records)
- Summary statistics
- Chunked processing with FASTQ-aligned chunk boundaries
- Multithreaded worker execution over chunks
- SIMD kernels for sequence/quality hot loops with scalar fallback
- Optional progress output
- Kernel microbenchmark mode (`--bench-kernels`)
- Unit and integration tests in the same source file

Primary implementation file: `src/main.zig`


## 2. High-Level Flow

Entry point: `main()`

1. Parse CLI args.
2. Handle `--help` / `--version`.
3. Parse run options:
   - input path
   - I/O mode preference
   - thread count (`--threads`, default CPU count)
   - chunk size (`--chunk-bytes`)
   - progress flag
   - optional benchmark mode (`--bench-kernels`)
   - optional full benchmark mode (`--bench`)
4. Open input file.
5. Load file bytes using selected/auto I/O mode.
6. Validate + compute stats over chunked FASTQ processing.
7. Render final summary and print it.
8. Exit with process code:
   - `0` success
   - `2` usage error
   - `3` I/O error
   - `4` validation failure

## 3. Core Data Types

### `ExitCode`

Maps semantic outcomes to process exit codes.

### `IoPreference`

User-selected preference from CLI:

- `auto`
- `mmap`
- `buffered`

### `RunOptions`

Parsed execution options:

- `input_path`
- `io_preference`
- `show_progress`
- `quiet`
- `threads` (`0` means auto-detect)
- `chunk_bytes`
- `profile`
- `bench_kernels`
- `bench`

### `ProcessingConfig`

Resolved runtime processing settings:

- `chunk_bytes`
- `threads`
- `profile`

### `InputMode`

Actual mode used after loading:

- `mmap`
- `buffered`

### `LoadedInput`

Represents loaded FASTQ content and ownership model.

Fields:

- `bytes`: read-only byte slice used by parser
- `mode`: chosen `InputMode`
- `file_size`: size from `stat`
- `backing`:
  - mapped memory slice
  - owned heap buffer

`deinit()` releases mapped memory or frees heap allocation.

### `ValidationIssue`, `ValidationError`, `ValidationResult`

Validation errors include:

- read index
- byte offset
- categorized issue

`ValidationResult` is either:

- `.ok` with `Stats`
- `.fail` with detailed `ValidationError`

### `Stats`

Accumulated metrics:

- total reads
- total bases
- GC base count
- quality score sum (Phred+33 decoded)
- low-quality read count

### `Line`

Helper type from `readLine()` containing:

- line bytes
- starting byte offset in the active slice

### Parallel Processing Types

- `Chunk`: byte start/end + read-index base for each chunk
- `WorkerContext`: worker assignment (chunk index range) + shared inputs
- `WorkerResult`: worker-local stats + worker-local earliest error
- `ProgressState`: mutex-protected progress accounting for thread-safe printing

## 4. CLI Behavior

Supported options:

- `-h`, `--help`
- `-V`, `--version`
- `--progress`
- `--quiet`
- `--threads <N>`
- `--threads=<N>`
- `--chunk-bytes <N>`
- `--chunk-bytes=<N>`
- `--profile <full|validate-stats|stats-only>`
- `--profile=<full|validate-stats|stats-only>`
- `--mode <strict|assume-valid>`
- `--mode=<strict|assume-valid>`
- `--max-errors <N>`
- `--max-errors=<N>`
- `--error-window-reads <N>`
- `--error-window-reads=<N>`
- `--extract-error-context`
- `--extract-error-context=<path>`
- `--extract-error-debug <path>`
- `--extract-error-debug=<path>`
- `--report-json <path>`
- `--report-json=<path>`
- `--against <path>`
- `--against=<path>`
- `--output <path|->`
- `--output=<path|->`
- `--repair-mode <drop-bad-records|truncate-to-last-good>`
- `--repair-mode=<drop-bad-records|truncate-to-last-good>`
- `--emit-bad-records <path>`
- `--emit-bad-records=<path>`
- `--seed <N>`, `--seed=<N>`
- `--fraction <F>`, `--fraction=<F>`
- `--n <N>`, `--n=<N>`
- `--io-mode <auto|mmap|buffered>`
- `--io-mode=<auto|mmap|buffered>`
- `--gzip-mode <stream|temp>`
- `--gzip-mode=<stream|temp>`
- `--bench-kernels`
- `--bench`
- `--profile-internal`
- `--json`
- `--json-schema-version <V>`, `--json-schema-version=<V>`
- `--gha-annotations`
- `--preset <strict-ci|fast-scan|qc-only>`, `--preset=<...>`
- `--config <path>`, `--config=<path>`
- `--ci`

Commands:

- `check`: strict validation + stats (trust-first default)
- `scan`: assume-valid fast scan + stats
- `stats`: fastest summary mode (assume-valid + stats profile)
- `repair`: rewrite FASTQ with configurable bad-record handling
- `sample`: deterministic FASTQ subsampling
- `explain`: interpret JSON report failures with actionable hints
- `compare`: diff key metrics/status fields across two reports (`--against`)
- `doctor`: environment diagnostics and recommended runtime defaults

Input path behavior:

- regular path: direct load (`mmap` preferred)
- `.gz` / `.bgz`: default in-memory stream decode then processed (`--gzip-mode temp` keeps temp-file bridge)
- `-`: reads FASTQ from stdin into buffered memory

Input path is required and must be exactly one positional argument, except `doctor` (input optional).

Unknown flags and malformed options return usage error (`2`).

`--quiet` behavior:

- suppresses the main human-readable summary block
- still prints benchmark block when `--bench` is used
- still prints a brief validation failure line to stderr on corruption
- in failure cases, stderr now includes structured diagnostics (read, line type, byte offset, bad char, caret context)

`--json` behavior:

- emits one JSON document with status, config, and metrics/error fields
- intended for CI and pipeline integration
- includes `schema_version` for compatibility checks
- `--json-schema-version <V>` enforces exact schema-version matching and fails fast on mismatch
- for `sample`, `--json` emits sampling summary JSON and requires `--output <path>` (cannot use `-`)
- for validation operations, failing runs auto-write `zdash_report.json` when `--json` is used and `--report-json` is not provided

`--ci` behavior:

- alias for CI-safe defaults (`strict-ci` preset + JSON schema pin + GitHub annotations + report path)

`--gha-annotations` behavior:

- on validation failures, emits GitHub Actions `::error ...` annotations with read/byte context

Preset behavior:

- `strict-ci`: `check` + strict mode + full profile + JSON enabled
- `fast-scan`: `scan` + assume-valid mode + `validate-stats` profile
- `qc-only`: `stats` + assume-valid mode + `stats-only` profile

Config behavior:

- `--config <path>` loads startup defaults from simple `key=value` lines
- CLI flags still override config values

Error collection behavior:

- default is fail-fast (`--max-errors 1`)
- `--max-errors N` collects up to `N` validation failures before returning
- `--extract-error-context[=<path>]` writes a FASTQ window around the first failure for debugging
- `--error-window-reads <N>` controls context FASTQ window size
- `--extract-error-debug <path>` writes a human-readable debug report with metadata and a local byte hex window
- `--report-json <path>` writes the same summary/error JSON payload to file

Repair behavior:

- `repair --repair-mode drop-bad-records`: writes only valid reads to output
- `repair --repair-mode truncate-to-last-good`: stops output at first bad read
- `--emit-bad-records <path>` writes rejected reads for debugging/recovery

Sample behavior:

- `sample` requires exactly one of `--fraction` or `--n`
- sampling is deterministic with `--seed`
- output defaults to stdout unless `--output` is set
- with `--json`, `sample` requires `--output <path>` so JSON output is not mixed with FASTQ bytes

`--profile-internal` behavior:

- prints internal aggregated timing counters for core functions
- reports approximate wall-time sums across threads (useful for hotspot ranking)
- does not provide hardware counters (IPC/cycles/instructions)
- includes framing + scheduling counters (`strictFraming`, `assumeFraming`, `newlineRefill`, `queueWait`) for optimization triage

### Exit codes

- `0`: success
- `2`: usage/argument error
- `3`: I/O or processing setup failure
- `4`: FASTQ validation failure (corruption)

## 5. I/O Model

Implemented in `loadInputFromPath(...)` and `loadInputData(...)`.

### Gzip/bgzip mode

- `--gzip-mode stream` (default): decode `.gz`/`.bgz` directly into in-memory FASTQ bytes.
- `--gzip-mode temp`: decode into a temporary FASTQ file, then use existing file-loading path.
- stream mode avoids temp-file write/read overhead and is usually much faster on large compressed inputs.

### Auto mode

- Tries memory mapping for regular non-empty files on non-Windows.
- Falls back to buffered full-file read if mmap fails.

### Forced mmap mode

- Requires mmap support and successful mapping.
- Returns error instead of falling back.

### Forced buffered mode

- Always reads file into owned memory.

### Buffered reader details

- Uses `seekTo(0)` before read.
- Uses file size as max bound when available.

## 6. FASTQ Parsing and Validation

Parsing model is strict FASTQ 4-line records:

1. header (must start with `@`)
2. sequence (only `A,C,G,T,N` plus lowercase)
3. separator (must start with `+`)
4. quality (ASCII range `33..126`)

Additional checks:

- sequence length must equal quality length
- truncated/incomplete records are detected

On first failure, processing stops and records:

- read index (1-based)
- absolute byte offset
- issue type

Current framing implementation note:

- chunk parsing uses `BatchedLineScanner`, which pre-scans newline offsets in SIMD-sized blocks and reuses them
- each parsed line is consumed directly as a slice (no per-line allocation)
- validation and stats still run in the same record pass
- `validateFastqChunk` now includes a speculative 4-line fast path that bulk-validates `@/+` framing and equal seq/qual lengths before kernel dispatch
- malformed or unusual records fall back to the existing strict line-by-line safety path for precise errors
- hot path now uses `analyzeSequenceHot` / `analyzeQualityHot` to bypass profiler wrappers in normal runs
- hot path also issues a lightweight `@prefetch` for upcoming record starts when enough newline offsets are buffered
- strict mode now runs primarily on newline-offset blocks (`fillTo(4)` + `peekNewline`) and only uses line-sliced parsing at chunk tails
- strict and assume-valid paths now use profile-specialized chunk loops to remove runtime validation/low-quality branch checks in hot code

## 7. Chunked Processing

`validateFastq(...)` delegates to `validateFastqWithChunkSizeAndThreads(...)`.

Current target chunk size:

- `TARGET_CHUNK_BYTES = 8 MiB`

Chunking rules:

- each chunk end is aligned to FASTQ record boundary
- boundary means line count modulo 4 equals zero
- if no boundary found past target, parser uses end-of-file

Per-chunk stats are merged into global totals.

Chunk metadata also tracks `read_index_base`, enabling deterministic absolute read index reporting across threaded workers.

## 8. Parallel Execution Model

Worker count:

- defaults to CPU count (bounded)
- configurable with `--threads`
- clamped to chunk count

Chunk sizing:

- configurable with `--chunk-bytes`
- defaults to `8 MiB`

Execution strategy:

- build chunk list once
- build chunk metadata in a single pass over file bytes (no per-chunk recount pass)
- workers process static contiguous chunk ranges (one contiguous region per worker)
- each worker processes only its assigned range
- each worker writes only to its own `WorkerResult` (no hot-path shared counters)
- main thread merges results in stable worker order

Assume-valid mode threading:

- uses a producer-consumer path with a lock-free chunk queue
- producer thread prefetches upcoming chunk starts and enqueues chunk indices
- workers pop chunk indices and process independently

Chunk boundary production:

- chunk boundaries are now generated by the producer during processing (single scan)
- removed full upfront chunk-list prepass from the multithreaded path

Error handling under concurrency:

- workers keep their local earliest error
- reducer picks global earliest error by byte offset, then read index
- this makes failure selection deterministic

### Profiling Snapshot (Real HG002, Feb 4, 2026)

Command:

- `./zig-out/bin/zdash --quiet --profile-internal --threads 10 data/hg002/SRR26901703_1.fastq`

Top internal hotspots (aggregated wall-time approximation):

- `validateFastqChunk`: ~62.1%
- `readLine` counter (newline batch refill): ~36.3%
- `analyzeSequenceKernel`: ~1.0%

Implication:

- framing and per-record control flow inside `validateFastqChunk` are still dominant
- SIMD math kernels are no longer the bottleneck on this dataset

Tuning note:

- on the same HG002 file, `--threads 12` materially outperformed `--threads 10` in local runs
- after batched newline scanning + speculative fast path, local HG002 runtime reached ~43.19s at `--threads 12`
- static contiguous scheduling further improved repeated HG002 runs into the ~37-40s range at `--threads 14`
- latest hot-path pass improved repeated HG002 runs to ~34.77-35.59s at `--threads 14`

Thread-safety checks:

- shared progress output is mutex-protected
- chunk assignment is lock-free via atomic index fetch

## 9. Statistics Computation

Computed during validation pass:

- total reads: increment per valid header line encountered
- total bases: sum of sequence lengths
- GC bases: count `G/g/C/c`
- quality sum: sum over `(quality_char - 33)`
- low-quality reads: read average Q below threshold

Current threshold:

- `LOW_QUALITY_AVG_Q_THRESHOLD = 20.0`

Profile behavior:

- `full`: strict validation + stats + low-quality counting
- `validate-stats`: strict validation + stats (no low-quality counting)
- `stats-only`: structural parsing + stats, skips base/quality character validation

Mode behavior:

- `strict`: full strict parser behavior with precise corruption reporting
- `assume-valid`: structural framing checks (`@`, `+`, seq/qual length) with fused seq+qual accumulation and no per-char validation

Strict validation kernels:

- sequence validation now uses ASCII case-folding (`ch & 0xDF`) to reduce SIMD compare count
- strict fast path validates blocks first and only runs scalar pinpoint scans when anomalies are detected
- strict parser now processes newline offsets in larger batches (`fillTo(1024)`) with an inner record loop
- newline offset storage uses a ring buffer (no compaction memmove in refill path)

Derived metrics:

- average quality = `quality_sum / total_bases`
- GC percent = `gc_bases * 100 / total_bases`
- low-quality read decision uses integer compare (`read_q_sum < threshold * read_len`)

## 10. SIMD Kernels

SIMD selection:

- uses `std.simd.suggestVectorLength(u8)` to pick width
- prefers 32-byte lanes when available, otherwise 16-byte lanes
- falls back to scalar on unsupported targets

Sequence kernel:

- SIMD-validates allowed bases via vector mask
- SIMD-counts `G/g/C/c` lanes
- scalar fallback for tails and invalid-byte pinpointing
- scalar path uses precomputed lookup tables for base validity and GC flags

Quality kernel:

- SIMD-validates ASCII `33..126`
- SIMD-accumulates `(q - 33)` score sums
- scalar fallback for tails and invalid-byte pinpointing
- scalar path uses precomputed lookup tables for quality validity and decoded scores

Parity guarantees:

- random-data tests compare scalar and SIMD outputs exactly
- behavior matches scalar semantics on early invalid bytes

## 11. Progress Output

Enabled with `--progress`.

Behavior:

- prints progress to stderr
- updates per completed chunk (from any worker)
- shows processed bytes and percentage
- emits final newline at 100%
- uses mutex-protected progress state for thread-safe output

## 12. Rendering and Output

Summary rendering is separated for testability:

- `renderSummary(...)` returns formatted string
- `printSummary(...)` writes rendered text to stdout

Success output includes:

- metrics
- low-quality count
- selected I/O mode
- selected profile
- thread count
- chunk bytes
- file bytes
- processed bytes
- status

Failure output includes:

- selected I/O mode
- file bytes
- first error details

When `--bench` is enabled, an additional benchmark block is printed with:

- warmup/measurement counts
- avg/best time
- avg/best throughput (GiB/s)
- reads/s
- host metadata (core count, arch, OS, storage type, Zig version)
- active config (threads, chunk bytes, I/O mode)
- profile

Benchmark tip:

- use `--quiet` during benchmark runs to reduce output overhead

## 13. Tests

Tests are currently in `src/main.zig`.

Covered:

- parser + stats correctness on valid fixture
- invalid sequence character
- sequence/quality mismatch
- truncated record
- golden summary output formatting
- chunk boundary correctness
- tiny chunk processing correctness
- mmap vs buffered integration parity
- robustness scenarios (tiny, large synthetic, malformed)
- deterministic totals across repeated parallel runs
- race-prone scenario (many threads + tiny chunks)
- SIMD sequence parity vs scalar (random datasets)
- SIMD quality parity vs scalar (random datasets)
- option parsing for `--quiet` and `--profile`
- option parsing guards for `sample --json` output-path requirements
- profile semantics (`validate-stats` and `stats-only`)

Run tests with:

```bash
zig build test
```

Run kernel microbench:

```bash
zig-out/bin/zdash --bench-kernels
```

Run full pipeline benchmark:

```bash
zig-out/bin/zdash --bench --threads=8 --chunk-bytes=8388608 input.fastq
```

## 14. Known Limitations / Next Engineering Steps

Current limitations:

- reads full file into memory in buffered mode
- gzip/bgzip in-memory decode currently buffers full decompressed payload (bounded-memory streaming still pending)
- parser and tests still co-located in `src/main.zig` (needs modular split)
- JSON schema is versioned, but formal change-log discipline must be maintained during releases

Recommended next direction:

- add benchmark guardrail thresholds per mode (`check/scan/stats`) in CI
- add bounded-memory guarantees for stdin and gzip paths
- split parser/validator/commands into dedicated source files

## 15. CI + Fixtures

- CI workflow: `.github/workflows/ci.yml` (Linux + macOS)
- fixtures live under `tests/fixtures`
- repair golden tests validate `drop-bad-records` and `truncate-to-last-good` behavior

## 16. Maintenance Rule

Project rule requested by owner:

- Update `documentation.md` after every meaningful code change so docs stay synchronized with implementation.
