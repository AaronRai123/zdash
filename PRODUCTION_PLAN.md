# Z-DASH Production Plan

Focused phases to move from current state to production-ready `v1.0`.

## Phase 1 - Reliability Core (Trust MVP)

Goal: make `zdash check` dependable for pipeline gatekeeping.

- [x] Add `--max-errors N` and aggregate error reporting.
- [x] Add `--extract-error-context` (write failing read +/- window).
- [ ] Freeze clear exit code contract and document it.
- [x] Add robust stdin/file error-path tests (including malformed/truncated edge cases).

## Phase 2 - Pipeline Economics

Goal: reduce developer/operator time in real workflows.

- [x] Add stable `--json` schema versioning and compatibility note.
- [x] Add `--report-json <path>` for artifact-friendly CI output.
- [x] Add deterministic `sample` command:
  - `sample --seed <N> --fraction <F>`
  - `sample --seed <N> --n <K>`
- [x] Add `check --fail-fast|--max-errors` ergonomics and examples.

## Phase 3 - Repair Workflows

Goal: convert breakage into actionable output.

- [x] Add `repair` command with modes:
  - `drop-bad-records`
  - `truncate-to-last-good`
  - `emit-bad-records.fastq`
- [x] Emit paired machine report (`report.json`) for repairs.
- [x] Add golden tests for each repair mode.

## Phase 4 - Compression + Streaming Hardening

Goal: production-grade ingestion path for real FASTQ feeds.

- [x] Replace shell `gzip -dc` bridge with integrated decode path.
- [ ] Add performance tests for `.fastq.gz` and stdin streaming.
- [ ] Add bounded-memory streaming guarantees (documented behavior).
- [ ] Add backpressure/throughput tuning knobs only if needed.

## Phase 5 - Release Hardening

Goal: ship `v1.0` with confidence.

- [x] CI matrix: macOS + Linux build/test.
- [ ] Add regression fixtures corpus (`tests/fixtures`).
- [x] Add benchmark guardrail script for key commands (`check/scan/stats`).
- [ ] Add versioned CLI/JSON changelog policy.
- [x] Write operator docs (quickstart, troubleshooting, expected failure modes).

## Exit Criteria for v1.0

- `check` trusted as a pipeline gate with actionable diagnostics.
- JSON output stable and documented.
- Repair workflows available and tested.
- `.fastq`, `.fastq.gz`, and stdin paths tested in CI.
- Clear release notes and compatibility guarantees.
