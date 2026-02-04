const std = @import("std");
const builtin = @import("builtin");

const VERSION = "0.1.0-dev";
const JSON_SCHEMA_VERSION = "1.0.0";
const LOW_QUALITY_AVG_Q_THRESHOLD: f64 = 20.0;
const TARGET_CHUNK_BYTES: usize = 8 * 1024 * 1024;
const MAX_THREADS: usize = 256;
const WORK_CHUNK_BATCH: usize = 8;
const NEWLINE_BATCH_CAP: usize = 1024;
const NEWLINE_BATCH_MASK: usize = NEWLINE_BATCH_CAP - 1;
const SIMD_SUGGESTED = std.simd.suggestVectorLength(u8) orelse 0;
const SIMD_WIDTH: usize = if (SIMD_SUGGESTED >= 32) 32 else if (SIMD_SUGGESTED >= 16) 16 else 0;
const SIMD_AVAILABLE = SIMD_WIDTH > 0;
const VALID_BASE_TABLE = buildValidBaseTable();
const GC_BASE_TABLE = buildGcBaseTable();
const VALID_QUALITY_TABLE = buildValidQualityTable();
const DECODED_QUALITY_TABLE = buildDecodedQualityTable();

comptime {
    if ((NEWLINE_BATCH_CAP & NEWLINE_BATCH_MASK) != 0) {
        @compileError("NEWLINE_BATCH_CAP must be a power of two");
    }
}

const ExitCode = enum(u8) {
    ok = 0,
    usage_error = 2,
    io_error = 3,
    validation_error = 4,
};

const InputMode = enum {
    mmap,
    buffered,
};

const IoPreference = enum {
    auto,
    mmap,
    buffered,
};

const GzipMode = enum {
    stream,
    temp,
};

const WorkProfile = enum {
    full,
    validate_stats,
    stats_only,
};

const Operation = enum {
    check,
    scan,
    stats,
    repair,
    sample,
    explain,
    compare,
};

const RunMode = enum {
    strict,
    assume_valid,
};

const RepairMode = enum {
    drop_bad_records,
    truncate_to_last_good,
};

const Preset = enum {
    strict_ci,
    fast_scan,
    qc_only,
};

const RunOptions = struct {
    operation: Operation = .check,
    operation_explicit: bool = false,
    input_path: []const u8,
    io_preference: IoPreference = .auto,
    gzip_mode: GzipMode = .stream,
    preset: ?Preset = null,
    config_path: ?[]const u8 = null,
    show_progress: bool = false,
    gha_annotations: bool = false,
    quiet: bool = false,
    json: bool = false,
    json_schema_version: ?[]const u8 = null,
    max_errors: usize = 1,
    error_window_reads: u64 = 5,
    extract_error_context_path: ?[]const u8 = null,
    extract_error_debug_path: ?[]const u8 = null,
    report_json_path: ?[]const u8 = null,
    compare_against_path: ?[]const u8 = null,
    output_path: ?[]const u8 = null,
    repair_mode: RepairMode = .drop_bad_records,
    emit_bad_records_path: ?[]const u8 = null,
    sample_seed: u64 = 42,
    sample_fraction: ?f64 = null,
    sample_n: ?usize = null,
    threads: usize = 0, // 0 means auto
    chunk_bytes: usize = TARGET_CHUNK_BYTES,
    profile: WorkProfile = .full,
    profile_explicit: bool = false,
    mode: RunMode = .strict,
    mode_explicit: bool = false,
    bench: bool = false,
    bench_kernels: bool = false,
    profile_internal: bool = false,
};

const ProcessingConfig = struct {
    chunk_bytes: usize,
    threads: usize,
    profile: WorkProfile,
    mode: RunMode,
};

const ValidationOutcome = struct {
    validation: ValidationResult,
    errors: []ValidationError = &.{},
    context_output_path: ?[]u8 = null,

    fn deinit(self: ValidationOutcome, allocator: std.mem.Allocator) void {
        if (self.errors.len > 0) allocator.free(self.errors);
        if (self.context_output_path) |p| allocator.free(p);
    }
};

const LoadedInput = struct {
    bytes: []const u8,
    mode: InputMode,
    file_size: u64,
    file_kind: std.fs.File.Kind,
    temp_path: ?[]u8 = null,
    backing: union(enum) {
        mapped: []align(std.heap.page_size_min) const u8,
        owned: []u8,
    },

    fn deinit(self: LoadedInput, allocator: std.mem.Allocator) void {
        switch (self.backing) {
            .mapped => |mapped| std.posix.munmap(mapped),
            .owned => |owned| allocator.free(owned),
        }
        if (self.temp_path) |temp_path| {
            std.fs.deleteFileAbsolute(temp_path) catch {};
            allocator.free(temp_path);
        }
    }
};

const ValidationIssue = union(enum) {
    truncated_record: struct { expected_line: u8 },
    invalid_header_prefix,
    invalid_separator_prefix,
    invalid_sequence_char: struct { ch: u8 },
    invalid_quality_char: struct { ch: u8 },
    sequence_quality_length_mismatch: struct { sequence_len: usize, quality_len: usize },
};

const ValidationError = struct {
    read_index: u64,
    byte_offset: usize,
    issue: ValidationIssue,
};

const ValidationResult = union(enum) {
    ok: Stats,
    fail: ValidationError,
};

const BenchResult = struct {
    warmup_runs: usize,
    measured_runs: usize,
    best_ns: u64,
    avg_ns: u64,
    validation: ValidationResult,
};

const Stats = struct {
    total_reads: u64 = 0,
    total_bases: u64 = 0,
    gc_bases: u64 = 0,
    quality_sum: u64 = 0,
    low_quality_reads: u64 = 0,
};

const Line = struct {
    bytes: []const u8,
    start: usize,
};

const BatchedLineScanner = struct {
    data: []const u8,
    line_start: usize = 0,
    search_pos: usize = 0,
    newline_offsets: [NEWLINE_BATCH_CAP]usize = undefined,
    newline_head: usize = 0,
    newline_tail: usize = 0,
    newline_len: usize = 0,

    fn init(data: []const u8) BatchedLineScanner {
        return .{
            .data = data,
            .line_start = 0,
            .search_pos = 0,
        };
    }

    fn available(self: *const BatchedLineScanner) usize {
        return self.newline_len;
    }

    fn fillTo(self: *BatchedLineScanner, min_available: usize) void {
        if (self.available() >= min_available) return;
        if (self.newline_len >= NEWLINE_BATCH_CAP) return;

        const _prof_start = profileStart();
        defer profileEnd(.newline_refill, _prof_start);

        var i = self.search_pos;
        if (SIMD_AVAILABLE and SIMD_WIDTH > 0 and self.newline_len < NEWLINE_BATCH_CAP) {
            const Vec = @Vector(SIMD_WIDTH, u8);
            const newline_vec: Vec = @splat(@as(u8, '\n'));
            while (i + SIMD_WIDTH <= self.data.len and self.newline_len < NEWLINE_BATCH_CAP) : (i += SIMD_WIDTH) {
                const vec: Vec = loadVector(SIMD_WIDTH, self.data, i);
                const mask = vec == newline_vec;
                if (!@reduce(.Or, mask)) continue;
                var lane: usize = 0;
                while (lane < SIMD_WIDTH) : (lane += 1) {
                    if (!mask[lane]) continue;
                    self.newline_offsets[self.newline_head] = i + lane;
                    self.newline_head = (self.newline_head + 1) & NEWLINE_BATCH_MASK;
                    self.newline_len += 1;
                    if (self.newline_len >= NEWLINE_BATCH_CAP) {
                        self.search_pos = i + lane + 1;
                        return;
                    }
                }
            }
        }

        while (i < self.data.len and self.newline_len < NEWLINE_BATCH_CAP) : (i += 1) {
            if (self.data[i] != '\n') continue;
            self.newline_offsets[self.newline_head] = i;
            self.newline_head = (self.newline_head + 1) & NEWLINE_BATCH_MASK;
            self.newline_len += 1;
        }
        self.search_pos = i;
    }

    fn nextLine(self: *BatchedLineScanner) ?Line {
        if (self.line_start >= self.data.len) return null;
        self.fillTo(1);

        if (self.newline_len > 0) {
            const raw_end = self.newline_offsets[self.newline_tail];
            self.newline_tail = (self.newline_tail + 1) & NEWLINE_BATCH_MASK;
            self.newline_len -= 1;

            const start = self.line_start;
            var line_end = raw_end;
            if (line_end > start and self.data[line_end - 1] == '\r') line_end -= 1;
            self.line_start = raw_end + 1;
            return .{
                .bytes = self.data[start..line_end],
                .start = start,
            };
        }

        const start = self.line_start;
        var line_end = self.data.len;
        if (line_end > start and self.data[line_end - 1] == '\r') line_end -= 1;
        self.line_start = self.data.len;
        return .{
            .bytes = self.data[start..line_end],
            .start = start,
        };
    }

    fn hasCompleteRecord(self: *BatchedLineScanner) bool {
        self.fillTo(4);
        return self.available() >= 4;
    }

    fn peekNewline(self: *const BatchedLineScanner, rel_index: usize) usize {
        return self.newline_offsets[(self.newline_tail + rel_index) & NEWLINE_BATCH_MASK];
    }

    fn consumeRecord(self: *BatchedLineScanner) void {
        const nl3 = self.newline_offsets[(self.newline_tail + 3) & NEWLINE_BATCH_MASK];
        self.newline_tail = (self.newline_tail + 4) & NEWLINE_BATCH_MASK;
        self.newline_len -= 4;
        self.line_start = nl3 + 1;
    }
};

const Chunk = struct {
    start: usize,
    end: usize,
    read_index_base: u64,
};

const WorkerResult = struct {
    stats: Stats = .{},
    first_error: ?ValidationError = null,
};

const SequenceKernelResult = struct {
    gc_bases: u64 = 0,
    invalid_index: ?usize = null,
};

const QualityKernelResult = struct {
    quality_sum: u64 = 0,
    invalid_index: ?usize = null,
};

const FusedRecordResult = struct {
    gc_bases: u64 = 0,
    quality_sum: u64 = 0,
};

const WorkerContext = struct {
    data: []const u8,
    queue: ?*ChunkQueue = null,
    profile: WorkProfile,
    mode: RunMode,
    result: *WorkerResult,
    progress: ?*ProgressState = null,
};

const ProducerContext = struct {
    data: []const u8,
    target_chunk_bytes: usize,
    queue: *ChunkQueue,
};

const ChunkQueue = struct {
    slots: []Chunk,
    head: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    tail: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    done: std.atomic.Value(bool) = std.atomic.Value(bool).init(false),

    fn init(allocator: std.mem.Allocator, capacity: usize) !ChunkQueue {
        return .{ .slots = try allocator.alloc(Chunk, capacity) };
    }

    fn deinit(self: *ChunkQueue, allocator: std.mem.Allocator) void {
        allocator.free(self.slots);
    }

    fn push(self: *ChunkQueue, chunk: Chunk) void {
        const capacity = self.slots.len;
        while (true) {
            const head = self.head.load(.acquire);
            const tail = self.tail.load(.acquire);
            if (head - tail < capacity) {
                self.slots[head % capacity] = chunk;
                self.head.store(head + 1, .release);
                return;
            }
            std.Thread.yield() catch {};
        }
    }

    fn markDone(self: *ChunkQueue) void {
        self.done.store(true, .release);
    }

    fn pop(self: *ChunkQueue) ?Chunk {
        while (true) {
            const tail = self.tail.load(.acquire);
            const head = self.head.load(.acquire);
            if (tail >= head) {
                const wait_prof_enabled = g_profiler.enabled;
                const wait_start = if (wait_prof_enabled) std.time.nanoTimestamp() else 0;
                if (self.done.load(.acquire)) {
                    if (wait_prof_enabled) profileEnd(.queue_wait, wait_start);
                    return null;
                }
                if (wait_prof_enabled) profileEnd(.queue_wait, wait_start);
                std.Thread.yield() catch {};
                continue;
            }
            if (self.tail.cmpxchgWeak(tail, tail + 1, .acq_rel, .acquire) == null) {
                return self.slots[tail % self.slots.len];
            }
        }
    }
};

const ProgressState = struct {
    enabled: bool,
    processed: usize = 0,
    total: usize = 0,
    last_percent: usize = 101,
    mutex: std.Thread.Mutex = .{},
};

const ProfileCounter = enum(usize) {
    validate_chunk,
    strict_framing,
    assume_framing,
    sequence_kernel,
    quality_kernel,
    newline_refill,
    queue_wait,
};

const InternalProfiler = struct {
    enabled: bool = false,
    nanos: [7]std.atomic.Value(u64) = .{
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
    },
    calls: [7]std.atomic.Value(u64) = .{
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
        std.atomic.Value(u64).init(0),
    },
};

var g_profiler = InternalProfiler{};

pub fn main() !void {
    const args = try std.process.argsAlloc(std.heap.page_allocator);
    defer std.process.argsFree(std.heap.page_allocator, args);

    if (args.len == 1) {
        try printUsage(true);
        std.process.exit(@intFromEnum(ExitCode.usage_error));
    }

    if (args.len == 2 and (std.mem.eql(u8, args[1], "--help") or std.mem.eql(u8, args[1], "-h"))) {
        try printUsage(false);
        return;
    }

    if (args.len == 2 and (std.mem.eql(u8, args[1], "--version") or std.mem.eql(u8, args[1], "-V"))) {
        try printVersion();
        return;
    }

    const run_options = parseRunOptions(args) catch {
        try printUsage(true);
        std.process.exit(@intFromEnum(ExitCode.usage_error));
    };
    if (run_options.operation == .explain) {
        runExplain(std.heap.page_allocator, run_options.input_path, run_options.json) catch |err| {
            try printStderr("error: explain failed: {s}\n", .{@errorName(err)});
            std.process.exit(@intFromEnum(ExitCode.io_error));
        };
        std.process.exit(@intFromEnum(ExitCode.ok));
    }
    if (run_options.operation == .compare) {
        runCompare(std.heap.page_allocator, run_options.input_path, run_options.compare_against_path.?, run_options.json) catch |err| {
            try printStderr("error: compare failed: {s}\n", .{@errorName(err)});
            std.process.exit(@intFromEnum(ExitCode.io_error));
        };
        std.process.exit(@intFromEnum(ExitCode.ok));
    }
    profilerReset(run_options.profile_internal);
    if (run_options.bench_kernels) {
        try runKernelBench();
        std.process.exit(@intFromEnum(ExitCode.ok));
    }
    const config = resolveProcessingConfig(run_options);
    const input_path = run_options.input_path;
    const loaded = loadInputFromPath(std.heap.page_allocator, input_path, run_options.io_preference, run_options.gzip_mode) catch |err| {
        try printStderr("error: failed to load '{s}': {s}\n", .{ input_path, @errorName(err) });
        std.process.exit(@intFromEnum(ExitCode.io_error));
    };
    defer loaded.deinit(std.heap.page_allocator);

    if (run_options.operation == .sample) {
        const summary = runSample(std.heap.page_allocator, loaded.bytes, run_options) catch |err| {
            try printStderr("error: sample failed: {s}\n", .{@errorName(err)});
            std.process.exit(@intFromEnum(ExitCode.io_error));
        };
        if (run_options.json) {
            const payload = .{
                .tool = "zdash",
                .version = VERSION,
                .schema_version = JSON_SCHEMA_VERSION,
                .operation = "sample",
                .status = "ok",
                .total_reads_seen = summary.total_reads_seen,
                .sampled_reads = summary.sampled_reads,
                .seed = summary.seed,
                .fraction = summary.fraction,
                .n = summary.n,
                .output_path = summary.output_path,
            };
            try printStdout("{f}\n", .{std.json.fmt(payload, .{})});
        }
        std.process.exit(@intFromEnum(ExitCode.ok));
    }

    if (run_options.operation == .repair) {
        const summary = runRepair(std.heap.page_allocator, loaded.bytes, run_options) catch |err| {
            try printStderr("error: repair failed: {s}\n", .{@errorName(err)});
            std.process.exit(@intFromEnum(ExitCode.io_error));
        };
        if (run_options.json) {
            const payload = .{
                .tool = "zdash",
                .version = VERSION,
                .schema_version = JSON_SCHEMA_VERSION,
                .operation = "repair",
                .status = "ok",
                .written_reads = summary.written_reads,
                .rejected_reads = summary.rejected_reads,
                .output_path = summary.output_path,
                .bad_records_path = summary.bad_records_path,
            };
            try printStdout("{f}\n", .{std.json.fmt(payload, .{})});
        } else if (!run_options.quiet) {
            try printStdout(
                "Repair complete: written={d}, rejected={d}, output={s}\n",
                .{ summary.written_reads, summary.rejected_reads, summary.output_path },
            );
        }
        std.process.exit(@intFromEnum(ExitCode.ok));
    }

    var bench: ?BenchResult = null;
    var outcome = ValidationOutcome{ .validation = .{ .ok = .{} } };
    defer outcome.deinit(std.heap.page_allocator);
    if (run_options.bench) {
        const validation = blk: {
            bench = runBench(std.heap.page_allocator, loaded.bytes, config) catch |err| {
                try printStderr("error: benchmark run failed: {s}\n", .{@errorName(err)});
                std.process.exit(@intFromEnum(ExitCode.io_error));
            };
            break :blk bench.?.validation;
        };
        outcome = .{ .validation = validation };
    } else {
        outcome = validateFastqDetailed(
            std.heap.page_allocator,
            loaded.bytes,
            config,
            run_options.show_progress,
            run_options.max_errors,
            run_options.extract_error_context_path,
            run_options.error_window_reads,
            run_options.extract_error_debug_path,
        ) catch |err| {
            try printStderr("error: processing failed: {s}\n", .{@errorName(err)});
            std.process.exit(@intFromEnum(ExitCode.io_error));
        };
    }
    if (run_options.json) {
        try printJsonSummary(input_path, loaded, outcome, config, run_options.operation);
    } else if (!run_options.quiet) {
        try printSummary(input_path, loaded, outcome.validation, config);
    } else if (outcome.validation == .fail) {
        for (outcome.errors) |err| try printValidationFailureDetailed(loaded.bytes, err);
        if (outcome.errors.len == 0) try printValidationFailureDetailed(loaded.bytes, outcome.validation.fail);
        if (outcome.context_output_path) |p| try printStderr("Error context written: {s}\n", .{p});
    }
    if (run_options.gha_annotations and outcome.validation == .fail) {
        for (outcome.errors) |err| try printGithubAnnotation(input_path, err);
        if (outcome.errors.len == 0) try printGithubAnnotation(input_path, outcome.validation.fail);
    }
    if (bench) |b| try printBenchSummary(loaded, config, b);
    if (run_options.profile_internal) try printInternalProfile();
    if (run_options.report_json_path) |report_path| {
        try writeValidationReportJson(report_path, input_path, loaded, outcome, config, run_options.operation);
    }

    switch (outcome.validation) {
        .ok => std.process.exit(@intFromEnum(ExitCode.ok)),
        .fail => std.process.exit(@intFromEnum(ExitCode.validation_error)),
    }
}

fn printUsage(to_stderr: bool) !void {
    if (to_stderr) {
        try printStderr(
            \\Usage:
            \\  zdash <check|scan|stats|repair|sample|explain|compare> [options] <input>
            \\  zdash [options] <input.fastq>
            \\
            \\Options:
            \\  -h, --help                Show help
            \\  -V, --version             Show version
            \\  --progress                Show byte progress while processing
            \\  --quiet                   Suppress main summary output
            \\  --gha-annotations         Emit GitHub Actions error annotations on validation failure
            \\  --json                    Emit machine-readable JSON output
            \\  --json-schema-version <v> Require an exact JSON schema version (current: 1.0.0)
            \\  --config <path>           Load defaults from config file (key=value lines)
            \\  --preset <strict-ci|fast-scan|qc-only>
            \\                            Apply common developer workflow defaults
            \\  --max-errors <N>          Stop after collecting N validation errors (default: 1)
            \\  --error-window-reads <N>  Reads before/after failure for context artifact (default: 5)
            \\  --extract-error-context[=<path>]
            \\                            Write failing read window to FASTQ file
            \\  --extract-error-debug <path>
            \\                            Write human debug report (context + byte hex window)
            \\  --report-json <path>      Write JSON report to file
            \\  --against <path>          Compare mode: second report/json path
            \\  --output <path|->         Output path for sample/repair (default: stdout)
            \\                            sample --json requires --output <path> (not '-')
            \\  --repair-mode <drop-bad-records|truncate-to-last-good>
            \\  --emit-bad-records <path> Write rejected reads during repair
            \\  --seed <N>                RNG seed for sampling (default: 42)
            \\  --fraction <F>            Sampling fraction (0.0 < F <= 1.0)
            \\  --n <N>                   Reservoir sample size
            \\  --threads <N>             Worker count (default: CPU count)
            \\  --chunk-bytes <N>         FASTQ chunk target size (default: 8388608)
            \\  --profile <full|validate-stats|stats-only>
            \\                            Work profile (default: full)
            \\  --mode <strict|assume-valid>
            \\                            Parser mode (default: strict)
            \\  --bench                   Run full pipeline benchmark mode
            \\  --io-mode <auto|mmap|buffered>
            \\                            Force I/O mode (default: auto)
            \\  --gzip-mode <stream|temp>
            \\                            gzip/bgzip path (default: stream; temp preserves legacy flow)
            \\  --bench-kernels           Run SIMD/scalar kernel microbenchmarks
            \\  --profile-internal        Print internal function timing counters
            \\
        , .{});
        return;
    }

    try printStdout(
        \\Usage:
        \\  zdash <check|scan|stats|repair|sample|explain|compare> [options] <input>
        \\  zdash [options] <input.fastq>
        \\
        \\Options:
        \\  -h, --help                Show help
        \\  -V, --version             Show version
        \\  --progress                Show byte progress while processing
        \\  --quiet                   Suppress main summary output
        \\  --gha-annotations         Emit GitHub Actions error annotations on validation failure
        \\  --json                    Emit machine-readable JSON output
        \\  --json-schema-version <v> Require an exact JSON schema version (current: 1.0.0)
        \\  --config <path>           Load defaults from config file (key=value lines)
        \\  --preset <strict-ci|fast-scan|qc-only>
        \\                            Apply common developer workflow defaults
        \\  --max-errors <N>          Stop after collecting N validation errors (default: 1)
        \\  --error-window-reads <N>  Reads before/after failure for context artifact (default: 5)
        \\  --extract-error-context[=<path>]
        \\                            Write failing read window to FASTQ file
        \\  --extract-error-debug <path>
        \\                            Write human debug report (context + byte hex window)
        \\  --report-json <path>      Write JSON report to file
        \\  --against <path>          Compare mode: second report/json path
        \\  --output <path|->         Output path for sample/repair (default: stdout)
        \\                            sample --json requires --output <path> (not '-')
        \\  --repair-mode <drop-bad-records|truncate-to-last-good>
        \\  --emit-bad-records <path> Write rejected reads during repair
        \\  --seed <N>                RNG seed for sampling (default: 42)
        \\  --fraction <F>            Sampling fraction (0.0 < F <= 1.0)
        \\  --n <N>                   Reservoir sample size
        \\  --threads <N>             Worker count (default: CPU count)
        \\  --chunk-bytes <N>         FASTQ chunk target size (default: 8388608)
        \\  --profile <full|validate-stats|stats-only>
        \\                            Work profile (default: full)
        \\  --mode <strict|assume-valid>
        \\                            Parser mode (default: strict)
        \\  --bench                   Run full pipeline benchmark mode
        \\  --io-mode <auto|mmap|buffered>
        \\                            Force I/O mode (default: auto)
        \\  --gzip-mode <stream|temp>
        \\                            gzip/bgzip path (default: stream; temp preserves legacy flow)
        \\  --bench-kernels           Run SIMD/scalar kernel microbenchmarks
        \\  --profile-internal        Print internal function timing counters
        \\
    , .{});
}

fn printVersion() !void {
    try printStdout("zdash {s}\n", .{VERSION});
}

fn printSummary(
    input_path: []const u8,
    loaded: LoadedInput,
    validation: ValidationResult,
    config: ProcessingConfig,
) !void {
    const rendered = try renderSummary(std.heap.page_allocator, input_path, loaded, validation, config);
    defer std.heap.page_allocator.free(rendered);
    try printStdout("{s}", .{rendered});
}

fn printJsonSummary(
    input_path: []const u8,
    loaded: LoadedInput,
    outcome: ValidationOutcome,
    config: ProcessingConfig,
    operation: Operation,
) !void {
    const json = try renderJsonSummary(std.heap.page_allocator, input_path, loaded, outcome, config, operation);
    defer std.heap.page_allocator.free(json);
    try printStdout("{s}\n", .{json});
}

fn renderJsonSummary(
    allocator: std.mem.Allocator,
    input_path: []const u8,
    loaded: LoadedInput,
    outcome: ValidationOutcome,
    config: ProcessingConfig,
    operation: Operation,
) ![]u8 {
    switch (outcome.validation) {
        .ok => |ok| {
            const payload = .{
                .tool = "zdash",
                .version = VERSION,
                .schema_version = JSON_SCHEMA_VERSION,
                .operation = operationLabel(operation),
                .status = "ok",
                .input = input_path,
                .io_mode = inputModeLabel(loaded.mode),
                .profile = workProfileLabel(config.profile),
                .mode = runModeLabel(config.mode),
                .threads = config.threads,
                .chunk_bytes = config.chunk_bytes,
                .file_bytes = loaded.file_size,
                .processed_bytes = loaded.bytes.len,
                .total_reads = ok.total_reads,
                .total_bases = ok.total_bases,
                .gc_bases = ok.gc_bases,
                .gc_percent = gcPercent(ok),
                .avg_quality = averageQuality(ok),
                .low_quality_reads = ok.low_quality_reads,
            };
            return std.fmt.allocPrint(allocator, "{f}", .{std.json.fmt(payload, .{})});
        },
        .fail => |fail| {
            const bad_char = validationBadChar(fail.issue);
            const line = validationIssueLine(fail.issue);
            const payload = .{
                .tool = "zdash",
                .version = VERSION,
                .schema_version = JSON_SCHEMA_VERSION,
                .operation = operationLabel(operation),
                .status = "failed",
                .input = input_path,
                .io_mode = inputModeLabel(loaded.mode),
                .profile = workProfileLabel(config.profile),
                .mode = runModeLabel(config.mode),
                .threads = config.threads,
                .chunk_bytes = config.chunk_bytes,
                .file_bytes = loaded.file_size,
                .max_errors = outcome.errors.len,
                .context_output_path = outcome.context_output_path,
                .@"error" = .{
                    .message = validationIssueMessage(fail.issue),
                    .line = line,
                    .read_index = fail.read_index,
                    .byte_offset = fail.byte_offset,
                    .bad_char = bad_char,
                },
            };
            return std.fmt.allocPrint(allocator, "{f}", .{std.json.fmt(payload, .{})});
        },
    }
}

fn renderSummary(
    allocator: std.mem.Allocator,
    input_path: []const u8,
    loaded: LoadedInput,
    validation: ValidationResult,
    config: ProcessingConfig,
) ![]u8 {
    switch (validation) {
        .ok => |ok| {
            const avg_q = averageQuality(ok);
            const gc_percent = gcPercent(ok);
            return std.fmt.allocPrint(allocator,
                \\Z-DASH: Zig DNA Shredder v{s}
                \\-----------------------------
                \\Input: {s}
                \\
                \\Summary:
                \\- Total Reads:  {d}
                \\- Avg Quality:  Q{d:.2}
                \\- GC Content:   {d:.2}%
                \\- Low Quality:  {d} reads (avg Q < {d:.2})
                \\- I/O Mode:     {s}
                \\- Profile:      {s}
                \\- Mode:         {s}
                \\- Threads:      {d}
                \\- Chunk Bytes:  {d}
                \\- File Bytes:   {d}
                \\- Proc Bytes:   {d}
                \\- Status:       PASSED (No corruption found)
                \\
            , .{
                VERSION,
                input_path,
                ok.total_reads,
                avg_q,
                gc_percent,
                ok.low_quality_reads,
                LOW_QUALITY_AVG_Q_THRESHOLD,
                inputModeLabel(loaded.mode),
                workProfileLabel(config.profile),
                runModeLabel(config.mode),
                config.threads,
                config.chunk_bytes,
                loaded.file_size,
                loaded.bytes.len,
            });
        },
        .fail => |fail| {
            return std.fmt.allocPrint(allocator,
                \\Z-DASH: Zig DNA Shredder v{s}
                \\-----------------------------
                \\Input: {s}
                \\
                \\Summary:
                \\- Total Reads:  N/A
                \\- Avg Quality:  N/A
                \\- GC Content:   N/A
                \\- I/O Mode:     {s}
                \\- Profile:      {s}
                \\- Mode:         {s}
                \\- Threads:      {d}
                \\- Chunk Bytes:  {d}
                \\- File Bytes:   {d}
                \\- Status:       FAILED (Corruption found)
                \\- First Error:  read #{d} at byte {d}: {s}
                \\
            , .{
                VERSION,
                input_path,
                inputModeLabel(loaded.mode),
                workProfileLabel(config.profile),
                runModeLabel(config.mode),
                config.threads,
                config.chunk_bytes,
                loaded.file_size,
                fail.read_index,
                fail.byte_offset,
                validationIssueMessage(fail.issue),
            });
        },
    }
}

fn printError(msg: []const u8) !void {
    try printStderr("{s}", .{msg});
}

fn parsePreset(value: []const u8) !Preset {
    if (std.mem.eql(u8, value, "strict-ci")) return .strict_ci;
    if (std.mem.eql(u8, value, "fast-scan")) return .fast_scan;
    if (std.mem.eql(u8, value, "qc-only")) return .qc_only;
    return error.InvalidPreset;
}

fn parseBool(value: []const u8) !bool {
    if (std.mem.eql(u8, value, "true")) return true;
    if (std.mem.eql(u8, value, "false")) return false;
    return error.InvalidBoolean;
}

fn trimConfigValue(value: []const u8) []const u8 {
    const trimmed = std.mem.trim(u8, value, " \t\r");
    if (trimmed.len >= 2) {
        if (trimmed[0] == '"' and trimmed[trimmed.len - 1] == '"') return trimmed[1 .. trimmed.len - 1];
        if (trimmed[0] == '\'' and trimmed[trimmed.len - 1] == '\'') return trimmed[1 .. trimmed.len - 1];
    }
    return trimmed;
}

fn applyPreset(run: *RunOptions, preset: Preset) void {
    run.preset = preset;
    run.operation_explicit = true;
    switch (preset) {
        .strict_ci => {
            run.operation = .check;
            run.mode = .strict;
            run.mode_explicit = true;
            run.profile = .full;
            run.profile_explicit = true;
            run.json = true;
        },
        .fast_scan => {
            run.operation = .scan;
            run.mode = .assume_valid;
            run.mode_explicit = true;
            run.profile = .validate_stats;
            run.profile_explicit = true;
        },
        .qc_only => {
            run.operation = .stats;
            run.mode = .assume_valid;
            run.mode_explicit = true;
            run.profile = .stats_only;
            run.profile_explicit = true;
        },
    }
}

fn applyConfigOption(run: *RunOptions, key: []const u8, value: []const u8) !void {
    if (std.mem.eql(u8, key, "preset")) {
        applyPreset(run, try parsePreset(value));
        return;
    }
    if (std.mem.eql(u8, key, "operation")) {
        if (std.mem.eql(u8, value, "check")) {
            run.operation = .check;
        } else if (std.mem.eql(u8, value, "scan")) {
            run.operation = .scan;
        } else if (std.mem.eql(u8, value, "stats")) {
            run.operation = .stats;
        } else if (std.mem.eql(u8, value, "repair")) {
            run.operation = .repair;
        } else if (std.mem.eql(u8, value, "sample")) {
            run.operation = .sample;
        } else if (std.mem.eql(u8, value, "explain")) {
            run.operation = .explain;
        } else if (std.mem.eql(u8, value, "compare")) {
            run.operation = .compare;
        } else {
            return error.InvalidOperation;
        }
        run.operation_explicit = true;
        return;
    }
    if (std.mem.eql(u8, key, "profile")) {
        run.profile = try parseWorkProfile(value);
        run.profile_explicit = true;
        return;
    }
    if (std.mem.eql(u8, key, "mode")) {
        run.mode = try parseRunMode(value);
        run.mode_explicit = true;
        return;
    }
    if (std.mem.eql(u8, key, "io_mode")) {
        run.io_preference = try parseIoPreference(value);
        return;
    }
    if (std.mem.eql(u8, key, "gzip_mode")) {
        run.gzip_mode = try parseGzipMode(value);
        return;
    }
    if (std.mem.eql(u8, key, "threads")) {
        run.threads = try parseThreadCount(value);
        return;
    }
    if (std.mem.eql(u8, key, "chunk_bytes")) {
        run.chunk_bytes = try parseChunkBytes(value);
        return;
    }
    if (std.mem.eql(u8, key, "max_errors")) {
        run.max_errors = try parseMaxErrors(value);
        return;
    }
    if (std.mem.eql(u8, key, "error_window_reads")) {
        run.error_window_reads = try std.fmt.parseInt(u64, value, 10);
        if (run.error_window_reads == 0) return error.InvalidWindowReads;
        return;
    }
    if (std.mem.eql(u8, key, "seed")) {
        run.sample_seed = try std.fmt.parseInt(u64, value, 10);
        return;
    }
    if (std.mem.eql(u8, key, "fraction")) {
        run.sample_fraction = try parseFraction(value);
        return;
    }
    if (std.mem.eql(u8, key, "n")) {
        run.sample_n = try parseMaxErrors(value);
        return;
    }
    if (std.mem.eql(u8, key, "json")) {
        run.json = try parseBool(value);
        return;
    }
    if (std.mem.eql(u8, key, "quiet")) {
        run.quiet = try parseBool(value);
        return;
    }
    if (std.mem.eql(u8, key, "progress")) {
        run.show_progress = try parseBool(value);
        return;
    }
    if (std.mem.eql(u8, key, "gha_annotations")) {
        run.gha_annotations = try parseBool(value);
        return;
    }
    if (std.mem.eql(u8, key, "json_schema_version")) {
        run.json_schema_version = value;
        return;
    }
    if (std.mem.eql(u8, key, "report_json")) {
        run.report_json_path = value;
        return;
    }
    if (std.mem.eql(u8, key, "output")) {
        run.output_path = value;
        return;
    }
    if (std.mem.eql(u8, key, "extract_error_context")) {
        run.extract_error_context_path = value;
        return;
    }
    if (std.mem.eql(u8, key, "extract_error_debug")) {
        run.extract_error_debug_path = value;
        return;
    }
    if (std.mem.eql(u8, key, "against")) {
        run.compare_against_path = value;
        return;
    }
    if (std.mem.eql(u8, key, "repair_mode")) {
        run.repair_mode = try parseRepairMode(value);
        return;
    }
    if (std.mem.eql(u8, key, "emit_bad_records")) {
        run.emit_bad_records_path = value;
        return;
    }
    return error.UnknownConfigKey;
}

fn loadRunOptionsConfig(allocator: std.mem.Allocator, run: *RunOptions, path: []const u8) !void {
    const bytes = try std.fs.cwd().readFileAlloc(allocator, path, 1024 * 1024);

    var it = std.mem.tokenizeScalar(u8, bytes, '\n');
    while (it.next()) |raw_line| {
        const line = std.mem.trim(u8, raw_line, " \t\r");
        if (line.len == 0 or line[0] == '#') continue;
        const eq = std.mem.indexOfScalar(u8, line, '=') orelse return error.InvalidConfigLine;
        const key = std.mem.trim(u8, line[0..eq], " \t\r");
        const value = trimConfigValue(line[eq + 1 ..]);
        if (key.len == 0) return error.InvalidConfigLine;
        try applyConfigOption(run, key, value);
    }
}

fn parseRunOptions(args: []const []const u8) !RunOptions {
    var run = RunOptions{ .input_path = "" };
    var config_path: ?[]const u8 = null;

    var pre_i: usize = 1;
    while (pre_i < args.len) : (pre_i += 1) {
        const arg = args[pre_i];
        if (std.mem.eql(u8, arg, "--config")) {
            if (pre_i + 1 >= args.len) {
                try printError("error: --config requires a path\n");
                return error.InvalidUsage;
            }
            pre_i += 1;
            config_path = args[pre_i];
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--config=")) {
            config_path = arg["--config=".len..];
            continue;
        }
    }
    if (config_path) |path| {
        run.config_path = path;
        loadRunOptionsConfig(std.heap.page_allocator, &run, path) catch |err| {
            try printStderr("error: failed to load config '{s}': {s}\n", .{ path, @errorName(err) });
            return error.InvalidUsage;
        };
    }

    var i: usize = 1;
    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "--progress")) {
            run.show_progress = true;
            continue;
        }

        if (std.mem.eql(u8, arg, "--quiet")) {
            run.quiet = true;
            continue;
        }
        if (std.mem.eql(u8, arg, "--gha-annotations")) {
            run.gha_annotations = true;
            continue;
        }
        if (std.mem.eql(u8, arg, "--json")) {
            run.json = true;
            continue;
        }
        if (std.mem.eql(u8, arg, "--json-schema-version")) {
            if (i + 1 >= args.len) {
                try printError("error: --json-schema-version requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.json_schema_version = args[i];
            continue;
        }
        if (std.mem.eql(u8, arg, "--preset")) {
            if (i + 1 >= args.len) {
                try printError("error: --preset requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            const preset = parsePreset(args[i]) catch {
                try printStderr("error: invalid preset '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            applyPreset(&run, preset);
            continue;
        }
        if (std.mem.eql(u8, arg, "--config")) {
            if (i + 1 >= args.len) {
                try printError("error: --config requires a path\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.config_path = args[i];
            continue;
        }
        if (std.mem.eql(u8, arg, "--max-errors")) {
            if (i + 1 >= args.len) {
                try printError("error: --max-errors requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.max_errors = parseMaxErrors(args[i]) catch {
                try printStderr("error: invalid max-errors '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            continue;
        }
        if (std.mem.eql(u8, arg, "--error-window-reads")) {
            if (i + 1 >= args.len) {
                try printError("error: --error-window-reads requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.error_window_reads = std.fmt.parseInt(u64, args[i], 10) catch {
                try printStderr("error: invalid error-window-reads '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            if (run.error_window_reads == 0) {
                try printError("error: --error-window-reads must be > 0\n");
                return error.InvalidUsage;
            }
            continue;
        }
        if (std.mem.eql(u8, arg, "--extract-error-context")) {
            if (i + 1 < args.len and !std.mem.startsWith(u8, args[i + 1], "-")) {
                i += 1;
                run.extract_error_context_path = args[i];
            } else {
                run.extract_error_context_path = "zdash_error_context.fastq";
            }
            continue;
        }
        if (std.mem.eql(u8, arg, "--extract-error-debug")) {
            if (i + 1 >= args.len) {
                try printError("error: --extract-error-debug requires a path\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.extract_error_debug_path = args[i];
            continue;
        }
        if (std.mem.eql(u8, arg, "--report-json")) {
            if (i + 1 >= args.len) {
                try printError("error: --report-json requires a path\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.report_json_path = args[i];
            continue;
        }
        if (std.mem.eql(u8, arg, "--against")) {
            if (i + 1 >= args.len) {
                try printError("error: --against requires a path\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.compare_against_path = args[i];
            continue;
        }
        if (std.mem.eql(u8, arg, "--output")) {
            if (i + 1 >= args.len) {
                try printError("error: --output requires a path\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.output_path = args[i];
            continue;
        }
        if (std.mem.eql(u8, arg, "--repair-mode")) {
            if (i + 1 >= args.len) {
                try printError("error: --repair-mode requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.repair_mode = parseRepairMode(args[i]) catch {
                try printStderr("error: invalid repair-mode '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            continue;
        }
        if (std.mem.eql(u8, arg, "--emit-bad-records")) {
            if (i + 1 >= args.len) {
                try printError("error: --emit-bad-records requires a path\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.emit_bad_records_path = args[i];
            continue;
        }
        if (std.mem.eql(u8, arg, "--seed")) {
            if (i + 1 >= args.len) {
                try printError("error: --seed requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.sample_seed = std.fmt.parseInt(u64, args[i], 10) catch {
                try printStderr("error: invalid seed '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            continue;
        }
        if (std.mem.eql(u8, arg, "--fraction")) {
            if (i + 1 >= args.len) {
                try printError("error: --fraction requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.sample_fraction = parseFraction(args[i]) catch {
                try printStderr("error: invalid fraction '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            continue;
        }
        if (std.mem.eql(u8, arg, "--n")) {
            if (i + 1 >= args.len) {
                try printError("error: --n requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.sample_n = parseMaxErrors(args[i]) catch {
                try printStderr("error: invalid n '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            continue;
        }

        if (std.mem.eql(u8, arg, "--bench")) {
            run.bench = true;
            continue;
        }

        if (std.mem.eql(u8, arg, "--bench-kernels")) {
            run.bench_kernels = true;
            continue;
        }

        if (std.mem.eql(u8, arg, "--profile-internal")) {
            run.profile_internal = true;
            continue;
        }

        if (std.mem.eql(u8, arg, "--io-mode")) {
            if (i + 1 >= args.len) {
                try printError("error: --io-mode requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.io_preference = parseIoPreference(args[i]) catch {
                try printStderr("error: invalid io mode '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            continue;
        }
        if (std.mem.eql(u8, arg, "--gzip-mode")) {
            if (i + 1 >= args.len) {
                try printError("error: --gzip-mode requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.gzip_mode = parseGzipMode(args[i]) catch {
                try printStderr("error: invalid gzip mode '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            continue;
        }

        if (std.mem.eql(u8, arg, "--threads")) {
            if (i + 1 >= args.len) {
                try printError("error: --threads requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.threads = parseThreadCount(args[i]) catch {
                try printStderr("error: invalid thread count '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            continue;
        }

        if (std.mem.eql(u8, arg, "--chunk-bytes")) {
            if (i + 1 >= args.len) {
                try printError("error: --chunk-bytes requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.chunk_bytes = parseChunkBytes(args[i]) catch {
                try printStderr("error: invalid chunk size '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            continue;
        }

        if (std.mem.eql(u8, arg, "--profile")) {
            if (i + 1 >= args.len) {
                try printError("error: --profile requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.profile = parseWorkProfile(args[i]) catch {
                try printStderr("error: invalid profile '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            run.profile_explicit = true;
            continue;
        }

        if (std.mem.eql(u8, arg, "--mode")) {
            if (i + 1 >= args.len) {
                try printError("error: --mode requires a value\n");
                return error.InvalidUsage;
            }
            i += 1;
            run.mode = parseRunMode(args[i]) catch {
                try printStderr("error: invalid mode '{s}'\n", .{args[i]});
                return error.InvalidUsage;
            };
            run.mode_explicit = true;
            continue;
        }

        if (std.mem.startsWith(u8, arg, "--io-mode=")) {
            const value = arg["--io-mode=".len..];
            run.io_preference = parseIoPreference(value) catch {
                try printStderr("error: invalid io mode '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--gzip-mode=")) {
            const value = arg["--gzip-mode=".len..];
            run.gzip_mode = parseGzipMode(value) catch {
                try printStderr("error: invalid gzip mode '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            continue;
        }

        if (std.mem.startsWith(u8, arg, "--chunk-bytes=")) {
            const value = arg["--chunk-bytes=".len..];
            run.chunk_bytes = parseChunkBytes(value) catch {
                try printStderr("error: invalid chunk size '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            continue;
        }

        if (std.mem.startsWith(u8, arg, "--threads=")) {
            const value = arg["--threads=".len..];
            run.threads = parseThreadCount(value) catch {
                try printStderr("error: invalid thread count '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            continue;
        }

        if (std.mem.startsWith(u8, arg, "--profile=")) {
            const value = arg["--profile=".len..];
            run.profile = parseWorkProfile(value) catch {
                try printStderr("error: invalid profile '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            run.profile_explicit = true;
            continue;
        }

        if (std.mem.startsWith(u8, arg, "--mode=")) {
            const value = arg["--mode=".len..];
            run.mode = parseRunMode(value) catch {
                try printStderr("error: invalid mode '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            run.mode_explicit = true;
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--preset=")) {
            const value = arg["--preset=".len..];
            const preset = parsePreset(value) catch {
                try printStderr("error: invalid preset '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            applyPreset(&run, preset);
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--json-schema-version=")) {
            run.json_schema_version = arg["--json-schema-version=".len..];
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--config=")) {
            run.config_path = arg["--config=".len..];
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--max-errors=")) {
            const value = arg["--max-errors=".len..];
            run.max_errors = parseMaxErrors(value) catch {
                try printStderr("error: invalid max-errors '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--error-window-reads=")) {
            const value = arg["--error-window-reads=".len..];
            run.error_window_reads = std.fmt.parseInt(u64, value, 10) catch {
                try printStderr("error: invalid error-window-reads '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            if (run.error_window_reads == 0) {
                try printError("error: --error-window-reads must be > 0\n");
                return error.InvalidUsage;
            }
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--extract-error-context=")) {
            const value = arg["--extract-error-context=".len..];
            run.extract_error_context_path = value;
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--extract-error-debug=")) {
            run.extract_error_debug_path = arg["--extract-error-debug=".len..];
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--report-json=")) {
            run.report_json_path = arg["--report-json=".len..];
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--against=")) {
            run.compare_against_path = arg["--against=".len..];
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--output=")) {
            run.output_path = arg["--output=".len..];
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--repair-mode=")) {
            const value = arg["--repair-mode=".len..];
            run.repair_mode = parseRepairMode(value) catch {
                try printStderr("error: invalid repair-mode '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--emit-bad-records=")) {
            run.emit_bad_records_path = arg["--emit-bad-records=".len..];
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--seed=")) {
            const value = arg["--seed=".len..];
            run.sample_seed = std.fmt.parseInt(u64, value, 10) catch {
                try printStderr("error: invalid seed '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--fraction=")) {
            const value = arg["--fraction=".len..];
            run.sample_fraction = parseFraction(value) catch {
                try printStderr("error: invalid fraction '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            continue;
        }
        if (std.mem.startsWith(u8, arg, "--n=")) {
            const value = arg["--n=".len..];
            run.sample_n = parseMaxErrors(value) catch {
                try printStderr("error: invalid n '{s}'\n", .{value});
                return error.InvalidUsage;
            };
            continue;
        }

        if (std.mem.startsWith(u8, arg, "-") and !std.mem.eql(u8, arg, "-")) {
            try printStderr("error: unknown option '{s}'\n", .{arg});
            return error.InvalidUsage;
        }

        if (run.input_path.len == 0 and std.mem.eql(u8, arg, "check")) {
            run.operation = .check;
            run.operation_explicit = true;
            continue;
        }
        if (run.input_path.len == 0 and std.mem.eql(u8, arg, "scan")) {
            run.operation = .scan;
            run.operation_explicit = true;
            continue;
        }
        if (run.input_path.len == 0 and std.mem.eql(u8, arg, "stats")) {
            run.operation = .stats;
            run.operation_explicit = true;
            continue;
        }
        if (run.input_path.len == 0 and std.mem.eql(u8, arg, "repair")) {
            run.operation = .repair;
            run.operation_explicit = true;
            continue;
        }
        if (run.input_path.len == 0 and std.mem.eql(u8, arg, "sample")) {
            run.operation = .sample;
            run.operation_explicit = true;
            continue;
        }
        if (run.input_path.len == 0 and std.mem.eql(u8, arg, "explain")) {
            run.operation = .explain;
            run.operation_explicit = true;
            continue;
        }
        if (run.input_path.len == 0 and std.mem.eql(u8, arg, "compare")) {
            run.operation = .compare;
            run.operation_explicit = true;
            continue;
        }

        if (run.input_path.len != 0) {
            try printError("error: expected exactly one FASTQ input path\n");
            return error.InvalidUsage;
        }
        run.input_path = arg;
    }

    if (run.input_path.len == 0) {
        if (run.bench_kernels) return run;
        try printError("error: expected FASTQ input path\n");
        return error.InvalidUsage;
    }
    if (run.operation == .sample) {
        const has_fraction = run.sample_fraction != null;
        const has_n = run.sample_n != null;
        if (has_fraction == has_n) {
            try printError("error: sample requires exactly one of --fraction or --n\n");
            return error.InvalidUsage;
        }
        if (run.json) {
            if (run.output_path == null) {
                try printError("error: sample --json requires --output <path>\n");
                return error.InvalidUsage;
            }
            if (std.mem.eql(u8, run.output_path.?, "-")) {
                try printError("error: sample --json requires --output to be a file path (not '-')\n");
                return error.InvalidUsage;
            }
        }
    }
    if (run.operation == .compare and run.compare_against_path == null) {
        try printError("error: compare requires --against <path>\n");
        return error.InvalidUsage;
    }
    if (run.operation == .explain and run.compare_against_path != null) {
        try printError("error: --against is only valid with compare\n");
        return error.InvalidUsage;
    }
    if (run.json_schema_version) |version| {
        if (!std.mem.eql(u8, version, JSON_SCHEMA_VERSION)) {
            try printStderr(
                "error: unsupported json schema version '{s}' (supported: {s})\n",
                .{ version, JSON_SCHEMA_VERSION },
            );
            return error.InvalidUsage;
        }
    }

    return run;
}

fn parseIoPreference(value: []const u8) !IoPreference {
    if (std.mem.eql(u8, value, "auto")) return .auto;
    if (std.mem.eql(u8, value, "mmap")) return .mmap;
    if (std.mem.eql(u8, value, "buffered")) return .buffered;
    return error.InvalidIoPreference;
}

fn parseGzipMode(value: []const u8) !GzipMode {
    if (std.mem.eql(u8, value, "stream")) return .stream;
    if (std.mem.eql(u8, value, "temp")) return .temp;
    return error.InvalidGzipMode;
}

fn parseThreadCount(value: []const u8) !usize {
    const parsed = std.fmt.parseInt(usize, value, 10) catch return error.InvalidThreadCount;
    if (parsed == 0 or parsed > MAX_THREADS) return error.InvalidThreadCount;
    return parsed;
}

fn parseChunkBytes(value: []const u8) !usize {
    const parsed = std.fmt.parseInt(usize, value, 10) catch return error.InvalidChunkSize;
    if (parsed == 0) return error.InvalidChunkSize;
    return parsed;
}

fn parseMaxErrors(value: []const u8) !usize {
    const parsed = std.fmt.parseInt(usize, value, 10) catch return error.InvalidMaxErrors;
    if (parsed == 0) return error.InvalidMaxErrors;
    return parsed;
}

fn parseFraction(value: []const u8) !f64 {
    const parsed = std.fmt.parseFloat(f64, value) catch return error.InvalidFraction;
    if (!(parsed > 0.0 and parsed <= 1.0)) return error.InvalidFraction;
    return parsed;
}

fn parseRepairMode(value: []const u8) !RepairMode {
    if (std.mem.eql(u8, value, "drop-bad-records")) return .drop_bad_records;
    if (std.mem.eql(u8, value, "truncate-to-last-good")) return .truncate_to_last_good;
    return error.InvalidRepairMode;
}

fn parseWorkProfile(value: []const u8) !WorkProfile {
    if (std.mem.eql(u8, value, "full")) return .full;
    if (std.mem.eql(u8, value, "validate-stats")) return .validate_stats;
    if (std.mem.eql(u8, value, "stats-only")) return .stats_only;
    return error.InvalidProfile;
}

fn parseRunMode(value: []const u8) !RunMode {
    if (std.mem.eql(u8, value, "strict")) return .strict;
    if (std.mem.eql(u8, value, "assume-valid")) return .assume_valid;
    return error.InvalidMode;
}

fn resolveProcessingConfig(run: RunOptions) ProcessingConfig {
    const cpu_count = std.Thread.getCpuCount() catch 1;
    const auto_threads = if (cpu_count == 0) @as(usize, 1) else @min(cpu_count, MAX_THREADS);
    const default_profile: WorkProfile = switch (run.operation) {
        .check => .full,
        .scan => .full,
        .stats => .stats_only,
        .repair => .full,
        .sample => .stats_only,
        .explain => .full,
        .compare => .full,
    };
    const default_mode: RunMode = switch (run.operation) {
        .check => .strict,
        .scan => .assume_valid,
        .stats => .assume_valid,
        .repair => .strict,
        .sample => .assume_valid,
        .explain => .strict,
        .compare => .strict,
    };
    return .{
        .chunk_bytes = run.chunk_bytes,
        .threads = if (run.threads == 0) auto_threads else run.threads,
        .profile = if (run.profile_explicit) run.profile else default_profile,
        .mode = if (run.mode_explicit) run.mode else default_mode,
    };
}

fn runKernelBench() !void {
    const allocator = std.heap.page_allocator;
    const bytes_per_buffer: usize = 8 * 1024 * 1024;
    const iterations: usize = 64;

    const sequence = try allocator.alloc(u8, bytes_per_buffer);
    defer allocator.free(sequence);
    const quality = try allocator.alloc(u8, bytes_per_buffer);
    defer allocator.free(quality);

    for (sequence, 0..) |*ch, i| {
        ch.* = switch (i % 5) {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            else => 'N',
        };
    }
    for (quality, 0..) |*ch, i| {
        ch.* = @as(u8, 33) + @as(u8, @intCast(i % 41));
    }

    var checksum: u64 = 0;
    var timer = try std.time.Timer.start();
    for (0..iterations) |_| {
        const s = analyzeSequenceScalar(sequence);
        const q = analyzeQualityScalar(quality);
        checksum +%= s.gc_bases + q.quality_sum;
    }
    const scalar_ns = timer.read();

    var simd_ns: u64 = 0;
    if (SIMD_AVAILABLE) {
        timer = try std.time.Timer.start();
        for (0..iterations) |_| {
            const s = analyzeSequenceSimd(sequence);
            const q = analyzeQualitySimd(quality);
            checksum +%= s.gc_bases + q.quality_sum;
        }
        simd_ns = timer.read();
    }

    const processed_bytes = @as(f64, @floatFromInt(bytes_per_buffer * iterations * 2));
    const gib = 1024.0 * 1024.0 * 1024.0;
    const scalar_gib_s = (processed_bytes / gib) / (@as(f64, @floatFromInt(scalar_ns)) / 1_000_000_000.0);

    try printStdout("Kernel Microbench\n", .{});
    try printStdout("-----------------\n", .{});
    try printStdout("SIMD width: {d} bytes\n", .{SIMD_WIDTH});
    try printStdout("Buffers: {d} bytes each, iterations: {d}\n", .{ bytes_per_buffer, iterations });
    try printStdout("Scalar: {d:.2} GiB/s ({d} ns)\n", .{ scalar_gib_s, scalar_ns });
    if (SIMD_AVAILABLE) {
        const simd_gib_s = (processed_bytes / gib) / (@as(f64, @floatFromInt(simd_ns)) / 1_000_000_000.0);
        try printStdout("SIMD:   {d:.2} GiB/s ({d} ns)\n", .{ simd_gib_s, simd_ns });
    } else {
        try printStdout("SIMD:   unavailable on this target\n", .{});
    }
    try printStdout("Checksum: {d}\n", .{checksum});
}


fn runBench(allocator: std.mem.Allocator, data: []const u8, config: ProcessingConfig) !BenchResult {
    const warmup_runs: usize = 1;
    const measured_runs: usize = 3;
    var baseline: ?ValidationResult = null;

    for (0..warmup_runs) |_| {
        const result = try validateFastq(allocator, data, config, false);
        if (baseline == null) baseline = result;
    }

    var best_ns: u64 = std.math.maxInt(u64);
    var total_ns: u128 = 0;
    var last_result: ValidationResult = baseline.?;
    for (0..measured_runs) |_| {
        var timer = try std.time.Timer.start();
        const result = try validateFastq(allocator, data, config, false);
        const elapsed = timer.read();
        if (!validationResultEqual(baseline.?, result)) return error.NonDeterministicResult;

        total_ns += elapsed;
        if (elapsed < best_ns) {
            best_ns = elapsed;
            last_result = result;
        }
    }

    return .{
        .warmup_runs = warmup_runs,
        .measured_runs = measured_runs,
        .best_ns = best_ns,
        .avg_ns = @as(u64, @intCast(total_ns / measured_runs)),
        .validation = last_result,
    };
}

fn validationResultEqual(a: ValidationResult, b: ValidationResult) bool {
    return switch (a) {
        .ok => |a_stats| switch (b) {
            .ok => |b_stats| std.meta.eql(a_stats, b_stats),
            .fail => false,
        },
        .fail => |a_fail| switch (b) {
            .ok => false,
            .fail => |b_fail| std.meta.eql(a_fail, b_fail),
        },
    };
}

fn printBenchSummary(loaded: LoadedInput, config: ProcessingConfig, bench: BenchResult) !void {
    const gib = 1024.0 * 1024.0 * 1024.0;
    const avg_secs = @as(f64, @floatFromInt(bench.avg_ns)) / 1_000_000_000.0;
    const best_secs = @as(f64, @floatFromInt(bench.best_ns)) / 1_000_000_000.0;
    const gib_total = @as(f64, @floatFromInt(loaded.bytes.len)) / gib;
    const avg_gib_s = if (avg_secs == 0) 0.0 else gib_total / avg_secs;
    const best_gib_s = if (best_secs == 0) 0.0 else gib_total / best_secs;
    const core_count = std.Thread.getCpuCount() catch 1;

    try printStdout("\nBenchmark:\n", .{});
    try printStdout("- Warmup Runs:  {d}\n", .{bench.warmup_runs});
    try printStdout("- Measured:     {d}\n", .{bench.measured_runs});
    try printStdout("- Avg Time:     {d:.3}s\n", .{avg_secs});
    try printStdout("- Best Time:    {d:.3}s\n", .{best_secs});
    try printStdout("- Avg Throughput:  {d:.2} GiB/s\n", .{avg_gib_s});
    try printStdout("- Best Throughput: {d:.2} GiB/s\n", .{best_gib_s});
    switch (bench.validation) {
        .ok => |stats| {
            const reads_per_s = if (avg_secs == 0) 0.0 else @as(f64, @floatFromInt(stats.total_reads)) / avg_secs;
            try printStdout("- Reads/s:      {d:.2}\n", .{reads_per_s});
        },
        .fail => try printStdout("- Reads/s:      N/A (validation failed)\n", .{}),
    }
    try printStdout("- CPU Cores:    {d}\n", .{core_count});
    try printStdout("- CPU Arch:     {s}\n", .{@tagName(builtin.cpu.arch)});
    try printStdout("- OS:           {s}\n", .{@tagName(builtin.os.tag)});
    try printStdout("- Storage:      {s}\n", .{storageLabel(loaded.file_kind)});
    try printStdout("- Zig:          {s}\n", .{builtin.zig_version_string});
    try printStdout("- Config:       threads={d}, chunk_bytes={d}, io_mode={s}\n", .{
        config.threads,
        config.chunk_bytes,
        inputModeLabel(loaded.mode),
    });
    try printStdout("- Profile:      {s}\n", .{workProfileLabel(config.profile)});
    try printStdout("- Mode:         {s}\n", .{runModeLabel(config.mode)});
}

fn profilerReset(enabled: bool) void {
    g_profiler.enabled = enabled;
    for (&g_profiler.nanos) |*n| n.store(0, .monotonic);
    for (&g_profiler.calls) |*c| c.store(0, .monotonic);
}

fn profileStart() i128 {
    if (!g_profiler.enabled) return 0;
    return std.time.nanoTimestamp();
}

fn profileEnd(counter: ProfileCounter, start_ns: i128) void {
    if (!g_profiler.enabled) return;
    const elapsed: u64 = @intCast(@max(@as(i128, 0), std.time.nanoTimestamp() - start_ns));
    _ = g_profiler.nanos[@intFromEnum(counter)].fetchAdd(elapsed, .monotonic);
    _ = g_profiler.calls[@intFromEnum(counter)].fetchAdd(1, .monotonic);
}

fn printInternalProfile() !void {
    var items = [_]struct { name: []const u8, ns: u64, calls: u64 }{
        .{ .name = "validateFastqChunk", .ns = g_profiler.nanos[@intFromEnum(ProfileCounter.validate_chunk)].load(.monotonic), .calls = g_profiler.calls[@intFromEnum(ProfileCounter.validate_chunk)].load(.monotonic) },
        .{ .name = "strictFraming", .ns = g_profiler.nanos[@intFromEnum(ProfileCounter.strict_framing)].load(.monotonic), .calls = g_profiler.calls[@intFromEnum(ProfileCounter.strict_framing)].load(.monotonic) },
        .{ .name = "assumeFraming", .ns = g_profiler.nanos[@intFromEnum(ProfileCounter.assume_framing)].load(.monotonic), .calls = g_profiler.calls[@intFromEnum(ProfileCounter.assume_framing)].load(.monotonic) },
        .{ .name = "analyzeSequenceKernel", .ns = g_profiler.nanos[@intFromEnum(ProfileCounter.sequence_kernel)].load(.monotonic), .calls = g_profiler.calls[@intFromEnum(ProfileCounter.sequence_kernel)].load(.monotonic) },
        .{ .name = "analyzeQualityKernel", .ns = g_profiler.nanos[@intFromEnum(ProfileCounter.quality_kernel)].load(.monotonic), .calls = g_profiler.calls[@intFromEnum(ProfileCounter.quality_kernel)].load(.monotonic) },
        .{ .name = "newlineRefill", .ns = g_profiler.nanos[@intFromEnum(ProfileCounter.newline_refill)].load(.monotonic), .calls = g_profiler.calls[@intFromEnum(ProfileCounter.newline_refill)].load(.monotonic) },
        .{ .name = "queueWait", .ns = g_profiler.nanos[@intFromEnum(ProfileCounter.queue_wait)].load(.monotonic), .calls = g_profiler.calls[@intFromEnum(ProfileCounter.queue_wait)].load(.monotonic) },
    };
    std.sort.insertion(
        @TypeOf(items[0]),
        &items,
        {},
        struct {
            fn lessThan(_: void, a: @TypeOf(items[0]), b: @TypeOf(items[0])) bool {
                return a.ns > b.ns;
            }
        }.lessThan,
    );

    try printStdout("\nInternal Profile (approx wall-time):\n", .{});
    var total: u64 = 0;
    for (items) |it| total += it.ns;
    for (items) |it| {
        if (it.ns == 0) continue;
        const pct = if (total == 0) 0.0 else (@as(f64, @floatFromInt(it.ns)) * 100.0) / @as(f64, @floatFromInt(total));
        try printStdout("- {s}: {d:.3}s ({d:.1}%), calls={d}\n", .{
            it.name,
            @as(f64, @floatFromInt(it.ns)) / 1_000_000_000.0,
            pct,
            it.calls,
        });
    }
}

fn storageLabel(kind: std.fs.File.Kind) []const u8 {
    return switch (kind) {
        .file => "regular-file",
        .directory => "directory",
        .character_device => "character-device",
        .block_device => "block-device",
        .named_pipe => "pipe",
        .sym_link => "symlink",
        .unix_domain_socket => "unix-socket",
        .whiteout => "whiteout",
        .door => "door",
        .event_port => "event-port",
        .unknown => "unknown",
    };
}

fn printStdout(comptime fmt: []const u8, args: anytype) !void {
    var buffer: [4096]u8 = undefined;
    var writer = std.fs.File.stdout().writer(&buffer);
    const out = &writer.interface;
    try out.print(fmt, args);
    try out.flush();
}

fn printStderr(comptime fmt: []const u8, args: anytype) !void {
    var buffer: [2048]u8 = undefined;
    var writer = std.fs.File.stderr().writer(&buffer);
    const out = &writer.interface;
    try out.print(fmt, args);
    try out.flush();
}

fn loadInputFromPath(
    allocator: std.mem.Allocator,
    input_path: []const u8,
    preference: IoPreference,
    gzip_mode: GzipMode,
) !LoadedInput {
    if (std.mem.eql(u8, input_path, "-")) {
        const bytes = try readStdinAll(allocator);
        return .{
            .bytes = bytes,
            .mode = .buffered,
            .file_size = bytes.len,
            .file_kind = .file,
            .backing = .{ .owned = bytes },
        };
    }

    if (isGzipPath(input_path)) {
        if (gzip_mode == .stream) {
            const bytes = try decompressGzipToMemory(allocator, input_path);
            return .{
                .bytes = bytes,
                .mode = .buffered,
                .file_size = bytes.len,
                .file_kind = .file,
                .backing = .{ .owned = bytes },
            };
        }

        const temp_path = try decompressGzipToTempFile(allocator, input_path);
        var file = try std.fs.cwd().openFile(temp_path, .{});
        defer file.close();
        var loaded = try loadInputData(allocator, file, .mmap);
        loaded.temp_path = temp_path;
        return loaded;
    }

    var file = try std.fs.cwd().openFile(input_path, .{});
    defer file.close();
    return try loadInputData(allocator, file, preference);
}

fn isGzipPath(path: []const u8) bool {
    return std.mem.endsWith(u8, path, ".gz") or std.mem.endsWith(u8, path, ".bgz");
}

fn readStdinAll(allocator: std.mem.Allocator) ![]u8 {
    var file = std.fs.File.stdin();
    return file.readToEndAlloc(allocator, std.math.maxInt(usize));
}

fn decompressGzipToMemory(allocator: std.mem.Allocator, input_path: []const u8) ![]u8 {
    var in_file = try std.fs.cwd().openFile(input_path, .{});
    defer in_file.close();
    var in_buf: [64 * 1024]u8 = undefined;
    var in_reader = in_file.reader(&in_buf);
    var flate_buf: [std.compress.flate.max_window_len]u8 = undefined;
    var decompress = std.compress.flate.Decompress.init(&in_reader.interface, .gzip, &flate_buf);

    var out: std.Io.Writer.Allocating = .init(allocator);
    errdefer out.deinit();
    _ = try decompress.reader.streamRemaining(&out.writer);
    return try out.toOwnedSlice();
}

fn decompressGzipToTempFile(allocator: std.mem.Allocator, input_path: []const u8) ![]u8 {
    const tmp_path = try std.fmt.allocPrint(allocator, "/tmp/zdash-{d}.fastq", .{std.time.nanoTimestamp()});
    errdefer allocator.free(tmp_path);

    var in_file = try std.fs.cwd().openFile(input_path, .{});
    defer in_file.close();
    var out_file = try std.fs.cwd().createFile(tmp_path, .{ .exclusive = true });
    defer out_file.close();
    var in_buf: [64 * 1024]u8 = undefined;
    var in_reader = in_file.reader(&in_buf);
    var out_buf: [64 * 1024]u8 = undefined;
    var out_writer = out_file.writer(&out_buf);
    var flate_buf: [std.compress.flate.max_window_len]u8 = undefined;
    var decompress = std.compress.flate.Decompress.init(&in_reader.interface, .gzip, &flate_buf);
    _ = try decompress.reader.streamRemaining(&out_writer.interface);
    try out_writer.interface.flush();

    return tmp_path;
}

const RecordSpan = struct {
    start: usize,
    end: usize,
};

const RepairSummary = struct {
    written_reads: u64,
    rejected_reads: u64,
    output_path: []const u8,
    bad_records_path: ?[]const u8,
};

const SampleSummary = struct {
    total_reads_seen: u64,
    sampled_reads: u64,
    output_path: ?[]const u8,
    seed: u64,
    fraction: ?f64,
    n: ?usize,
};

fn runSample(allocator: std.mem.Allocator, data: []const u8, run: RunOptions) !SampleSummary {
    var records = std.ArrayList(RecordSpan).empty;
    defer records.deinit(allocator);

    var cursor: usize = 0;
    while (cursor < data.len) {
        const start = cursor;
        const h = nextLineFast(data, cursor) orelse break;
        cursor = h.next;
        const s = nextLineFast(data, cursor) orelse break;
        cursor = s.next;
        const p = nextLineFast(data, cursor) orelse break;
        cursor = p.next;
        const q = nextLineFast(data, cursor) orelse break;
        cursor = q.next;
        try records.append(allocator, .{ .start = start, .end = cursor });
    }

    var prng = std.Random.DefaultPrng.init(run.sample_seed);
    const random = prng.random();
    var picked = std.ArrayList(RecordSpan).empty;
    defer picked.deinit(allocator);

    if (run.sample_fraction) |fraction| {
        for (records.items) |rec| {
            if (random.float(f64) < fraction) try picked.append(allocator, rec);
        }
    } else if (run.sample_n) |n| {
        const k = @min(n, records.items.len);
        try picked.ensureTotalCapacity(allocator, k);
        var i: usize = 0;
        while (i < records.items.len) : (i += 1) {
            const rec = records.items[i];
            if (i < k) {
                try picked.append(allocator, rec);
            } else {
                const j = random.uintLessThan(usize, i + 1);
                if (j < k) picked.items[j] = rec;
            }
        }
    }

    std.sort.insertion(RecordSpan, picked.items, {}, struct {
        fn less(_: void, a: RecordSpan, b: RecordSpan) bool {
            return a.start < b.start;
        }
    }.less);

    if (run.output_path) |output_path| {
        if (std.mem.eql(u8, output_path, "-")) {
            for (picked.items) |rec| try printStdout("{s}", .{data[rec.start..rec.end]});
        } else {
            var f = try std.fs.cwd().createFile(output_path, .{ .truncate = true });
            defer f.close();
            for (picked.items) |rec| try f.writeAll(data[rec.start..rec.end]);
        }
    } else {
        for (picked.items) |rec| try printStdout("{s}", .{data[rec.start..rec.end]});
    }

    return .{
        .total_reads_seen = records.items.len,
        .sampled_reads = picked.items.len,
        .output_path = run.output_path,
        .seed = run.sample_seed,
        .fraction = run.sample_fraction,
        .n = run.sample_n,
    };
}

fn runRepair(allocator: std.mem.Allocator, data: []const u8, run: RunOptions) !RepairSummary {
    const output_path = run.output_path orelse "zdash_repaired.fastq";
    var out = try std.fs.cwd().createFile(output_path, .{ .truncate = true });
    defer out.close();
    var bad_out: ?std.fs.File = null;
    if (run.emit_bad_records_path) |p| bad_out = try std.fs.cwd().createFile(p, .{ .truncate = true });
    defer if (bad_out) |*f| f.close();

    var written: u64 = 0;
    var rejected: u64 = 0;
    var cursor: usize = 0;
    while (cursor < data.len) {
        const start = cursor;
        const header = nextLineFast(data, cursor) orelse break;
        cursor = header.next;
        const seq = nextLineFast(data, cursor) orelse break;
        cursor = seq.next;
        const sep = nextLineFast(data, cursor) orelse break;
        cursor = sep.next;
        const qual = nextLineFast(data, cursor) orelse break;
        cursor = qual.next;
        const end = cursor;

        var ok = true;
        if (header.line.bytes.len == 0 or header.line.bytes[0] != '@') ok = false;
        if (sep.line.bytes.len == 0 or sep.line.bytes[0] != '+') ok = false;
        if (seq.line.bytes.len != qual.line.bytes.len) ok = false;
        if (ok) {
            const s = analyzeSequenceHot(seq.line.bytes, true);
            const q = analyzeQualityHot(qual.line.bytes, true);
            if (s.invalid_index != null or q.invalid_index != null) ok = false;
        }

        if (ok) {
            try out.writeAll(data[start..end]);
            written += 1;
        } else {
            rejected += 1;
            if (bad_out) |*f| try f.writeAll(data[start..end]);
            if (run.repair_mode == .truncate_to_last_good) break;
        }
    }

    _ = allocator;
    return .{
        .written_reads = written,
        .rejected_reads = rejected,
        .output_path = output_path,
        .bad_records_path = run.emit_bad_records_path,
    };
}

fn writeValidationReportJson(
    report_path: []const u8,
    input_path: []const u8,
    loaded: LoadedInput,
    outcome: ValidationOutcome,
    config: ProcessingConfig,
    operation: Operation,
) !void {
    const json = try renderJsonSummary(std.heap.page_allocator, input_path, loaded, outcome, config, operation);
    defer std.heap.page_allocator.free(json);
    var f = try std.fs.cwd().createFile(report_path, .{ .truncate = true });
    defer f.close();
    try f.writeAll(json);
    try f.writeAll("\n");
}

fn jsonGet(value: std.json.Value, key: []const u8) ?std.json.Value {
    if (value != .object) return null;
    return value.object.get(key);
}

fn jsonString(value: std.json.Value) ?[]const u8 {
    return switch (value) {
        .string => |s| s,
        else => null,
    };
}

fn jsonNumber(value: std.json.Value) ?f64 {
    return switch (value) {
        .integer => |n| @as(f64, @floatFromInt(n)),
        .float => |f| f,
        else => null,
    };
}

fn inferFailureHint(message: []const u8) []const u8 {
    if (std.mem.indexOf(u8, message, "truncated") != null) return "Likely incomplete upload/write or stream truncation. Re-fetch or regenerate the file.";
    if (std.mem.indexOf(u8, message, "invalid base") != null) return "Sequence has non-ACGTN bytes. Validate upstream conversion and text encoding.";
    if (std.mem.indexOf(u8, message, "invalid quality") != null) return "Quality line includes non-Phred+33 printable bytes. Check pipeline encoding assumptions.";
    if (std.mem.indexOf(u8, message, "length mismatch") != null) return "Sequence and quality lengths differ. Upstream read framing is corrupted.";
    if (std.mem.indexOf(u8, message, "header") != null or std.mem.indexOf(u8, message, "separator") != null) return "FASTQ framing markers (@/+) are malformed. Check record boundaries.";
    return "Inspect the first failing read and upstream preprocessing for record damage.";
}

fn runExplain(allocator: std.mem.Allocator, report_path: []const u8, emit_json: bool) !void {
    const bytes = try std.fs.cwd().readFileAlloc(allocator, report_path, 8 * 1024 * 1024);
    defer allocator.free(bytes);
    var parsed = try std.json.parseFromSlice(std.json.Value, allocator, bytes, .{});
    defer parsed.deinit();

    const root = parsed.value;
    const status = if (jsonGet(root, "status")) |v| jsonString(v) orelse "unknown" else "unknown";
    const operation = if (jsonGet(root, "operation")) |v| jsonString(v) orelse "unknown" else "unknown";
    const schema = if (jsonGet(root, "schema_version")) |v| jsonString(v) orelse "unknown" else "unknown";
    const input = if (jsonGet(root, "input")) |v| jsonString(v) orelse "unknown" else "unknown";

    if (!std.mem.eql(u8, status, "failed")) {
        if (emit_json) {
            const payload = .{
                .status = "ok",
                .operation = operation,
                .schema_version = schema,
                .input = input,
                .summary = "Report indicates success; no validation failures to explain.",
            };
            try printStdout("{f}\n", .{std.json.fmt(payload, .{})});
        } else {
            try printStdout("Explain: report status is ok (operation={s}, input={s}, schema={s}).\n", .{ operation, input, schema });
        }
        return;
    }

    const err_obj = jsonGet(root, "error") orelse {
        return error.InvalidJson;
    };
    const message = if (jsonGet(err_obj, "message")) |v| jsonString(v) orelse "unknown error" else "unknown error";
    const line = if (jsonGet(err_obj, "line")) |v| jsonString(v) orelse "unknown" else "unknown";
    const read_index = if (jsonGet(err_obj, "read_index")) |v| jsonNumber(v) orelse 0 else 0;
    const byte_offset = if (jsonGet(err_obj, "byte_offset")) |v| jsonNumber(v) orelse 0 else 0;
    const hint = inferFailureHint(message);

    if (emit_json) {
        const payload = .{
            .status = "failed",
            .operation = operation,
            .schema_version = schema,
            .input = input,
            .failure = .{
                .message = message,
                .line = line,
                .read_index = read_index,
                .byte_offset = byte_offset,
            },
            .hint = hint,
        };
        try printStdout("{f}\n", .{std.json.fmt(payload, .{})});
    } else {
        try printStdout(
            "Explain ({s}):\n- Input: {s}\n- Schema: {s}\n- Failure: {s}\n- Line: {s}\n- Read: {d}\n- Byte: {d}\n- Hint: {s}\n",
            .{ operation, input, schema, message, line, @as(u64, @intFromFloat(read_index)), @as(u64, @intFromFloat(byte_offset)), hint },
        );
    }
}

fn runCompare(allocator: std.mem.Allocator, before_path: []const u8, after_path: []const u8, emit_json: bool) !void {
    const before_bytes = try std.fs.cwd().readFileAlloc(allocator, before_path, 8 * 1024 * 1024);
    defer allocator.free(before_bytes);
    const after_bytes = try std.fs.cwd().readFileAlloc(allocator, after_path, 8 * 1024 * 1024);
    defer allocator.free(after_bytes);

    var before_parsed = try std.json.parseFromSlice(std.json.Value, allocator, before_bytes, .{});
    defer before_parsed.deinit();
    var after_parsed = try std.json.parseFromSlice(std.json.Value, allocator, after_bytes, .{});
    defer after_parsed.deinit();

    const b = before_parsed.value;
    const a = after_parsed.value;

    const b_status = if (jsonGet(b, "status")) |v| jsonString(v) orelse "unknown" else "unknown";
    const a_status = if (jsonGet(a, "status")) |v| jsonString(v) orelse "unknown" else "unknown";
    const b_reads = if (jsonGet(b, "total_reads")) |v| jsonNumber(v) orelse 0 else 0;
    const a_reads = if (jsonGet(a, "total_reads")) |v| jsonNumber(v) orelse 0 else 0;
    const b_bases = if (jsonGet(b, "total_bases")) |v| jsonNumber(v) orelse 0 else 0;
    const a_bases = if (jsonGet(a, "total_bases")) |v| jsonNumber(v) orelse 0 else 0;
    const b_gc = if (jsonGet(b, "gc_percent")) |v| jsonNumber(v) orelse 0 else 0;
    const a_gc = if (jsonGet(a, "gc_percent")) |v| jsonNumber(v) orelse 0 else 0;
    const b_q = if (jsonGet(b, "avg_quality")) |v| jsonNumber(v) orelse 0 else 0;
    const a_q = if (jsonGet(a, "avg_quality")) |v| jsonNumber(v) orelse 0 else 0;

    if (emit_json) {
        const payload = .{
            .before = .{ .path = before_path, .status = b_status },
            .after = .{ .path = after_path, .status = a_status },
            .delta = .{
                .total_reads = a_reads - b_reads,
                .total_bases = a_bases - b_bases,
                .gc_percent = a_gc - b_gc,
                .avg_quality = a_q - b_q,
            },
        };
        try printStdout("{f}\n", .{std.json.fmt(payload, .{})});
        return;
    }

    try printStdout(
        "Compare:\n- Before: {s} (status={s})\n- After:  {s} (status={s})\n- Delta reads: {d}\n- Delta bases: {d}\n- Delta GC%%: {d:.4}\n- Delta AvgQ: {d:.4}\n",
        .{
            before_path,
            b_status,
            after_path,
            a_status,
            @as(i64, @intFromFloat(a_reads - b_reads)),
            @as(i64, @intFromFloat(a_bases - b_bases)),
            a_gc - b_gc,
            a_q - b_q,
        },
    );
}

fn loadInputData(allocator: std.mem.Allocator, file: std.fs.File, preference: IoPreference) !LoadedInput {
    const stat = try file.stat();
    const size_u64 = stat.size;
    const mmap_supported = builtin.os.tag != .windows and stat.kind == .file and size_u64 > 0;

    if (preference == .buffered) {
        const owned = try readEntireFileFallback(allocator, file, size_u64);
            return .{
                .bytes = owned,
                .mode = .buffered,
                .file_size = size_u64,
                .file_kind = stat.kind,
                .backing = .{ .owned = owned },
            };
        }

    if (mmap_supported) {
        const map_len = std.math.cast(usize, size_u64) orelse return error.FileTooLarge;
        const mapped = std.posix.mmap(
            null,
            map_len,
            std.posix.PROT.READ,
            .{ .TYPE = .PRIVATE },
            file.handle,
            0,
        ) catch |err| {
            if (preference == .mmap) return err;
            const owned = try readEntireFileFallback(allocator, file, size_u64);
            return .{
                .bytes = owned,
                .mode = .buffered,
                .file_size = size_u64,
                .file_kind = stat.kind,
                .backing = .{ .owned = owned },
            };
        };
        // Hint the OS that we will scan this mapping linearly.
        std.posix.madvise(mapped.ptr, mapped.len, std.posix.MADV.SEQUENTIAL) catch {};
        return .{
            .bytes = mapped,
            .mode = .mmap,
            .file_size = size_u64,
            .file_kind = stat.kind,
            .backing = .{ .mapped = mapped },
        };
    }

    if (preference == .mmap) return error.MemoryMappingNotSupported;

    const owned = try readEntireFileFallback(allocator, file, size_u64);
    return .{
        .bytes = owned,
        .mode = .buffered,
        .file_size = size_u64,
        .file_kind = stat.kind,
        .backing = .{ .owned = owned },
    };
}

fn readEntireFileFallback(allocator: std.mem.Allocator, file: std.fs.File, size_hint_u64: u64) ![]u8 {
    try file.seekTo(0);
    const max_bytes = if (size_hint_u64 > 0)
        (std.math.cast(usize, size_hint_u64) orelse return error.FileTooLarge)
    else
        std.math.maxInt(usize);
    return file.readToEndAlloc(allocator, max_bytes);
}

fn validateFastq(
    allocator: std.mem.Allocator,
    data: []const u8,
    config: ProcessingConfig,
    show_progress: bool,
) !ValidationResult {
    return validateFastqWithChunkSizeAndThreads(
        allocator,
        data,
        config.chunk_bytes,
        config.threads,
        config.profile,
        config.mode,
        show_progress,
    );
}

fn validateFastqDetailed(
    allocator: std.mem.Allocator,
    data: []const u8,
    config: ProcessingConfig,
    show_progress: bool,
    max_errors: usize,
    extract_error_context_path: ?[]const u8,
    error_window_reads: u64,
    extract_error_debug_path: ?[]const u8,
) !ValidationOutcome {
    if (max_errors <= 1 and extract_error_context_path == null and extract_error_debug_path == null) {
        const validation = try validateFastq(allocator, data, config, show_progress);
        return .{ .validation = validation };
    }

    var errors = std.ArrayList(ValidationError).empty;
    defer errors.deinit(allocator);
    var stats = Stats{};
    collectValidationErrors(allocator, data, 0, config.profile, config.mode, max_errors, &stats, &errors);

    if (errors.items.len == 0) return .{ .validation = .{ .ok = stats } };

    var out = ValidationOutcome{
        .validation = .{ .fail = errors.items[0] },
        .errors = try allocator.dupe(ValidationError, errors.items),
    };
    if (extract_error_context_path) |path| {
        out.context_output_path = try writeErrorContextFastq(allocator, data, errors.items[0].read_index, path, error_window_reads);
    }
    if (extract_error_debug_path) |path| {
        try writeErrorDebugReport(data, errors.items[0], path, out.context_output_path);
    }
    return out;
}

fn collectValidationErrors(
    allocator: std.mem.Allocator,
    data: []const u8,
    read_index_offset: u64,
    profile: WorkProfile,
    mode: RunMode,
    max_errors: usize,
    stats: *Stats,
    errors: *std.ArrayList(ValidationError),
) void {
    var scanner = BatchedLineScanner.init(data);
    const validate_chars = profile != .stats_only and mode == .strict;
    const track_low_quality = profile == .full;
    const q_threshold: u64 = @as(u64, @intFromFloat(LOW_QUALITY_AVG_Q_THRESHOLD));
    var records_seen: u64 = 0;

    while (scanner.line_start < data.len and errors.items.len < max_errors) {
        scanner.fillTo(4);
        if (scanner.available() < 4) {
            const read_index = read_index_offset + records_seen + 1;
            _ = errors.append(allocator, .{
                .read_index = read_index,
                .byte_offset = data.len,
                .issue = .{ .truncated_record = .{ .expected_line = 2 } },
            }) catch {};
            break;
        }

        records_seen += 1;
        const absolute_read_index = read_index_offset + records_seen;
        const header_start = scanner.line_start;
        const nl0 = scanner.peekNewline(0);
        const nl1 = scanner.peekNewline(1);
        const nl2 = scanner.peekNewline(2);
        const nl3 = scanner.peekNewline(3);
        const sequence_start = nl0 + 1;
        const plus_start = nl1 + 1;
        const quality_start = nl2 + 1;

        var had_error = false;
        if (header_start >= nl0 or data[header_start] != '@') {
            _ = errors.append(allocator, .{
                .read_index = absolute_read_index,
                .byte_offset = header_start,
                .issue = .invalid_header_prefix,
            }) catch {};
            had_error = true;
        }
        if (!had_error and (plus_start >= nl2 or data[plus_start] != '+')) {
            _ = errors.append(allocator, .{
                .read_index = absolute_read_index,
                .byte_offset = plus_start,
                .issue = .invalid_separator_prefix,
            }) catch {};
            had_error = true;
        }

        var sequence_end = nl1;
        if (sequence_end > sequence_start and data[sequence_end - 1] == '\r') sequence_end -= 1;
        var quality_end = nl3;
        if (quality_end > quality_start and data[quality_end - 1] == '\r') quality_end -= 1;
        const sequence = data[sequence_start..sequence_end];
        const quality = data[quality_start..quality_end];

        if (!had_error and sequence.len != quality.len) {
            _ = errors.append(allocator, .{
                .read_index = absolute_read_index,
                .byte_offset = quality_start,
                .issue = .{ .sequence_quality_length_mismatch = .{ .sequence_len = sequence.len, .quality_len = quality.len } },
            }) catch {};
            had_error = true;
        }

        if (!had_error) {
            const seq_kernel = analyzeSequenceHot(sequence, validate_chars);
            if (validate_chars and seq_kernel.invalid_index != null) {
                const invalid_index = seq_kernel.invalid_index.?;
                _ = errors.append(allocator, .{
                    .read_index = absolute_read_index,
                    .byte_offset = sequence_start + invalid_index,
                    .issue = .{ .invalid_sequence_char = .{ .ch = sequence[invalid_index] } },
                }) catch {};
                had_error = true;
            }

            if (!had_error) {
                const qual_kernel = analyzeQualityHot(quality, validate_chars);
                if (validate_chars and qual_kernel.invalid_index != null) {
                    const invalid_index = qual_kernel.invalid_index.?;
                    _ = errors.append(allocator, .{
                        .read_index = absolute_read_index,
                        .byte_offset = quality_start + invalid_index,
                        .issue = .{ .invalid_quality_char = .{ .ch = quality[invalid_index] } },
                    }) catch {};
                    had_error = true;
                } else {
                    stats.gc_bases += seq_kernel.gc_bases;
                    stats.total_bases += sequence.len;
                    stats.quality_sum += qual_kernel.quality_sum;
                    if (track_low_quality and quality.len > 0 and qual_kernel.quality_sum < q_threshold * quality.len) {
                        stats.low_quality_reads += 1;
                    }
                    stats.total_reads += 1;
                }
            }
        }
        scanner.consumeRecord();
    }
}

fn validateFastqWithChunkSizeAndThreads(
    allocator: std.mem.Allocator,
    data: []const u8,
    target_chunk_bytes: usize,
    requested_threads: usize,
    profile: WorkProfile,
    mode: RunMode,
    show_progress: bool,
) !ValidationResult {
    const safe_target_chunk_bytes = if (target_chunk_bytes == 0) @as(usize, 1) else target_chunk_bytes;
    if (data.len == 0) return .{ .ok = .{} };
    const worker_count = @max(@as(usize, 1), requested_threads);
    if (worker_count == 1) {
        return processChunksSingleThread(data, profile, mode, show_progress);
    }
    return try processChunksMultiThread(allocator, data, worker_count, safe_target_chunk_bytes, profile, mode, show_progress);
}

fn processChunksSingleThread(
    data: []const u8,
    profile: WorkProfile,
    mode: RunMode,
    show_progress: bool,
) ValidationResult {
    const result = validateFastqChunk(data, 0, 0, profile, mode);
    if (show_progress) {
        var progress = ProgressState{ .enabled = true, .total = data.len };
        addProgress(&progress, data.len);
    }
    return result;
}

fn processChunksMultiThread(
    allocator: std.mem.Allocator,
    data: []const u8,
    worker_count: usize,
    target_chunk_bytes: usize,
    profile: WorkProfile,
    mode: RunMode,
    show_progress: bool,
) !ValidationResult {
    var results = try allocator.alloc(WorkerResult, worker_count);
    defer allocator.free(results);
    @memset(results, .{});

    var contexts = try allocator.alloc(WorkerContext, worker_count);
    defer allocator.free(contexts);

    var threads = try allocator.alloc(std.Thread, worker_count - 1);
    defer allocator.free(threads);
    var queue = try ChunkQueue.init(allocator, @max(@as(usize, 64), worker_count * WORK_CHUNK_BATCH * 8));
    defer queue.deinit(allocator);
    var producer: ?std.Thread = null;

    var progress = ProgressState{ .enabled = true, .total = data.len };
    const progress_ptr: ?*ProgressState = if (show_progress) &progress else null;
    for (0..worker_count) |i| {
        contexts[i] = .{
            .data = data,
            .queue = &queue,
            .profile = profile,
            .mode = mode,
            .result = &results[i],
            .progress = progress_ptr,
        };
    }

    for (1..worker_count) |i| {
        threads[i - 1] = try std.Thread.spawn(.{}, workerMain, .{&contexts[i]});
    }
    var producer_ctx = ProducerContext{
        .data = data,
        .target_chunk_bytes = target_chunk_bytes,
        .queue = &queue,
    };
    producer = try std.Thread.spawn(.{}, producerMain, .{&producer_ctx});
    workerMain(&contexts[0]);
    for (threads) |thread| thread.join();
    if (producer) |p| p.join();

    var merged_stats = Stats{};
    var first_error: ?ValidationError = null;
    for (results) |result| {
        addStats(&merged_stats, result.stats);
        if (result.first_error) |err| {
            first_error = pickFirstError(first_error, err);
        }
    }

    if (first_error) |err| return .{ .fail = err };
    return .{ .ok = merged_stats };
}

fn workerMain(context: *WorkerContext) void {
    var local = WorkerResult{};
    if (context.queue) |queue| {
        while (queue.pop()) |chunk| {
            const chunk_slice = context.data[chunk.start..chunk.end];
            const result = validateFastqChunk(chunk_slice, chunk.read_index_base, chunk.start, context.profile, context.mode);
            switch (result) {
                .ok => |chunk_stats| addStats(&local.stats, chunk_stats),
                .fail => |err| local.first_error = pickFirstError(local.first_error, err),
            }
            if (context.progress) |progress| addProgress(progress, chunk_slice.len);
        }
    }
    context.result.* = local;
}

fn producerMain(context: *ProducerContext) void {
    const data = context.data;
    const target = if (context.target_chunk_bytes == 0) @as(usize, 1) else context.target_chunk_bytes;
    var chunk_start: usize = 0;
    var chunk_read_index_base: u64 = 0;
    var cursor: usize = 0;
    var line_in_record: u8 = 0;
    var reads_so_far: u64 = 0;
    var last_record_boundary: ?usize = null;
    var last_record_boundary_reads: u64 = 0;

    while (cursor < data.len) {
        while (cursor < data.len and data[cursor] != '\n') : (cursor += 1) {}
        if (cursor < data.len) cursor += 1;

        line_in_record += 1;
        if (line_in_record == 4) {
            line_in_record = 0;
            reads_so_far += 1;
            last_record_boundary = cursor;
            last_record_boundary_reads = reads_so_far;
        }

        if (cursor - chunk_start >= target) {
            if (last_record_boundary) |boundary| {
                if (boundary > chunk_start and boundary < data.len) {
                    if (chunk_start < data.len) @prefetch(data.ptr + chunk_start, .{});
                    context.queue.push(.{
                        .start = chunk_start,
                        .end = boundary,
                        .read_index_base = chunk_read_index_base,
                    });
                    chunk_start = boundary;
                    chunk_read_index_base = last_record_boundary_reads;
                    last_record_boundary = null;
                }
            }
        }
    }

    if (chunk_start < data.len) {
        if (chunk_start < data.len) @prefetch(data.ptr + chunk_start, .{});
        context.queue.push(.{
            .start = chunk_start,
            .end = data.len,
            .read_index_base = chunk_read_index_base,
        });
    }
    context.queue.markDone();
}

fn addStats(total: *Stats, chunk_stats: Stats) void {
    total.total_reads += chunk_stats.total_reads;
    total.total_bases += chunk_stats.total_bases;
    total.gc_bases += chunk_stats.gc_bases;
    total.quality_sum += chunk_stats.quality_sum;
    total.low_quality_reads += chunk_stats.low_quality_reads;
}

fn pickFirstError(existing: ?ValidationError, candidate: ValidationError) ValidationError {
    if (existing) |prior| {
        if (candidate.byte_offset < prior.byte_offset) return candidate;
        if (candidate.byte_offset == prior.byte_offset and candidate.read_index < prior.read_index) return candidate;
        return prior;
    }
    return candidate;
}

fn buildChunks(allocator: std.mem.Allocator, data: []const u8, target_chunk_bytes: usize) ![]Chunk {
    var chunks = std.ArrayList(Chunk).empty;
    defer chunks.deinit(allocator);

    if (data.len == 0) return try chunks.toOwnedSlice(allocator);

    var chunk_start: usize = 0;
    var chunk_read_index_base: u64 = 0;
    var cursor: usize = 0;
    var line_in_record: u8 = 0;
    var reads_so_far: u64 = 0;
    var last_record_boundary: ?usize = null;
    var last_record_boundary_reads: u64 = 0;

    while (cursor < data.len) {
        while (cursor < data.len and data[cursor] != '\n') : (cursor += 1) {}
        if (cursor < data.len) cursor += 1;

        line_in_record += 1;
        if (line_in_record == 4) {
            line_in_record = 0;
            reads_so_far += 1;
            last_record_boundary = cursor;
            last_record_boundary_reads = reads_so_far;
        }

        if (cursor - chunk_start >= target_chunk_bytes) {
            if (last_record_boundary) |boundary| {
                if (boundary > chunk_start and boundary < data.len) {
                    try chunks.append(allocator, .{
                        .start = chunk_start,
                        .end = boundary,
                        .read_index_base = chunk_read_index_base,
                    });
                    chunk_start = boundary;
                    chunk_read_index_base = last_record_boundary_reads;
                    last_record_boundary = null;
                }
            }
        }
    }

    if (chunk_start < data.len) {
        try chunks.append(allocator, .{
            .start = chunk_start,
            .end = data.len,
            .read_index_base = chunk_read_index_base,
        });
    }

    return try chunks.toOwnedSlice(allocator);
}

fn addProgress(progress: *ProgressState, delta_bytes: usize) void {
    if (!progress.enabled or progress.total == 0) return;
    progress.mutex.lock();
    defer progress.mutex.unlock();
    progress.processed += delta_bytes;
    const percent = @min((progress.processed * 100) / progress.total, 100);
    if (percent == progress.last_percent and progress.processed != progress.total) return;
    progress.last_percent = percent;
    var buffer: [256]u8 = undefined;
    var writer = std.fs.File.stderr().writer(&buffer);
    const out = &writer.interface;
    out.print("\r[progress] {d: >3}% ({d}/{d} bytes)", .{ percent, progress.processed, progress.total }) catch return;
    if (progress.processed == progress.total) {
        out.print("\n", .{}) catch return;
    }
    out.flush() catch return;
}

fn validateFastqChunk(
    chunk: []const u8,
    read_index_offset: u64,
    byte_offset_base: usize,
    profile: WorkProfile,
    mode: RunMode,
) ValidationResult {
    const _prof_start = profileStart();
    defer profileEnd(.validate_chunk, _prof_start);
    return switch (mode) {
        .strict => validateFastqChunkStrict(chunk, read_index_offset, byte_offset_base, profile),
        .assume_valid => validateFastqChunkAssumeValid(chunk, read_index_offset, byte_offset_base, profile),
    };
}

fn validateFastqChunkStrict(
    chunk: []const u8,
    read_index_offset: u64,
    byte_offset_base: usize,
    profile: WorkProfile,
) ValidationResult {
    return switch (profile) {
        .full => validateFastqChunkStrictImpl(true, true, chunk, read_index_offset, byte_offset_base),
        .validate_stats => validateFastqChunkStrictImpl(true, false, chunk, read_index_offset, byte_offset_base),
        .stats_only => validateFastqChunkStrictImpl(false, false, chunk, read_index_offset, byte_offset_base),
    };
}

fn validateFastqChunkStrictImpl(
    comptime validate_chars: bool,
    comptime track_low_quality: bool,
    chunk: []const u8,
    read_index_offset: u64,
    byte_offset_base: usize,
) ValidationResult {
    var stats = Stats{};
    const q_threshold: u64 = @as(u64, @intFromFloat(LOW_QUALITY_AVG_Q_THRESHOLD));
    var scanner = BatchedLineScanner.init(chunk);
    const newline_batch_target: usize = 1024;

    while (scanner.line_start < chunk.len) {
        scanner.fillTo(newline_batch_target);
        if (scanner.available() < 4) {
            return validateFastqChunkStrictTailImpl(validate_chars, track_low_quality,
                chunk,
                scanner.line_start,
                read_index_offset,
                byte_offset_base,
                q_threshold,
                &stats,
            );
        }

        while (scanner.available() >= 4) {
            const frame_prof_enabled = g_profiler.enabled;
            const frame_start = if (frame_prof_enabled) std.time.nanoTimestamp() else 0;
            const absolute_read_index = read_index_offset + stats.total_reads + 1;
            const header_start = scanner.line_start;
            const nl0 = scanner.peekNewline(0);
            const nl1 = scanner.peekNewline(1);
            const nl2 = scanner.peekNewline(2);
            const nl3 = scanner.peekNewline(3);
            const sequence_start = nl0 + 1;
            const plus_start = nl1 + 1;
            const quality_start = nl2 + 1;

            if (header_start >= nl0 or chunk[header_start] != '@') {
                return .{ .fail = .{
                    .read_index = absolute_read_index,
                    .byte_offset = byte_offset_base + header_start,
                    .issue = .invalid_header_prefix,
                } };
            }
            if (plus_start >= nl2 or chunk[plus_start] != '+') {
                return .{ .fail = .{
                    .read_index = absolute_read_index,
                    .byte_offset = byte_offset_base + plus_start,
                    .issue = .invalid_separator_prefix,
                } };
            }

            var sequence_end = nl1;
            if (sequence_end > sequence_start and chunk[sequence_end - 1] == '\r') sequence_end -= 1;
            var quality_end = nl3;
            if (quality_end > quality_start and chunk[quality_end - 1] == '\r') quality_end -= 1;
            const sequence = chunk[sequence_start..sequence_end];
            const quality = chunk[quality_start..quality_end];

            if (sequence.len != quality.len) {
                return .{ .fail = .{
                    .read_index = absolute_read_index,
                    .byte_offset = byte_offset_base + quality_start,
                    .issue = .{
                        .sequence_quality_length_mismatch = .{
                            .sequence_len = sequence.len,
                            .quality_len = quality.len,
                        },
                    },
                } };
            }

            if (scanner.available() >= 8) {
                const prefetch_start = scanner.peekNewline(4) + 1;
                if (prefetch_start < chunk.len) @prefetch(chunk.ptr + prefetch_start, .{});
            }
            if (frame_prof_enabled) profileEnd(.strict_framing, frame_start);

            const seq_kernel = analyzeSequenceHot(sequence, validate_chars);
            if (validate_chars and seq_kernel.invalid_index != null) {
                const invalid_index = seq_kernel.invalid_index.?;
                return .{ .fail = .{
                    .read_index = absolute_read_index,
                    .byte_offset = byte_offset_base + sequence_start + invalid_index,
                    .issue = .{ .invalid_sequence_char = .{ .ch = sequence[invalid_index] } },
                } };
            }
            stats.gc_bases += seq_kernel.gc_bases;
            stats.total_bases += sequence.len;

            const qual_kernel = analyzeQualityHot(quality, validate_chars);
            if (validate_chars and qual_kernel.invalid_index != null) {
                const invalid_index = qual_kernel.invalid_index.?;
                return .{ .fail = .{
                    .read_index = absolute_read_index,
                    .byte_offset = byte_offset_base + quality_start + invalid_index,
                    .issue = .{ .invalid_quality_char = .{ .ch = quality[invalid_index] } },
                } };
            }
            stats.quality_sum += qual_kernel.quality_sum;

            if (track_low_quality and quality.len > 0 and qual_kernel.quality_sum < q_threshold * quality.len) {
                stats.low_quality_reads += 1;
            }

            stats.total_reads += 1;
            scanner.consumeRecord();
        }
    }

    return .{ .ok = stats };
}

fn validateFastqChunkStrictTailImpl(
    comptime validate_chars: bool,
    comptime track_low_quality: bool,
    chunk: []const u8,
    start: usize,
    read_index_offset: u64,
    byte_offset_base: usize,
    q_threshold: u64,
    stats: *Stats,
) ValidationResult {
    var cursor = start;
    while (true) {
        const header = nextLineFast(chunk, cursor) orelse break;
        cursor = header.next;
        const absolute_read_index = read_index_offset + stats.total_reads + 1;
        if (header.line.bytes.len == 0 or header.line.bytes[0] != '@') {
            return .{ .fail = .{
                .read_index = absolute_read_index,
                .byte_offset = byte_offset_base + header.line.start,
                .issue = .invalid_header_prefix,
            } };
        }

        const sequence = nextLineFast(chunk, cursor) orelse {
            return .{ .fail = .{
                .read_index = absolute_read_index,
                .byte_offset = byte_offset_base + chunk.len,
                .issue = .{ .truncated_record = .{ .expected_line = 2 } },
            } };
        };
        cursor = sequence.next;

        const separator = nextLineFast(chunk, cursor) orelse {
            return .{ .fail = .{
                .read_index = absolute_read_index,
                .byte_offset = byte_offset_base + chunk.len,
                .issue = .{ .truncated_record = .{ .expected_line = 3 } },
            } };
        };
        cursor = separator.next;
        if (separator.line.bytes.len == 0 or separator.line.bytes[0] != '+') {
            return .{ .fail = .{
                .read_index = absolute_read_index,
                .byte_offset = byte_offset_base + separator.line.start,
                .issue = .invalid_separator_prefix,
            } };
        }

        const quality = nextLineFast(chunk, cursor) orelse {
            return .{ .fail = .{
                .read_index = absolute_read_index,
                .byte_offset = byte_offset_base + chunk.len,
                .issue = .{ .truncated_record = .{ .expected_line = 4 } },
            } };
        };
        cursor = quality.next;

        if (sequence.line.bytes.len != quality.line.bytes.len) {
            return .{ .fail = .{
                .read_index = absolute_read_index,
                .byte_offset = byte_offset_base + quality.line.start,
                .issue = .{
                    .sequence_quality_length_mismatch = .{
                        .sequence_len = sequence.line.bytes.len,
                        .quality_len = quality.line.bytes.len,
                    },
                },
            } };
        }

        const seq_kernel = analyzeSequenceHot(sequence.line.bytes, validate_chars);
        if (validate_chars and seq_kernel.invalid_index != null) {
            const invalid_index = seq_kernel.invalid_index.?;
            return .{ .fail = .{
                .read_index = absolute_read_index,
                .byte_offset = byte_offset_base + sequence.line.start + invalid_index,
                .issue = .{ .invalid_sequence_char = .{ .ch = sequence.line.bytes[invalid_index] } },
            } };
        }
        stats.gc_bases += seq_kernel.gc_bases;
        stats.total_bases += sequence.line.bytes.len;

        const qual_kernel = analyzeQualityHot(quality.line.bytes, validate_chars);
        if (validate_chars and qual_kernel.invalid_index != null) {
            const invalid_index = qual_kernel.invalid_index.?;
            return .{ .fail = .{
                .read_index = absolute_read_index,
                .byte_offset = byte_offset_base + quality.line.start + invalid_index,
                .issue = .{ .invalid_quality_char = .{ .ch = quality.line.bytes[invalid_index] } },
            } };
        }
        stats.quality_sum += qual_kernel.quality_sum;
        if (track_low_quality and quality.line.bytes.len > 0 and qual_kernel.quality_sum < q_threshold * quality.line.bytes.len) {
            stats.low_quality_reads += 1;
        }
        stats.total_reads += 1;
    }

    return .{ .ok = stats.* };
}

fn validateFastqChunkAssumeValid(
    chunk: []const u8,
    read_index_offset: u64,
    byte_offset_base: usize,
    profile: WorkProfile,
) ValidationResult {
    return switch (profile) {
        .full => validateFastqChunkAssumeValidImpl(true, chunk, read_index_offset, byte_offset_base),
        .validate_stats => validateFastqChunkAssumeValidImpl(false, chunk, read_index_offset, byte_offset_base),
        .stats_only => validateFastqChunkAssumeValidImpl(false, chunk, read_index_offset, byte_offset_base),
    };
}

fn validateFastqChunkAssumeValidImpl(
    comptime track_low_quality: bool,
    chunk: []const u8,
    read_index_offset: u64,
    byte_offset_base: usize,
) ValidationResult {
    var stats = Stats{};
    const q_threshold: u64 = @as(u64, @intFromFloat(LOW_QUALITY_AVG_Q_THRESHOLD));
    var scanner = BatchedLineScanner.init(chunk);
    const newline_batch_target: usize = 1024;

    while (scanner.line_start < chunk.len) {
        scanner.fillTo(newline_batch_target);
        if (scanner.available() < 4) {
            const read_index = read_index_offset + stats.total_reads + 1;
            return .{ .fail = .{
                .read_index = read_index,
                .byte_offset = byte_offset_base + chunk.len,
                .issue = .{ .truncated_record = .{ .expected_line = 2 } },
            } };
        }

        while (scanner.available() >= 4) {
            const frame_prof_enabled = g_profiler.enabled;
            const frame_start = if (frame_prof_enabled) std.time.nanoTimestamp() else 0;
            const header_start = scanner.line_start;
            const nl0 = scanner.peekNewline(0);
            const nl1 = scanner.peekNewline(1);
            const nl2 = scanner.peekNewline(2);
            const nl3 = scanner.peekNewline(3);
            const sequence_start = nl0 + 1;
            const plus_start = nl1 + 1;
            const quality_start = nl2 + 1;

            if (header_start >= nl0 or chunk[header_start] != '@') {
                const read_index = read_index_offset + stats.total_reads + 1;
                return .{ .fail = .{
                    .read_index = read_index,
                    .byte_offset = byte_offset_base + header_start,
                    .issue = .invalid_header_prefix,
                } };
            }
            if (plus_start >= nl2 or chunk[plus_start] != '+') {
                const read_index = read_index_offset + stats.total_reads + 1;
                return .{ .fail = .{
                    .read_index = read_index,
                    .byte_offset = byte_offset_base + plus_start,
                    .issue = .invalid_separator_prefix,
                } };
            }

            var sequence_end = nl1;
            if (sequence_end > sequence_start and chunk[sequence_end - 1] == '\r') sequence_end -= 1;
            var quality_end = nl3;
            if (quality_end > quality_start and chunk[quality_end - 1] == '\r') quality_end -= 1;
            const sequence = chunk[sequence_start..sequence_end];
            const quality = chunk[quality_start..quality_end];

            if (sequence.len != quality.len) {
                const read_index = read_index_offset + stats.total_reads + 1;
                return .{ .fail = .{
                    .read_index = read_index,
                    .byte_offset = byte_offset_base + quality_start,
                    .issue = .{
                        .sequence_quality_length_mismatch = .{
                            .sequence_len = sequence.len,
                            .quality_len = quality.len,
                        },
                    },
                } };
            }

            if (scanner.available() >= 8) {
                const prefetch_start = scanner.peekNewline(4) + 1;
                if (prefetch_start < chunk.len) @prefetch(chunk.ptr + prefetch_start, .{});
            }
            if (frame_prof_enabled) profileEnd(.assume_framing, frame_start);

            const fused = analyzeRecordFused(sequence, quality);
            stats.total_reads += 1;
            stats.total_bases += sequence.len;
            stats.gc_bases += fused.gc_bases;
            stats.quality_sum += fused.quality_sum;
            if (track_low_quality and quality.len > 0 and fused.quality_sum < q_threshold * quality.len) {
                stats.low_quality_reads += 1;
            }

            scanner.consumeRecord();
        }
    }

    return .{ .ok = stats };
}

inline fn analyzeRecordFused(sequence: []const u8, quality: []const u8) FusedRecordResult {
    if (SIMD_AVAILABLE and sequence.len >= SIMD_WIDTH and quality.len >= SIMD_WIDTH) {
        return analyzeRecordFusedSimd(sequence, quality);
    }
    var gc: u64 = 0;
    var qsum: u64 = 0;
    for (sequence, quality) |base, q| {
        gc += GC_BASE_TABLE[base];
        qsum += DECODED_QUALITY_TABLE[q];
    }
    return .{ .gc_bases = gc, .quality_sum = qsum };
}

fn analyzeRecordFusedSimd(sequence: []const u8, quality: []const u8) FusedRecordResult {
    if (comptime !SIMD_AVAILABLE) return analyzeRecordFused(sequence, quality);
    const VecU8 = @Vector(SIMD_WIDTH, u8);
    const VecU16 = @Vector(SIMD_WIDTH, u16);
    const one: VecU8 = @splat(1);
    const zero: VecU8 = @splat(0);
    const vC: VecU8 = @splat(@as(u8, 'C'));
    const vG: VecU8 = @splat(@as(u8, 'G'));
    const vc: VecU8 = @splat(@as(u8, 'c'));
    const vg: VecU8 = @splat(@as(u8, 'g'));
    const q_offset: VecU16 = @splat(33);

    var gc: u64 = 0;
    var qsum: u64 = 0;
    const n = @min(sequence.len, quality.len);
    var i: usize = 0;
    while (i + SIMD_WIDTH <= n) : (i += SIMD_WIDTH) {
        const s: VecU8 = loadVector(SIMD_WIDTH, sequence, i);
        const q: VecU8 = loadVector(SIMD_WIDTH, quality, i);

        const gc_mask = (s == vG) | (s == vg) | (s == vC) | (s == vc);
        const gc_counts: VecU8 = @select(u8, gc_mask, one, zero);
        gc += @as(u64, @reduce(.Add, gc_counts));

        const decoded: VecU16 = @as(VecU16, q) - q_offset;
        qsum += @as(u64, @reduce(.Add, decoded));
    }

    for (sequence[i..n], quality[i..n]) |base, q| {
        gc += GC_BASE_TABLE[base];
        qsum += DECODED_QUALITY_TABLE[q];
    }
    return .{ .gc_bases = gc, .quality_sum = qsum };
}

inline fn analyzeSequenceHot(sequence: []const u8, validate_chars: bool) SequenceKernelResult {
    if (validate_chars) {
        if (g_profiler.enabled) return analyzeSequenceKernel(sequence);
        if (SIMD_AVAILABLE and sequence.len >= SIMD_WIDTH) return analyzeSequenceSimd(sequence);
        return analyzeSequenceScalar(sequence);
    }
    return analyzeSequenceGcOnly(sequence);
}

inline fn analyzeQualityHot(quality: []const u8, validate_chars: bool) QualityKernelResult {
    if (validate_chars) {
        if (g_profiler.enabled) return analyzeQualityKernel(quality);
        if (SIMD_AVAILABLE and quality.len >= SIMD_WIDTH) return analyzeQualitySimd(quality);
        return analyzeQualityScalar(quality);
    }
    return analyzeQualitySumOnly(quality);
}

fn nextLineFast(data: []const u8, start: usize) ?struct { line: Line, next: usize } {
    if (start >= data.len) return null;

    const maybe_newline = std.mem.indexOfScalarPos(u8, data, start, '\n');
    const raw_end = maybe_newline orelse data.len;
    var line_end = raw_end;
    if (line_end > start and data[line_end - 1] == '\r') line_end -= 1;

    return .{
        .line = .{
            .bytes = data[start..line_end],
            .start = start,
        },
        .next = if (maybe_newline != null) raw_end + 1 else raw_end,
    };
}

fn findChunkEnd(data: []const u8, start: usize, target_chunk_bytes: usize) usize {
    const remaining = data.len - start;
    if (remaining <= target_chunk_bytes) return data.len;

    var cursor = start;
    var line_in_record: u8 = 0;
    var last_record_boundary: ?usize = null;

    while (cursor < data.len) {
        while (cursor < data.len and data[cursor] != '\n') : (cursor += 1) {}
        if (cursor < data.len) cursor += 1;

        line_in_record += 1;
        if (line_in_record == 4) line_in_record = 0;
        if (line_in_record == 0) {
            last_record_boundary = cursor;
            if (cursor - start >= target_chunk_bytes) break;
        }
    }

    if (last_record_boundary) |boundary| {
        if (boundary > start) return boundary;
    }
    return data.len;
}

fn analyzeSequenceKernel(sequence: []const u8) SequenceKernelResult {
    const _prof_start = profileStart();
    defer profileEnd(.sequence_kernel, _prof_start);
    if (SIMD_AVAILABLE and sequence.len >= SIMD_WIDTH) {
        return analyzeSequenceSimd(sequence);
    }
    return analyzeSequenceScalar(sequence);
}

fn analyzeSequenceGcOnly(sequence: []const u8) SequenceKernelResult {
    if (SIMD_AVAILABLE and sequence.len >= SIMD_WIDTH) {
        return analyzeSequenceGcOnlySimd(sequence);
    }
    var out = SequenceKernelResult{};
    for (sequence) |ch| {
        const up = ch & 0xDF;
        out.gc_bases += @intFromBool(up == 'C' or up == 'G');
    }
    return out;
}

fn analyzeSequenceScalar(sequence: []const u8) SequenceKernelResult {
    var out = SequenceKernelResult{};
    for (sequence, 0..) |ch, i| {
        const up = ch & 0xDF;
        if (up != 'A' and up != 'C' and up != 'G' and up != 'T' and up != 'N') {
            out.invalid_index = i;
            return out;
        }
        out.gc_bases += @intFromBool(up == 'C' or up == 'G');
    }
    return out;
}

fn analyzeSequenceSimd(sequence: []const u8) SequenceKernelResult {
    if (comptime !SIMD_AVAILABLE) return analyzeSequenceScalar(sequence);

    const Vec = @Vector(SIMD_WIDTH, u8);
    const one: Vec = @splat(1);
    const zero: Vec = @splat(0);
    const case_mask: Vec = @splat(0xDF);
    const vA: Vec = @splat(@as(u8, 'A'));
    const vC: Vec = @splat(@as(u8, 'C'));
    const vG: Vec = @splat(@as(u8, 'G'));
    const vT: Vec = @splat(@as(u8, 'T'));
    const vN: Vec = @splat(@as(u8, 'N'));

    var out = SequenceKernelResult{};
    var i: usize = 0;
    while (i + SIMD_WIDTH <= sequence.len) : (i += SIMD_WIDTH) {
        const vec: Vec = loadVector(SIMD_WIDTH, sequence, i);
        const upper = vec & case_mask;
        const valid_mask = (upper == vA) | (upper == vC) | (upper == vG) | (upper == vT) | (upper == vN);
        if (!@reduce(.And, valid_mask)) {
            var j: usize = 0;
            while (j < SIMD_WIDTH) : (j += 1) {
                const ch = sequence[i + j];
                const up = ch & 0xDF;
                if (up != 'A' and up != 'C' and up != 'G' and up != 'T' and up != 'N') {
                    out.invalid_index = i + j;
                    return out;
                }
                out.gc_bases += @intFromBool(up == 'C' or up == 'G');
            }
            continue;
        }

        const gc_mask = (upper == vG) | (upper == vC);
        const gc_counts: Vec = @select(u8, gc_mask, one, zero);
        out.gc_bases += @as(u64, @reduce(.Add, gc_counts));
    }

    const tail = analyzeSequenceScalar(sequence[i..]);
    if (tail.invalid_index) |idx| out.invalid_index = i + idx;
    out.gc_bases += tail.gc_bases;
    return out;
}

fn analyzeSequenceGcOnlySimd(sequence: []const u8) SequenceKernelResult {
    if (comptime !SIMD_AVAILABLE) return analyzeSequenceGcOnly(sequence);
    const Vec = @Vector(SIMD_WIDTH, u8);
    const one: Vec = @splat(1);
    const zero: Vec = @splat(0);
    const case_mask: Vec = @splat(0xDF);
    const vC: Vec = @splat(@as(u8, 'C'));
    const vG: Vec = @splat(@as(u8, 'G'));

    var out = SequenceKernelResult{};
    var i: usize = 0;
    while (i + SIMD_WIDTH <= sequence.len) : (i += SIMD_WIDTH) {
        const vec: Vec = loadVector(SIMD_WIDTH, sequence, i);
        const upper = vec & case_mask;
        const gc_mask = (upper == vG) | (upper == vC);
        const gc_counts: Vec = @select(u8, gc_mask, one, zero);
        out.gc_bases += @as(u64, @reduce(.Add, gc_counts));
    }
    for (sequence[i..]) |ch| {
        const up = ch & 0xDF;
        out.gc_bases += @intFromBool(up == 'C' or up == 'G');
    }
    return out;
}

fn analyzeQualityKernel(quality: []const u8) QualityKernelResult {
    const _prof_start = profileStart();
    defer profileEnd(.quality_kernel, _prof_start);
    if (SIMD_AVAILABLE and quality.len >= SIMD_WIDTH) {
        return analyzeQualitySimd(quality);
    }
    return analyzeQualityScalar(quality);
}

fn analyzeQualitySumOnly(quality: []const u8) QualityKernelResult {
    if (SIMD_AVAILABLE and quality.len >= SIMD_WIDTH) {
        return analyzeQualitySumOnlySimd(quality);
    }
    var out = QualityKernelResult{};
    for (quality) |ch| out.quality_sum += DECODED_QUALITY_TABLE[ch];
    return out;
}

fn analyzeQualityScalar(quality: []const u8) QualityKernelResult {
    var out = QualityKernelResult{};
    for (quality, 0..) |ch, i| {
        if (!VALID_QUALITY_TABLE[ch]) {
            out.invalid_index = i;
            return out;
        }
        out.quality_sum += DECODED_QUALITY_TABLE[ch];
    }
    return out;
}

fn analyzeQualitySimd(quality: []const u8) QualityKernelResult {
    if (comptime !SIMD_AVAILABLE) return analyzeQualityScalar(quality);

    const VecU8 = @Vector(SIMD_WIDTH, u8);
    const VecU16 = @Vector(SIMD_WIDTH, u16);
    const v_min: VecU8 = @splat(33);
    const v_max: VecU8 = @splat(126);
    const q_offset: VecU16 = @splat(33);

    var out = QualityKernelResult{};
    var i: usize = 0;
    while (i + SIMD_WIDTH <= quality.len) : (i += SIMD_WIDTH) {
        const vec: VecU8 = loadVector(SIMD_WIDTH, quality, i);
        const valid_mask = (vec >= v_min) & (vec <= v_max);
        if (!@reduce(.And, valid_mask)) {
            var j: usize = 0;
            while (j < SIMD_WIDTH) : (j += 1) {
                const ch = quality[i + j];
                if (!VALID_QUALITY_TABLE[ch]) {
                    out.invalid_index = i + j;
                    return out;
                }
                out.quality_sum += DECODED_QUALITY_TABLE[ch];
            }
            continue;
        }

        const decoded: VecU16 = @as(VecU16, vec) - q_offset;
        out.quality_sum += @as(u64, @reduce(.Add, decoded));
    }

    const tail = analyzeQualityScalar(quality[i..]);
    if (tail.invalid_index) |idx| out.invalid_index = i + idx;
    out.quality_sum += tail.quality_sum;
    return out;
}

fn analyzeQualitySumOnlySimd(quality: []const u8) QualityKernelResult {
    if (comptime !SIMD_AVAILABLE) return analyzeQualitySumOnly(quality);
    const VecU8 = @Vector(SIMD_WIDTH, u8);
    const VecU16 = @Vector(SIMD_WIDTH, u16);
    const q_offset: VecU16 = @splat(33);

    var out = QualityKernelResult{};
    var i: usize = 0;
    while (i + SIMD_WIDTH <= quality.len) : (i += SIMD_WIDTH) {
        const vec: VecU8 = loadVector(SIMD_WIDTH, quality, i);
        const decoded: VecU16 = @as(VecU16, vec) - q_offset;
        out.quality_sum += @as(u64, @reduce(.Add, decoded));
    }
    for (quality[i..]) |ch| out.quality_sum += DECODED_QUALITY_TABLE[ch];
    return out;
}

fn loadVector(comptime lanes: usize, bytes: []const u8, start: usize) @Vector(lanes, u8) {
    const ptr: *align(1) const @Vector(lanes, u8) = @ptrCast(bytes.ptr + start);
    return ptr.*;
}

fn readLine(data: []const u8, cursor: *usize) ?Line {
    const _prof_start = profileStart();
    defer profileEnd(.newline_refill, _prof_start);
    if (cursor.* >= data.len) return null;

    const start = cursor.*;
    var end = cursor.*;
    while (end < data.len and data[end] != '\n') : (end += 1) {}

    cursor.* = if (end < data.len) end + 1 else end;

    var line_end = end;
    if (line_end > start and data[line_end - 1] == '\r') {
        line_end -= 1;
    }

    return .{
        .bytes = data[start..line_end],
        .start = start,
    };
}

fn isValidSequenceBase(ch: u8) bool {
    return VALID_BASE_TABLE[ch];
}

fn buildValidBaseTable() [256]bool {
    var table = [_]bool{false} ** 256;
    table['A'] = true;
    table['C'] = true;
    table['G'] = true;
    table['T'] = true;
    table['N'] = true;
    table['a'] = true;
    table['c'] = true;
    table['g'] = true;
    table['t'] = true;
    table['n'] = true;
    return table;
}

fn buildGcBaseTable() [256]u8 {
    var table = [_]u8{0} ** 256;
    table['C'] = 1;
    table['G'] = 1;
    table['c'] = 1;
    table['g'] = 1;
    return table;
}

fn buildValidQualityTable() [256]bool {
    var table = [_]bool{false} ** 256;
    var ch: usize = 33;
    while (ch <= 126) : (ch += 1) {
        table[ch] = true;
    }
    return table;
}

fn buildDecodedQualityTable() [256]u64 {
    var table = [_]u64{0} ** 256;
    var ch: usize = 33;
    while (ch <= 126) : (ch += 1) {
        table[ch] = ch - 33;
    }
    return table;
}

fn validationIssueMessage(issue: ValidationIssue) []const u8 {
    return switch (issue) {
        .truncated_record => |v| switch (v.expected_line) {
            2 => "record ended early (missing sequence line)",
            3 => "record ended early (missing separator line)",
            4 => "record ended early (missing quality line)",
            else => "record ended early",
        },
        .invalid_header_prefix => "header line must start with '@'",
        .invalid_separator_prefix => "separator line must start with '+'",
        .invalid_sequence_char => "sequence contains invalid base (allowed: A,C,G,T,N, lowercase)",
        .invalid_quality_char => "quality contains invalid ASCII (allowed range: 33..126)",
        .sequence_quality_length_mismatch => "sequence and quality line lengths do not match",
    };
}

fn validationIssueLine(issue: ValidationIssue) []const u8 {
    return switch (issue) {
        .invalid_header_prefix => "header",
        .invalid_sequence_char => "sequence",
        .invalid_separator_prefix => "separator",
        .invalid_quality_char => "quality",
        .sequence_quality_length_mismatch => "quality",
        .truncated_record => |v| switch (v.expected_line) {
            2 => "sequence",
            3 => "separator",
            4 => "quality",
            else => "unknown",
        },
    };
}

fn validationBadChar(issue: ValidationIssue) ?u8 {
    return switch (issue) {
        .invalid_sequence_char => |v| v.ch,
        .invalid_quality_char => |v| v.ch,
        else => null,
    };
}

fn printValidationFailureDetailed(data: []const u8, fail: ValidationError) !void {
    try printStderr("ERROR: {s}\n", .{validationIssueMessage(fail.issue)});
    try printStderr("Read: {d}\n", .{fail.read_index});
    try printStderr("Line: {s}\n", .{validationIssueLine(fail.issue)});
    try printStderr("Byte offset: {d}\n", .{fail.byte_offset});
    if (validationBadChar(fail.issue)) |bad| {
        if (std.ascii.isPrint(bad)) {
            try printStderr("Bad char: '{c}' (0x{x})\n", .{ bad, bad });
        } else {
            try printStderr("Bad char: 0x{x}\n", .{bad});
        }
    }

    if (lineContextAt(data, fail.byte_offset)) |ctx| {
        try printStderr("Context: {s}\n", .{ctx.line});
        if (ctx.caret > 0) {
            var i: usize = 0;
            while (i < ctx.caret) : (i += 1) try printStderr(" ", .{});
        }
        try printStderr("^\n", .{});
    }
}

fn printGithubAnnotation(input_path: []const u8, fail: ValidationError) !void {
    var title_buf: [128]u8 = undefined;
    const title = std.fmt.bufPrint(&title_buf, "FASTQ validation error (read {d})", .{fail.read_index}) catch "FASTQ validation error";
    var message_buf: [512]u8 = undefined;
    const message = std.fmt.bufPrint(&message_buf, "{s} (line={s}, byte={d})", .{
        validationIssueMessage(fail.issue),
        validationIssueLine(fail.issue),
        fail.byte_offset,
    }) catch "validation failed";
    var escaped_buf: [768]u8 = undefined;
    const escaped = fillEscapeGithubAnnotation(&escaped_buf, message);
    try printStderr(
        "::error file={s},line={d},title={s}::{s}\n",
        .{ input_path, fail.read_index, title, escaped },
    );
}

fn fillEscapeGithubAnnotation(out: []u8, input: []const u8) []const u8 {
    // Keep this lightweight: we only escape control markers used by GitHub command parser.
    var w: usize = 0;
    for (input) |ch| {
        if (w + 3 >= out.len) break;
        switch (ch) {
            '%' => {
                out[w] = '%';
                out[w + 1] = '2';
                out[w + 2] = '5';
                w += 3;
            },
            '\n' => {
                out[w] = '%';
                out[w + 1] = '0';
                out[w + 2] = 'A';
                w += 3;
            },
            '\r' => {
                out[w] = '%';
                out[w + 1] = '0';
                out[w + 2] = 'D';
                w += 3;
            },
            else => {
                out[w] = ch;
                w += 1;
            },
        }
    }
    return out[0..w];
}

fn lineContextAt(data: []const u8, byte_offset: usize) ?struct { line: []const u8, caret: usize } {
    if (byte_offset >= data.len) return null;
    var start = byte_offset;
    while (start > 0 and data[start - 1] != '\n') : (start -= 1) {}
    var end = byte_offset;
    while (end < data.len and data[end] != '\n') : (end += 1) {}
    if (end > start and data[end - 1] == '\r') end -= 1;
    return .{
        .line = data[start..end],
        .caret = byte_offset - start,
    };
}

fn writeErrorContextFastq(
    allocator: std.mem.Allocator,
    data: []const u8,
    read_index: u64,
    output_path: []const u8,
    window_reads: u64,
) ![]u8 {
    const target = if (read_index == 0) @as(u64, 1) else read_index;
    const start_read = if (target > window_reads) target - window_reads else 1;
    const end_read = target + window_reads;

    var cursor: usize = 0;
    var current_read: u64 = 0;
    var write_start: ?usize = null;
    var write_end: usize = data.len;

    while (cursor < data.len) {
        const r_start = cursor;
        const l1 = nextLineFast(data, cursor) orelse break;
        cursor = l1.next;
        const l2 = nextLineFast(data, cursor) orelse break;
        cursor = l2.next;
        const l3 = nextLineFast(data, cursor) orelse break;
        cursor = l3.next;
        const l4 = nextLineFast(data, cursor) orelse break;
        cursor = l4.next;

        current_read += 1;
        if (current_read == start_read) write_start = r_start;
        if (current_read == end_read) {
            write_end = cursor;
            break;
        }
    }

    const slice_start = write_start orelse 0;
    const slice = data[slice_start..write_end];
    var file = try std.fs.cwd().createFile(output_path, .{ .truncate = true });
    defer file.close();
    try file.writeAll(slice);
    return try allocator.dupe(u8, output_path);
}

fn writeErrorDebugReport(
    data: []const u8,
    fail: ValidationError,
    output_path: []const u8,
    context_output_path: ?[]const u8,
) !void {
    var file = try std.fs.cwd().createFile(output_path, .{ .truncate = true });
    defer file.close();
    var out_buf: [4096]u8 = undefined;
    var file_writer = file.writer(&out_buf);
    const writer = &file_writer.interface;

    try writer.print("Z-DASH Error Debug Report\n", .{});
    try writer.print("=========================\n", .{});
    try writer.print("message: {s}\n", .{validationIssueMessage(fail.issue)});
    try writer.print("line: {s}\n", .{validationIssueLine(fail.issue)});
    try writer.print("read_index: {d}\n", .{fail.read_index});
    try writer.print("byte_offset: {d}\n", .{fail.byte_offset});
    if (context_output_path) |path| try writer.print("context_fastq: {s}\n", .{path});
    if (validationBadChar(fail.issue)) |bad| {
        if (std.ascii.isPrint(bad)) {
            try writer.print("bad_char: '{c}' (0x{x})\n", .{ bad, bad });
        } else {
            try writer.print("bad_char: 0x{x}\n", .{bad});
        }
    }
    if (lineContextAt(data, fail.byte_offset)) |ctx| {
        try writer.print("line_context: {s}\n", .{ctx.line});
    }

    const radius: usize = 16;
    const start = fail.byte_offset -| radius;
    const end = @min(data.len, fail.byte_offset + radius + 1);
    try writer.print("hex_window_start: {d}\n", .{start});
    try writer.print("hex_window_end: {d}\n", .{end});
    try writer.print("hex_window: ", .{});
    for (data[start..end]) |b| {
        try writer.print("{X:0>2} ", .{b});
    }
    try writer.print("\n", .{});
    try writer.flush();
}

fn inputModeLabel(mode: InputMode) []const u8 {
    return switch (mode) {
        .mmap => "mmap",
        .buffered => "buffered",
    };
}

fn workProfileLabel(profile: WorkProfile) []const u8 {
    return switch (profile) {
        .full => "full",
        .validate_stats => "validate-stats",
        .stats_only => "stats-only",
    };
}

fn runModeLabel(mode: RunMode) []const u8 {
    return switch (mode) {
        .strict => "strict",
        .assume_valid => "assume-valid",
    };
}

fn operationLabel(op: Operation) []const u8 {
    return switch (op) {
        .check => "check",
        .scan => "scan",
        .stats => "stats",
        .repair => "repair",
        .sample => "sample",
        .explain => "explain",
        .compare => "compare",
    };
}

fn averageQuality(stats: Stats) f64 {
    if (stats.total_bases == 0) return 0.0;
    return @as(f64, @floatFromInt(stats.quality_sum)) / @as(f64, @floatFromInt(stats.total_bases));
}

fn gcPercent(stats: Stats) f64 {
    if (stats.total_bases == 0) return 0.0;
    return (@as(f64, @floatFromInt(stats.gc_bases)) * 100.0) / @as(f64, @floatFromInt(stats.total_bases));
}

fn sumQualityScores(quality: []const u8) u64 {
    var sum: u64 = 0;
    for (quality) |ch| {
        sum += @as(u64, ch) - 33;
    }
    return sum;
}

test "parser and stats on valid FASTQ" {
    const input =
        \\@r1
        \\ACGTN
        \\+
        \\IIIII
        \\@r2
        \\GGCC
        \\+
        \\!!!!
        \\
    ;

    const got = try validateFastq(
        std.testing.allocator,
        input,
        .{ .chunk_bytes = TARGET_CHUNK_BYTES, .threads = 1, .profile = .full, .mode = .strict },
        false,
    );
    switch (got) {
        .ok => |stats| {
            try std.testing.expectEqual(@as(u64, 2), stats.total_reads);
            try std.testing.expectEqual(@as(u64, 9), stats.total_bases);
            try std.testing.expectEqual(@as(u64, 6), stats.gc_bases);
            try std.testing.expectEqual(@as(u64, 200), stats.quality_sum);
            try std.testing.expectEqual(@as(u64, 1), stats.low_quality_reads);
            try std.testing.expectApproxEqAbs(@as(f64, 22.2222), averageQuality(stats), 0.001);
            try std.testing.expectApproxEqAbs(@as(f64, 66.6666), gcPercent(stats), 0.001);
        },
        .fail => try std.testing.expect(false),
    }
}

test "validation fails on invalid sequence character" {
    const input =
        \\@r1
        \\ACGX
        \\+
        \\IIII
        \\
    ;
    const got = try validateFastq(
        std.testing.allocator,
        input,
        .{ .chunk_bytes = TARGET_CHUNK_BYTES, .threads = 1, .profile = .full, .mode = .strict },
        false,
    );
    switch (got) {
        .ok => try std.testing.expect(false),
        .fail => |fail| {
            try std.testing.expectEqual(@as(u64, 1), fail.read_index);
            try std.testing.expectEqual(@as(usize, 7), fail.byte_offset);
            try std.testing.expect(fail.issue == .invalid_sequence_char);
        },
    }
}

test "validation fails on sequence quality length mismatch" {
    const input =
        \\@r1
        \\ACGT
        \\+
        \\III
        \\
    ;
    const got = try validateFastq(
        std.testing.allocator,
        input,
        .{ .chunk_bytes = TARGET_CHUNK_BYTES, .threads = 1, .profile = .full, .mode = .strict },
        false,
    );
    switch (got) {
        .ok => try std.testing.expect(false),
        .fail => |fail| {
            try std.testing.expect(fail.issue == .sequence_quality_length_mismatch);
        },
    }
}

test "validation fails on truncated record" {
    const input =
        \\@r1
        \\ACGT
        \\+
    ;
    const got = try validateFastq(
        std.testing.allocator,
        input,
        .{ .chunk_bytes = TARGET_CHUNK_BYTES, .threads = 1, .profile = .full, .mode = .strict },
        false,
    );
    switch (got) {
        .ok => try std.testing.expect(false),
        .fail => |fail| {
            try std.testing.expect(fail.issue == .truncated_record);
        },
    }
}

test "summary formatting golden output for pass case" {
    const stats = Stats{
        .total_reads = 2,
        .total_bases = 9,
        .gc_bases = 6,
        .quality_sum = 200,
        .low_quality_reads = 1,
    };
    const fake_loaded = LoadedInput{
        .bytes = "dummy-bytes",
        .mode = .buffered,
        .file_size = 1234,
        .file_kind = .file,
        .backing = .{ .owned = @constCast("dummy-bytes") },
    };
    const output = try renderSummary(
        std.testing.allocator,
        "sample.fastq",
        fake_loaded,
        .{ .ok = stats },
        .{ .chunk_bytes = TARGET_CHUNK_BYTES, .threads = 8, .profile = .full, .mode = .strict },
    );
    defer std.testing.allocator.free(output);

    const expected =
        \\Z-DASH: Zig DNA Shredder v0.1.0-dev
        \\-----------------------------
        \\Input: sample.fastq
        \\
        \\Summary:
        \\- Total Reads:  2
        \\- Avg Quality:  Q22.22
        \\- GC Content:   66.67%
        \\- Low Quality:  1 reads (avg Q < 20.00)
        \\- I/O Mode:     buffered
        \\- Profile:      full
        \\- Mode:         strict
        \\- Threads:      8
        \\- Chunk Bytes:  8388608
        \\- File Bytes:   1234
        \\- Proc Bytes:   11
        \\- Status:       PASSED (No corruption found)
        \\
    ;

    try std.testing.expectEqualStrings(expected, output);
}

test "chunk boundary aligns to FASTQ record starts" {
    const input =
        \\@r1
        \\ACGTN
        \\+
        \\IIIII
        \\@r2
        \\GGCC
        \\+
        \\!!!!
        \\
    ;
    const end = findChunkEnd(input, 0, 10);
    try std.testing.expectEqual(@as(usize, 18), end);
    var cursor = end;
    try std.testing.expectEqualStrings("@r2", readLine(input, &cursor).?.bytes);
}

test "validation works with tiny chunk size" {
    const input =
        \\@r1
        \\ACGTN
        \\+
        \\IIIII
        \\@r2
        \\GGCC
        \\+
        \\!!!!
        \\
    ;
    const got = try validateFastqWithChunkSizeAndThreads(std.testing.allocator, input, 1, 1, .full, .strict, false);
    switch (got) {
        .ok => |stats| {
            try std.testing.expectEqual(@as(u64, 2), stats.total_reads);
            try std.testing.expectEqual(@as(u64, 9), stats.total_bases);
        },
        .fail => try std.testing.expect(false),
    }
}

test "integration: mmap and buffered modes produce same stats" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const sample =
        \\@r1
        \\ACGTN
        \\+
        \\IIIII
        \\@r2
        \\GGCC
        \\+
        \\!!!!
        \\
    ;
    try tmp.dir.writeFile(.{ .sub_path = "input.fastq", .data = sample });

    var file_mmap = try tmp.dir.openFile("input.fastq", .{});
    defer file_mmap.close();
    const loaded_mmap = try loadInputData(std.testing.allocator, file_mmap, .mmap);
    defer loaded_mmap.deinit(std.testing.allocator);
    const mmap_result = try validateFastq(
        std.testing.allocator,
        loaded_mmap.bytes,
        .{ .chunk_bytes = TARGET_CHUNK_BYTES, .threads = 2, .profile = .full, .mode = .strict },
        false,
    );

    var file_buffered = try tmp.dir.openFile("input.fastq", .{});
    defer file_buffered.close();
    const loaded_buffered = try loadInputData(std.testing.allocator, file_buffered, .buffered);
    defer loaded_buffered.deinit(std.testing.allocator);
    const buffered_result = try validateFastq(
        std.testing.allocator,
        loaded_buffered.bytes,
        .{ .chunk_bytes = TARGET_CHUNK_BYTES, .threads = 2, .profile = .full, .mode = .strict },
        false,
    );

    switch (mmap_result) {
        .ok => |mmap_stats| switch (buffered_result) {
            .ok => |buffered_stats| {
                try std.testing.expectEqual(mmap_stats.total_reads, buffered_stats.total_reads);
                try std.testing.expectEqual(mmap_stats.total_bases, buffered_stats.total_bases);
                try std.testing.expectEqual(mmap_stats.gc_bases, buffered_stats.gc_bases);
                try std.testing.expectEqual(mmap_stats.quality_sum, buffered_stats.quality_sum);
                try std.testing.expectEqual(mmap_stats.low_quality_reads, buffered_stats.low_quality_reads);
            },
            .fail => try std.testing.expect(false),
        },
        .fail => try std.testing.expect(false),
    }
}

test "robustness: tiny, huge, and malformed files" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    // tiny
    const tiny =
        \\@r1
        \\A
        \\+
        \\I
        \\
    ;
    try tmp.dir.writeFile(.{ .sub_path = "tiny.fastq", .data = tiny });
    var tiny_file = try tmp.dir.openFile("tiny.fastq", .{});
    defer tiny_file.close();
    const tiny_loaded = try loadInputData(std.testing.allocator, tiny_file, .auto);
    defer tiny_loaded.deinit(std.testing.allocator);
    const tiny_result = try validateFastq(
        std.testing.allocator,
        tiny_loaded.bytes,
        .{ .chunk_bytes = TARGET_CHUNK_BYTES, .threads = 2, .profile = .full, .mode = .strict },
        false,
    );
    switch (tiny_result) {
        .ok => |stats| try std.testing.expectEqual(@as(u64, 1), stats.total_reads),
        .fail => try std.testing.expect(false),
    }

    // huge-ish synthetic workload
    const record = "@r\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n";
    const record_count: usize = 50_000;
    const total_len = record.len * record_count;
    const big_data = try std.testing.allocator.alloc(u8, total_len);
    defer std.testing.allocator.free(big_data);
    for (0..record_count) |i| {
        const off = i * record.len;
        @memcpy(big_data[off .. off + record.len], record);
    }
    try tmp.dir.writeFile(.{ .sub_path = "huge.fastq", .data = big_data });
    var huge_file = try tmp.dir.openFile("huge.fastq", .{});
    defer huge_file.close();
    const huge_loaded = try loadInputData(std.testing.allocator, huge_file, .auto);
    defer huge_loaded.deinit(std.testing.allocator);
    const huge_result = try validateFastqWithChunkSizeAndThreads(std.testing.allocator, huge_loaded.bytes, 1024, 8, .full, .strict, false);
    switch (huge_result) {
        .ok => |stats| {
            try std.testing.expectEqual(@as(u64, record_count), stats.total_reads);
            try std.testing.expectEqual(@as(u64, record_count * 16), stats.total_bases);
        },
        .fail => try std.testing.expect(false),
    }

    // malformed
    const malformed =
        \\@r1
        \\ACGT
        \\+
        \\III
        \\
    ;
    try tmp.dir.writeFile(.{ .sub_path = "bad.fastq", .data = malformed });
    var bad_file = try tmp.dir.openFile("bad.fastq", .{});
    defer bad_file.close();
    const bad_loaded = try loadInputData(std.testing.allocator, bad_file, .auto);
    defer bad_loaded.deinit(std.testing.allocator);
    const bad_result = try validateFastq(
        std.testing.allocator,
        bad_loaded.bytes,
        .{ .chunk_bytes = TARGET_CHUNK_BYTES, .threads = 2, .profile = .full, .mode = .strict },
        false,
    );
    switch (bad_result) {
        .ok => try std.testing.expect(false),
        .fail => |fail| try std.testing.expect(fail.issue == .sequence_quality_length_mismatch),
    }
}

test "deterministic totals across repeated parallel runs" {
    const input =
        \\@r1
        \\ACGTN
        \\+
        \\IIIII
        \\@r2
        \\GGCC
        \\+
        \\!!!!
        \\@r3
        \\TTTT
        \\+
        \\IIII
        \\
    ;

    var expected: ?Stats = null;
    for (0..20) |_| {
        const result = try validateFastqWithChunkSizeAndThreads(std.testing.allocator, input, 1, 8, .full, .strict, false);
        switch (result) {
            .ok => |stats| {
                if (expected) |e| {
                    try std.testing.expectEqual(e.total_reads, stats.total_reads);
                    try std.testing.expectEqual(e.total_bases, stats.total_bases);
                    try std.testing.expectEqual(e.gc_bases, stats.gc_bases);
                    try std.testing.expectEqual(e.quality_sum, stats.quality_sum);
                    try std.testing.expectEqual(e.low_quality_reads, stats.low_quality_reads);
                } else {
                    expected = stats;
                }
            },
            .fail => try std.testing.expect(false),
        }
    }
}

test "race-prone chunking: many threads with tiny chunks still succeeds" {
    const record = "@r\nACGT\n+\nIIII\n";
    const count: usize = 2000;
    const total_len = record.len * count;
    const data = try std.testing.allocator.alloc(u8, total_len);
    defer std.testing.allocator.free(data);
    for (0..count) |i| {
        const off = i * record.len;
        @memcpy(data[off .. off + record.len], record);
    }

    const result = try validateFastqWithChunkSizeAndThreads(std.testing.allocator, data, 1, 32, .full, .strict, false);
    switch (result) {
        .ok => |stats| {
            try std.testing.expectEqual(@as(u64, count), stats.total_reads);
            try std.testing.expectEqual(@as(u64, count * 4), stats.total_bases);
        },
        .fail => try std.testing.expect(false),
    }
}

test "SIMD sequence kernel parity with scalar" {
    var prng = std.Random.DefaultPrng.init(0x5eed1234);
    const random = prng.random();

    var buf: [512]u8 = undefined;
    const alphabet = [_]u8{ 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n', 'X', '?' };

    for (0..200) |_| {
        const len = random.intRangeLessThan(usize, 0, buf.len + 1);
        for (buf[0..len]) |*ch| {
            ch.* = alphabet[random.uintLessThan(usize, alphabet.len)];
        }
        const scalar = analyzeSequenceScalar(buf[0..len]);
        const simd = analyzeSequenceKernel(buf[0..len]);
        try std.testing.expectEqual(scalar.gc_bases, simd.gc_bases);
        try std.testing.expectEqual(scalar.invalid_index, simd.invalid_index);
    }
}

test "SIMD quality kernel parity with scalar" {
    var prng = std.Random.DefaultPrng.init(0xabcddcba);
    const random = prng.random();

    var buf: [512]u8 = undefined;
    for (0..200) |_| {
        const len = random.intRangeLessThan(usize, 0, buf.len + 1);
        for (buf[0..len]) |*ch| {
            const bucket = random.uintLessThan(u8, 10);
            ch.* = if (bucket < 8) @as(u8, 33) + random.uintLessThan(u8, 94) else random.int(u8);
        }
        const scalar = analyzeQualityScalar(buf[0..len]);
        const simd = analyzeQualityKernel(buf[0..len]);
        try std.testing.expectEqual(scalar.quality_sum, simd.quality_sum);
        try std.testing.expectEqual(scalar.invalid_index, simd.invalid_index);
    }
}

test "parse run options supports quiet and chunk bytes" {
    const args = [_][]const u8{ "zdash", "--quiet", "--json", "--threads=4", "--chunk-bytes=1024", "in.fastq" };
    const run = try parseRunOptions(&args);
    try std.testing.expect(run.quiet);
    try std.testing.expect(run.json);
    try std.testing.expectEqual(@as(usize, 4), run.threads);
    try std.testing.expectEqual(@as(usize, 1024), run.chunk_bytes);
    try std.testing.expectEqualStrings("in.fastq", run.input_path);
}

test "parse run options supports gzip mode" {
    const args = [_][]const u8{ "zdash", "--gzip-mode=temp", "in.fastq.gz" };
    const run = try parseRunOptions(&args);
    try std.testing.expect(run.gzip_mode == .temp);
}

test "parse run options supports max-errors and context path" {
    const args = [_][]const u8{ "zdash", "check", "--max-errors=5", "--extract-error-context=errs.fastq", "in.fastq" };
    const run = try parseRunOptions(&args);
    try std.testing.expectEqual(@as(usize, 5), run.max_errors);
    try std.testing.expectEqualStrings("errs.fastq", run.extract_error_context_path.?);
}

test "parse run options supports profile" {
    const args = [_][]const u8{ "zdash", "--profile=stats-only", "in.fastq" };
    const run = try parseRunOptions(&args);
    try std.testing.expect(run.profile == .stats_only);
}

test "parse run options supports mode" {
    const args = [_][]const u8{ "zdash", "--mode=assume-valid", "in.fastq" };
    const run = try parseRunOptions(&args);
    try std.testing.expect(run.mode == .assume_valid);
}

test "parse run options supports preset" {
    const args = [_][]const u8{ "zdash", "--preset=fast-scan", "in.fastq" };
    const run = try parseRunOptions(&args);
    try std.testing.expect(run.operation == .scan);
    try std.testing.expect(run.mode == .assume_valid);
    try std.testing.expect(run.profile == .validate_stats);
    try std.testing.expect(run.profile_explicit);
}

test "parse run options rejects unsupported json schema version" {
    const args = [_][]const u8{ "zdash", "--json", "--json-schema-version=9.9.9", "in.fastq" };
    try std.testing.expectError(error.InvalidUsage, parseRunOptions(&args));
}

test "parse run options requires --against for compare" {
    const args = [_][]const u8{ "zdash", "compare", "before.json" };
    try std.testing.expectError(error.InvalidUsage, parseRunOptions(&args));
}

test "parse run options supports compare with against and debug context options" {
    const args = [_][]const u8{
        "zdash",
        "compare",
        "--against=after.json",
        "--extract-error-debug=debug.txt",
        "--error-window-reads=9",
        "before.json",
    };
    const run = try parseRunOptions(&args);
    try std.testing.expect(run.operation == .compare);
    try std.testing.expectEqualStrings("after.json", run.compare_against_path.?);
    try std.testing.expectEqualStrings("debug.txt", run.extract_error_debug_path.?);
    try std.testing.expectEqual(@as(u64, 9), run.error_window_reads);
}

test "config option parser applies preset and booleans" {
    var run = RunOptions{ .input_path = "" };
    try applyConfigOption(&run, "preset", "qc-only");
    try applyConfigOption(&run, "json", "true");
    try applyConfigOption(&run, "gha_annotations", "true");
    try std.testing.expect(run.operation == .stats);
    try std.testing.expect(run.mode == .assume_valid);
    try std.testing.expect(run.profile == .stats_only);
    try std.testing.expect(run.json);
    try std.testing.expect(run.gha_annotations);
}

test "parse run options supports subcommand defaults" {
    const args = [_][]const u8{ "zdash", "stats", "in.fastq" };
    const run = try parseRunOptions(&args);
    try std.testing.expect(run.operation == .stats);
    const cfg = resolveProcessingConfig(run);
    try std.testing.expect(cfg.mode == .assume_valid);
    try std.testing.expect(cfg.profile == .stats_only);
}

test "parse run options supports sample and repair flags" {
    const sample_args = [_][]const u8{ "zdash", "sample", "--seed=7", "--n=10", "--output=out.fastq", "in.fastq" };
    const sample = try parseRunOptions(&sample_args);
    try std.testing.expect(sample.operation == .sample);
    try std.testing.expectEqual(@as(u64, 7), sample.sample_seed);
    try std.testing.expectEqual(@as(usize, 10), sample.sample_n.?);
    try std.testing.expectEqualStrings("out.fastq", sample.output_path.?);

    const repair_args = [_][]const u8{ "zdash", "repair", "--repair-mode=truncate-to-last-good", "--emit-bad-records=bad.fastq", "in.fastq" };
    const repair = try parseRunOptions(&repair_args);
    try std.testing.expect(repair.operation == .repair);
    try std.testing.expect(repair.repair_mode == .truncate_to_last_good);
    try std.testing.expectEqualStrings("bad.fastq", repair.emit_bad_records_path.?);
}

test "parse run options rejects sample json without file output" {
    const args = [_][]const u8{ "zdash", "sample", "--json", "--n=2", "in.fastq" };
    try std.testing.expectError(error.InvalidUsage, parseRunOptions(&args));
}

test "parse run options rejects sample json stdout output" {
    const args = [_][]const u8{ "zdash", "sample", "--json", "--n=2", "--output=-", "in.fastq" };
    try std.testing.expectError(error.InvalidUsage, parseRunOptions(&args));
}

test "validate-stats profile disables low-quality counting" {
    const input =
        \\@r1
        \\ACGT
        \\+
        \\!!!!
        \\
    ;
    const result = try validateFastqWithChunkSizeAndThreads(std.testing.allocator, input, 1024, 1, .validate_stats, .strict, false);
    switch (result) {
        .ok => |stats| try std.testing.expectEqual(@as(u64, 0), stats.low_quality_reads),
        .fail => try std.testing.expect(false),
    }
}

test "stats-only profile skips base and quality character validation" {
    const input =
        \\@r1
        \\XXXX
        \\+
        \\~~~~
        \\
    ;
    const strict_result = try validateFastqWithChunkSizeAndThreads(std.testing.allocator, input, 1024, 1, .full, .strict, false);
    try std.testing.expect(strict_result == .fail);

    const relaxed_result = try validateFastqWithChunkSizeAndThreads(std.testing.allocator, input, 1024, 1, .stats_only, .strict, false);
    try std.testing.expect(relaxed_result == .ok);
}

test "assume-valid mode skips base and quality character validation" {
    const input =
        \\@r1
        \\XXXX
        \\+
        \\~~~~
        \\
    ;
    const strict_result = try validateFastqWithChunkSizeAndThreads(std.testing.allocator, input, 1024, 1, .full, .strict, false);
    try std.testing.expect(strict_result == .fail);

    const assume_valid_result = try validateFastqWithChunkSizeAndThreads(std.testing.allocator, input, 1024, 1, .full, .assume_valid, false);
    try std.testing.expect(assume_valid_result == .ok);
}

test "max-errors collects multiple failures from fixture" {
    var bad_f = try std.fs.cwd().openFile("tests/fixtures/bad_two_errors.fastq", .{});
    defer bad_f.close();
    const bad = try bad_f.readToEndAlloc(std.testing.allocator, std.math.maxInt(usize));
    defer std.testing.allocator.free(bad);
    const cfg = ProcessingConfig{
        .chunk_bytes = TARGET_CHUNK_BYTES,
        .threads = 1,
        .profile = .full,
        .mode = .strict,
    };
    const outcome = try validateFastqDetailed(std.testing.allocator, bad, cfg, false, 2, null, 5, null);
    defer outcome.deinit(std.testing.allocator);
    try std.testing.expect(outcome.validation == .fail);
    try std.testing.expectEqual(@as(usize, 2), outcome.errors.len);
    try std.testing.expectEqual(@as(u64, 1), outcome.errors[0].read_index);
    try std.testing.expectEqual(@as(u64, 2), outcome.errors[1].read_index);
}

test "repair modes match golden fixtures" {
    var in_f = try std.fs.cwd().openFile("tests/fixtures/repair_input.fastq", .{});
    defer in_f.close();
    const input = try in_f.readToEndAlloc(std.testing.allocator, std.math.maxInt(usize));
    defer std.testing.allocator.free(input);

    var drop_f = try std.fs.cwd().openFile("tests/fixtures/repair_expected_drop.fastq", .{});
    defer drop_f.close();
    const expect_drop = try drop_f.readToEndAlloc(std.testing.allocator, std.math.maxInt(usize));
    defer std.testing.allocator.free(expect_drop);

    var trunc_f = try std.fs.cwd().openFile("tests/fixtures/repair_expected_truncate.fastq", .{});
    defer trunc_f.close();
    const expect_trunc = try trunc_f.readToEndAlloc(std.testing.allocator, std.math.maxInt(usize));
    defer std.testing.allocator.free(expect_trunc);

    const out_drop = try std.fmt.allocPrint(std.testing.allocator, "/tmp/zdash-repair-drop-{d}.fastq", .{std.time.nanoTimestamp()});
    defer std.testing.allocator.free(out_drop);
    defer std.fs.deleteFileAbsolute(out_drop) catch {};
    const out_trunc = try std.fmt.allocPrint(std.testing.allocator, "/tmp/zdash-repair-trunc-{d}.fastq", .{std.time.nanoTimestamp()});
    defer std.testing.allocator.free(out_trunc);
    defer std.fs.deleteFileAbsolute(out_trunc) catch {};

    _ = try runRepair(std.testing.allocator, input, .{
        .input_path = "ignored.fastq",
        .operation = .repair,
        .output_path = out_drop,
        .repair_mode = .drop_bad_records,
    });
    _ = try runRepair(std.testing.allocator, input, .{
        .input_path = "ignored.fastq",
        .operation = .repair,
        .output_path = out_trunc,
        .repair_mode = .truncate_to_last_good,
    });

    var f_drop = try std.fs.openFileAbsolute(out_drop, .{});
    defer f_drop.close();
    const got_drop = try f_drop.readToEndAlloc(std.testing.allocator, std.math.maxInt(usize));
    defer std.testing.allocator.free(got_drop);

    var f_trunc = try std.fs.openFileAbsolute(out_trunc, .{});
    defer f_trunc.close();
    const got_trunc = try f_trunc.readToEndAlloc(std.testing.allocator, std.math.maxInt(usize));
    defer std.testing.allocator.free(got_trunc);

    try std.testing.expectEqualStrings(expect_drop, got_drop);
    try std.testing.expectEqualStrings(expect_trunc, got_trunc);
}
