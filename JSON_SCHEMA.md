# Z-DASH JSON Schema Contract

This document defines the stability policy for `--json` and `--report-json`.

## Versioning

- Current schema version: `1.0.0`
- Field is emitted as: `"schema_version": "1.0.0"`

Compatibility policy:

- Patch (`1.0.x`): additive/non-breaking changes only.
- Minor (`1.x.0`): additive changes allowed; existing fields remain stable.
- Major (`x.0.0`): breaking changes allowed (field renames/removals/type changes).

Runtime contract check:

- `--json-schema-version <v>` requires an exact schema version match.
- Current supported value is `1.0.0`.
- Mismatches fail fast with usage error (`2`).

## Common Top-Level Fields

- `tool` (string)
- `version` (string, binary version)
- `schema_version` (string, JSON contract version)
- `operation` (string: `check|scan|stats|repair|sample|explain|compare`)
- `status` (string: `ok|failed`)

## Validation/Stats Success Payload

Additional fields:

- `input`, `io_mode`, `profile`, `mode`
- `threads`, `chunk_bytes`, `file_bytes`, `processed_bytes`
- `total_reads`, `total_bases`, `gc_bases`, `gc_percent`
- `avg_quality`, `low_quality_reads`

## Validation Failure Payload

Additional fields:

- `input`, `io_mode`, `profile`, `mode`
- `threads`, `chunk_bytes`, `file_bytes`
- `max_errors` (number of collected errors)
- `context_output_path` (string or null)
- `error` object:
  - `message`
  - `line`
  - `read_index`
  - `byte_offset`
  - `bad_char` (number or null)

## Repair Success Payload

Additional fields:

- `written_reads`
- `rejected_reads`
- `output_path`
- `bad_records_path`

## Sample Success Payload

Additional fields:

- `total_reads_seen`
- `sampled_reads`
- `seed`
- `fraction` (number or null)
- `n` (number or null)
- `output_path` (string path; required when `--json` is used)

## Changelog Discipline

- Record schema-visible changes in `JSON_SCHEMA_CHANGELOG.md` before release.
