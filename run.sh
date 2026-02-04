#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

export ZIG_GLOBAL_CACHE_DIR="${ZIG_GLOBAL_CACHE_DIR:-/tmp/zig-global-cache}"
mkdir -p "$ZIG_GLOBAL_CACHE_DIR"

zig build -Doptimize=ReleaseFast >/dev/null

BIN="$ROOT_DIR/zig-out/bin/zdash"

if [[ "${1:-}" == "doctor" ]]; then
  shift
  exec "$BIN" doctor "$@"
fi

INPUT="${1:-tests/fixtures/valid_small.fastq}"
if [[ $# -gt 0 ]]; then
  shift
fi

if [[ ! -f "$INPUT" && "$INPUT" != "-" ]]; then
  echo "input not found: $INPUT" >&2
  exit 2
fi

echo "Running: $BIN --ci $INPUT $*"
exec "$BIN" --ci "$INPUT" "$@"
