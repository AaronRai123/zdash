#!/usr/bin/env bash
set -euo pipefail

REPO="${REPO:-zdash/zdash}"
PREFIX="${PREFIX:-$HOME/.local/bin}"

OS="$(uname -s)"
ARCH="$(uname -m)"

case "$OS" in
  Darwin) OS_TAG="macos" ;;
  Linux) OS_TAG="linux" ;;
  *)
    echo "Unsupported OS: $OS" >&2
    exit 1
    ;;
esac

case "$ARCH" in
  arm64|aarch64) ARCH_TAG="aarch64" ;;
  x86_64|amd64) ARCH_TAG="x86_64" ;;
  *)
    echo "Unsupported architecture: $ARCH" >&2
    exit 1
    ;;
esac

ASSET_SUFFIX="${OS_TAG}-${ARCH_TAG}"
API_URL="https://api.github.com/repos/${REPO}/releases/latest"

JSON="$(curl -fsSL "$API_URL")"

TAG="$(python3 -c 'import json,sys; print(json.load(sys.stdin)["tag_name"])' <<<"$JSON")"
ASSET_URL="$(python3 -c 'import json,sys
data=json.load(sys.stdin)
suffix=sys.argv[1]
for a in data["assets"]:
    if a["name"].endswith(f"{suffix}.tar.gz"):
        print(a["browser_download_url"])
        break
else:
    raise SystemExit(1)
' "$ASSET_SUFFIX" <<<"$JSON")"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

ARCHIVE="$TMP_DIR/zdash.tar.gz"
curl -fL "$ASSET_URL" -o "$ARCHIVE"
tar -xzf "$ARCHIVE" -C "$TMP_DIR"

mkdir -p "$PREFIX"
BIN_PATH="$(find "$TMP_DIR" -type f -name zdash | head -n1)"
install -m 0755 "$BIN_PATH" "$PREFIX/zdash"

echo "Installed zdash ${TAG} to $PREFIX/zdash"
echo "If needed: export PATH=\"$PREFIX:\$PATH\""
