#!/usr/bin/env python3
import argparse
from pathlib import Path


def parse_checksums(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        sha = parts[0]
        filename = parts[-1]
        out[filename] = sha
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Homebrew formula for zdash from release checksums.")
    parser.add_argument("--repo", required=True, help="GitHub repo slug, e.g. zdash/zdash")
    parser.add_argument("--version", required=True, help="Release version without leading v, e.g. 0.1.0")
    parser.add_argument("--checksums", required=True, help="Path to checksums.txt")
    parser.add_argument("--out", required=True, help="Output formula path")
    args = parser.parse_args()

    checksums = parse_checksums(Path(args.checksums))
    mac_name = f"zdash-{args.version}-macos-aarch64.tar.gz"
    linux_name = f"zdash-{args.version}-linux-x86_64.tar.gz"

    if mac_name not in checksums:
        raise SystemExit(f"missing checksum for {mac_name}")
    if linux_name not in checksums:
        raise SystemExit(f"missing checksum for {linux_name}")

    formula = f"""class Zdash < Formula
  desc "Trust-first FASTQ validator and stats tool"
  homepage "https://github.com/{args.repo}"
  version "{args.version}"
  license "MIT"

  on_macos do
    if Hardware::CPU.arm?
      url "https://github.com/{args.repo}/releases/download/v#{{version}}/{mac_name}"
      sha256 "{checksums[mac_name]}"
    end
  end

  on_linux do
    if Hardware::CPU.intel?
      url "https://github.com/{args.repo}/releases/download/v#{{version}}/{linux_name}"
      sha256 "{checksums[linux_name]}"
    end
  end

  def install
    bin.install "zdash"
    doc.install "README.md"
    doc.install "LICENSE"
  end

  test do
    (testpath/"s.fastq").write <<~EOS
      @r1
      ACGT
      +
      IIII
    EOS
    output = shell_output("#{{bin}}/zdash check --json #{{testpath}}/s.fastq")
    assert_match "\\"status\\":\\"ok\\"", output
  end
end
"""
    Path(args.out).write_text(formula)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
