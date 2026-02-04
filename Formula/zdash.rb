class Zdash < Formula
  desc "Trust-first FASTQ validator and stats tool"
  homepage "https://github.com/AaronRai123/zdash"
  version "0.1.0"
  license "MIT"

  on_macos do
    if Hardware::CPU.arm?
      url "https://github.com/AaronRai123/zdash/releases/download/v#{version}/zdash-#{version}-macos-aarch64.tar.gz"
      sha256 "REPLACE_WITH_RELEASE_SHA256"
    end
  end

  on_linux do
    if Hardware::CPU.intel?
      url "https://github.com/AaronRai123/zdash/releases/download/v#{version}/zdash-#{version}-linux-x86_64.tar.gz"
      sha256 "REPLACE_WITH_RELEASE_SHA256"
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
    output = shell_output("#{bin}/zdash check --json #{testpath}/s.fastq")
    assert_match "\"status\":\"ok\"", output
  end
end
