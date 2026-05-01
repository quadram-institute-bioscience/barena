# Package
version       = "0.1.1"
author        = "Andrea Telatin"
description   = "Filter FASTQ reads using Kraken2 classification"
license       = "MIT"
srcDir        = "src"
binDir        = "bin"
bin           = @["barena"]

# Dependencies
requires "nim >= 2.0.0"
requires "readfx >= 0.3.0"
requires "argparse >= 4.0.0"

task test, "Build and run the test suite":
  exec "nim c -d:release -r tests/test_barena.nim"
