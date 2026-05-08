---
title: Installation
layout: default
nav_order: 2
---

# Installation
{: .no_toc }

<details open markdown="block">
  <summary>Contents</summary>
  {: .text-delta }
- TOC
{:toc}
</details>

---

## Pre-built binaries

The simplest way to install barena is to download a pre-built binary from the [GitHub Releases page](https://github.com/quadram-institute-bioscience/barena/releases).

```bash
# Linux x86_64
curl -L https://github.com/quadram-institute-bioscience/barena/releases/latest/download/barena-linux-x86_64 \
     -o barena && chmod +x barena
sudo mv barena /usr/local/bin/
```

Verify the installation:

```bash
barena --version
```

---

## Build from source

barena is written in [Nim](https://nim-lang.org/). You need Nim ≥ 2.0 and the `nimble` package manager.

### Install Nim

```bash
# Using choosenim (recommended)
curl https://nim-lang.org/choosenim/init.sh -sSf | sh
```

### Clone and build

```bash
git clone https://github.com/quadram-institute-bioscience/barena.git
cd barena
nimble build -d:release
```

The compiled binary lands in `./bin/barena`. Copy it to a directory on your `$PATH`:

```bash
sudo cp bin/barena /usr/local/bin/
```

### Run the test suite

```bash
nimble test
```

---

## Dependencies

| Dependency | Version | Purpose |
|---|---|---|
| `nim` | ≥ 2.0.0 | Compiler |
| [`readfx`](https://github.com/telatin/readfx) | ≥ 0.3.0 | FASTQ/FASTA parsing |
| [`argparse`](https://github.com/iffy/nim-argparse) | ≥ 4.0.0 | CLI argument parsing |

These are declared in `barena.nimble` and fetched automatically by `nimble build`.

---

## Conda / Bioconda

{: .note }
A Bioconda package is planned. Watch the [GitHub repository](https://github.com/quadram-institute-bioscience/barena) for updates.
