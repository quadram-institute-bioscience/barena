---
title: Home
layout: home
nav_order: 1
---

# barena

**Filter FASTQ reads by taxon using Kraken2 classification output**

[![Tests](https://github.com/quadram-institute-bioscience/barena/actions/workflows/tests.yml/badge.svg)](https://github.com/quadram-institute-bioscience/barena/actions/workflows/tests.yml)
[![GitHub release](https://img.shields.io/github/v/release/quadram-institute-bioscience/barena)](https://github.com/quadram-institute-bioscience/barena/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

---

`barena` takes a FASTQ file (single-end or paired-end) and a Kraken2 raw output file, and splits reads into **classified** and **unclassified** streams based on taxonomy.

![barena logo]({{ site.baseurl }}/assets/barena-kraken.svg){: style="max-width: 480px; display: block; margin: 2rem auto;"}

## Key features

- **Fast** — written in Nim; consistently 4–20× faster than KrakenTools on real datasets
- **Streaming** — pipe Kraken2 output directly without writing intermediate files
- **Flexible output** — output directory, per-read-pair paths, or placeholder templates
- **Taxonomy-aware** — optionally include all descendants of a taxon with `--child`
- **Paired-end aware** — handles R1/R2 in lock-step

## Quick start

```bash
# Extract all human reads (taxon 9606) from paired Illumina data
barena -1 R1.fq.gz -2 R2.fq.gz -k kraken.out -t 9606 -o human_reads/
```

```bash
# Stream from Kraken2 on the fly
kraken2 --db /db reads.fq | barena -1 reads.fq -t 9606 -o filtered/
```

```bash
# Include all descendants of taxon 9606 (requires kraken2-inspect output)
barena -1 R1.fq.gz -2 R2.fq.gz -k kraken.out \
       -t 9606 --child -d inspect.txt -o human_reads/
```

---

{: .note }
barena assumes that the order of reads in the FASTQ file matches the order of lines in the Kraken2 output. This is always true when both files come from the same Kraken2 run.

## Compared to KrakenTools

| | barena | KrakenTools |
|---|:---:|:---:|
| Language | Nim | Python |
| ONT `.fastq` (~3.5 M reads) | ~12 s | ~255 s |
| Illumina paired `.fastq.gz` (~53 M reads) | ~49 s | ~1033 s |
| Streaming from Kraken2 | ✅ | ❌ |
| Descendant filtering | ✅ | ✅ |
| Paired-end | ✅ | ✅ |

See the [benchmarks page]({{ site.baseurl }}/benchmarks) for full details.
