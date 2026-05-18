# Barena

[![Tests](https://github.com/quadram-institute-bioscience/barena/actions/workflows/tests.yml/badge.svg)](https://github.com/quadram-institute-bioscience/barena/actions/workflows/tests.yml)

![logo](assets/barena-kraken.svg)

**Barena** is a tool to filter FASTQ reads by taxon using Kraken2 classification output.

## Usage

### Basic filtering
Filter reads matching taxon 9606:
```bash
barena -1 reads_R1.fq.gz -2 reads_R2.fq.gz -k kraken.out -t 9606 -o output_dir
```

### Filtering with descendants
To include all children of a taxon, provide the `kraken2-inspect` output:
```bash
barena -1 reads.fq -k kraken.out -t 9606 --child -d inspect.txt -c human_reads.fq.gz
```

### Streaming from Kraken2
```bash
kraken2 --db DB reads.fq | barena -1 reads.fq -t 9606 -o human_only
```

### Splitting into multiple taxon outputs
Provide a two-column whitespace-separated file with `TaxID` and `OutputFile`:
```text
TaxID   OutputFile
626929  Bacteroides
1578    Lactobacillus.fastq.gz
9605    Homo_#.fq
```

```bash
barena -1 reads_R1.fq.gz -2 reads_R2.fq.gz -k kraken.out --split targets.tsv -o split_reads
```

Reads that match more than one split target are written to every matching output. With
`--child`, each split TaxID is expanded independently.

When the `OutputFile` basename has no dot, Barena treats it as a stem and adds
`_R1.fastq` / `_R2.fastq` for paired reads, or `.fastq` for single-end reads. If the
basename contains a dot, it is treated as an explicit filename or `#` template.

## Options
- `-1, --r1`: FASTQ file (R1)
- `-2, --r2`: FASTQ file (R2)
- `-k, --kraken`: Kraken2 raw output (default: stdin)
- `-d, --inspect-db`: Kraken2-inspect output (required for `--child`)
- `-t, --taxon`: Taxon ID to match (can be repeated)
- `--split`: TaxID/OutputFile table for per-taxon split outputs
- `--child`: Include descendants of the specified taxons
- `-o, --outdir`: Output directory; normal mode uses default classified/unclassified names, split mode roots `OutputFile` paths
