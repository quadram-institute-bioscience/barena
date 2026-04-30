# Barena

![logo](assets/barena.png)

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

## Options
- `-1, --r1`: FASTQ file (R1)
- `-2, --r2`: FASTQ file (R2)
- `-k, --kraken`: Kraken2 raw output (default: stdin)
- `-d, --inspect-db`: Kraken2-inspect output (required for `--child`)
- `-t, --taxon`: Taxon ID to match (can be repeated)
- `--child`: Include descendants of the specified taxons
- `-o, --outdir`: Output directory for classified/unclassified reads
