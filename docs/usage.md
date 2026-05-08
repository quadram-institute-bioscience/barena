---
title: Usage
layout: default
nav_order: 3
---

# Usage
{: .no_toc }

<details open markdown="block">
  <summary>Contents</summary>
  {: .text-delta }
- TOC
{:toc}
</details>

---

## Synopsis

```
barena [options]
```

barena reads a FASTQ file (or R1/R2 pair) and a Kraken2 raw output file **in lock-step** — the _n_-th FASTQ record is matched to the _n_-th Kraken2 line. Reads matching the requested taxon(s) are written to the classified output; all others go to the unclassified output.

---

## Options reference

### Input

| Flag | Description |
|---|---|
| `-1`, `--r1 FILE` | FASTQ file — single-end reads, or R1 of a pair. **(required)** |
| `-2`, `--r2 FILE` | FASTQ file — R2 of a pair; enables paired-end mode. |
| `-k`, `--kraken FILE` | Kraken2 raw output file. Defaults to **stdin** if omitted. |
| `-d`, `--inspect-db FILE` | `kraken2-inspect` output. Required when `--child` is set. |

### Output

Specify either `-o` (directory mode) or individual path flags (`-c`/`-C`/`-u`/`-U`).

| Flag | Description |
|---|---|
| `-o`, `--outdir DIR` | Write `classified_R1.fastq.gz` etc. into this directory. Overrides all other output flags. |
| `-c`, `--class-out FILE` | Classified R1 path. Use `#` as a placeholder replaced with `R1`/`R2`. |
| `-C`, `--class-out-r2 FILE` | Classified R2 path (explicit, no placeholder needed). |
| `-u`, `--unclass-out FILE` | Unclassified R1 path. Use `#` as a placeholder. |
| `-U`, `--unclass-out-r2 FILE` | Unclassified R2 path (explicit). |

{: .note }
When `-o` is used, all output files are **gzip-compressed** (`.fastq.gz`). When using `-c`/`-u`, the compression is inferred from the file extension.

### Filter

| Flag | Description |
|---|---|
| `-t`, `--taxon ID` | Taxonomy ID to match. Repeat to match multiple taxa. |
| `--child` | Include all descendants of the given taxid(s). Requires `-d`. |
| `--include-unclassified` | Route Kraken2-unclassified reads (taxId = 0) to the classified output. |

### Behaviour

| Flag | Description |
|---|---|
| `--strict` | Abort if read names differ between the Kraken2 output and the FASTQ. |
| `--verbose` | Print taxonomy stats and progress to stderr every 20,000 reads. |
| `--version` | Print version and exit. |

---

## Output modes

### Directory mode (`-o`)

The easiest option. barena creates the directory if it does not exist and writes:

```
outdir/
  classified_R1.fastq.gz
  classified_R2.fastq.gz    ← paired-end only
  unclassified_R1.fastq.gz
  unclassified_R2.fastq.gz  ← paired-end only
```

```bash
barena -1 R1.fq.gz -2 R2.fq.gz -k kraken.out -t 9606 -o human/
```

### Classified-only (`-c`)

Write only the matching reads; discard non-matching reads.

```bash
barena -1 reads.fq -k kraken.out -t 9606 -c human.fq.gz
```

### Placeholder mode

Use `#` in the path and barena replaces it with `R1` / `R2`:

```bash
barena -1 R1.fq.gz -2 R2.fq.gz -k kraken.out -t 9606 -c human_#.fq.gz
# writes: human_R1.fq.gz  human_R2.fq.gz
```

### Separate R1 / R2 paths

```bash
barena -1 R1.fq.gz -2 R2.fq.gz -k kraken.out -t 9606 \
       -c human_R1.fq.gz -C human_R2.fq.gz
```

---

## Filtering modes

### Exact taxon match

Match only reads assigned to **exactly** the given taxid:

```bash
barena -1 R1.fq.gz -k kraken.out -t 9606 -c human.fq.gz
```

### Multiple taxa

```bash
barena -1 R1.fq.gz -k kraken.out -t 9606 -t 10090 -c mammalian.fq.gz
```

### Include descendants (`--child`)

Match the taxid and **all its taxonomic descendants**. Requires the `kraken2-inspect` output:

```bash
# Generate the inspect file once per database
kraken2-inspect --db /path/to/db > inspect.txt

barena -1 R1.fq.gz -k kraken.out -t 9606 --child -d inspect.txt -c human.fq.gz
```

### Keep unclassified reads

Route Kraken2-unclassified reads (taxId = 0) into the classified output:

```bash
barena -1 R1.fq.gz -k kraken.out -t 9606 --include-unclassified -c host_and_unknown.fq.gz
```

---

## Streaming from Kraken2

barena reads Kraken2 output from **stdin** when `-k` is omitted, enabling on-the-fly filtering without writing intermediate files:

```bash
kraken2 --db /db reads.fq | barena -1 reads.fq -t 9606 -c filtered.fq.gz
```

{: .warning }
In streaming mode the FASTQ file is read **twice** (once by Kraken2, once by barena). Ensure the file is seekable (i.e. not itself a stream).
