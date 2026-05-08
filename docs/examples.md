---
title: Examples
layout: default
nav_order: 4
---

# Examples
{: .no_toc }

<details open markdown="block">
  <summary>Contents</summary>
  {: .text-delta }
- TOC
{:toc}
</details>

---

## Host depletion (Illumina paired-end)

Remove human reads (taxid 9606) from a metagenome, keeping only non-human reads:

```bash
barena \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -k sample.kraken \
  -t 9606 \
  -u non_human_#.fastq.gz
```

Output: `non_human_R1.fastq.gz` and `non_human_R2.fastq.gz`.

---

## Targeted extraction (ONT long reads)

Extract all reads classified as *Bacillus spizizenii* (taxid 96241):

```bash
barena \
  -1 zymo.fastq \
  -k zymo.kraken \
  -t 96241 \
  -c bacillus.fastq.gz
```

---

## Extract a whole bacterial order and its descendants

Use `--child` to capture all taxa under *Lachnospirales* (taxid 45405):

```bash
# One-time: generate inspect file for the database
kraken2-inspect --db /db --threads 8 > db_inspect.txt

barena \
  -1 gut_R1.fastq.gz \
  -2 gut_R2.fastq.gz \
  -k gut.kraken \
  -t 45405 --child \
  -d db_inspect.txt \
  -o lachnospirales_reads/
```

---

## Split classified and unclassified simultaneously

Save both matching (e.g. bacterial) and non-matching reads in one pass:

```bash
barena \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -k sample.kraken \
  -t 2 --child \
  -d inspect.txt \
  -c bacteria_#.fastq.gz \
  -u non_bacterial_#.fastq.gz
```

---

## Streaming pipeline

Avoid writing the Kraken2 output file entirely:

```bash
kraken2 \
  --db /db \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --output - \
  --report /dev/null \
| barena \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -t 9606 \
  -c human_#.fastq.gz
```

{: .tip }
Streaming is useful when disk space is limited. The Kraken2 raw output for a 50 M read paired dataset is ~9 GB; streaming avoids writing it at all.

---

## Keep unclassified reads alongside matches

Useful when you want everything that is **not** clearly non-host:

```bash
barena \
  -1 R1.fastq.gz \
  -k run.kraken \
  -t 9606 \
  --include-unclassified \
  -c host_and_unknown.fastq.gz \
  -u confirmed_non_host.fastq.gz
```

---

## Multiple taxids

Filter reads from several species in a single run:

```bash
barena \
  -1 sample.fastq.gz \
  -k sample.kraken \
  -t 1314 -t 1313 -t 1311 \
  -c streptococcus.fastq.gz
```
