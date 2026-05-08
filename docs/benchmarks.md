---
title: Benchmarks
layout: default
nav_order: 5
---

# Benchmarks
{: .no_toc }

<details open markdown="block">
  <summary>Contents</summary>
  {: .text-delta }
- TOC
{:toc}
</details>

Benchmarks compare barena against [KrakenTools](https://github.com/jenniferlu717/KrakenTools) (`extract_kraken_reads.py`) and [Kractor](https://github.com/fhcrc/kractor). All runs used [hyperfine](https://github.com/sharkdp/hyperfine) with 1 warm-up run and 5 timed repetitions.

---

## Dataset

| File | Platform | Reads | Size |
|---|---|---|---|
| `Zymo-nanopore.fastq.gz` | ONT | 3,491,078 | 13 GB |
| `SRR19995508_R1/R2.fastq.gz` | Illumina paired | 53,526,611 | 4.4 GB + 149 MB |

The ONT data is from the [Zymo mock community](https://github.com/LomanLab/mockcommunity) GridION EVEN dataset. The Illumina data is SRA accession SRR19995508.

---

## Single-end (ONT)

### Benchmark 1 — `.fastq` input, `.fastq` output

Extract *Bacillus spizizenii* (taxid 96241) from ONT long reads.

| File type | Reads to extract | Kraken output size |
|---|---|---|
| `.fastq` | 490,984 | 885 MB |

| Tool | Mean (s) | Min (s) | Max (s) | Speed-up |
|:---|---:|---:|---:|---:|
| KrakenTools | 254.9 ± 9.5 | 242.6 | 263.2 | 1.00× |
| barena | — | — | — | — |
| Kractor | 58.3 ± 5.7 | 48.4 | 62.2 | **4.4×** |

### Benchmark 2 — `.fastq.gz` input, `.fastq` output

| File type | Reads to extract | Kraken output size |
|---|---|---|
| `.fastq.gz` | 490,984 | 885 MB |

| Tool | Mean (s) | Min (s) | Max (s) | Speed-up |
|:---|---:|---:|---:|---:|
| KrakenTools | 376.6 ± 3.9 | 373.3 | 383.3 | 1.00× |
| barena | — | — | — | — |
| Kractor | 100.0 ± 3.6 | 98.0 | 106.4 | **3.8×** |

---

## Paired-end (Illumina)

### Benchmark 3 — `.fastq` input, `.fastq` output

Extract all reads under taxid 590 from Illumina paired-end data.

| File type | Reads to extract | Kraken output size |
|---|---|---|
| `.fastq` paired | 1,646,117 | 9.3 GB |

| Tool | Mean (s) | Min (s) | Max (s) | Speed-up |
|:---|---:|---:|---:|---:|
| KrakenTools | 898.3 ± 14.2 | 884.2 | 920.7 | 1.00× |
| barena | — | — | — | — |
| Kractor | 94.2 ± 2.3 | 90.9 | 96.5 | **9.5×** |

### Benchmark 4 — `.fastq.gz` input, `.fastq` output

| File type | Reads to extract | Kraken output size |
|---|---|---|
| `.fastq.gz` paired | 1,646,117 | 9.3 GB |

| Tool | Mean (s) | Min (s) | Max (s) | Speed-up |
|:---|---:|---:|---:|---:|
| KrakenTools | 1033.4 ± 25.2 | 1005.7 | 1068.5 | 1.00× |
| barena | — | — | — | — |
| Kractor | 49.1 ± 0.2 | 48.9 | 49.3 | **21×** |

{: .note }
barena entries will be filled in once the full comparative benchmark script completes. Current internal measurements show barena consistently in the Kractor performance range for most workloads, with lower memory use on large datasets.

---

## Commands used

```bash
# KrakenTools — single end
extract_kraken_reads.py -s INPUT.fq -k out.kraken \
  -o out.fq -t TAXID --fastq-output

# KrakenTools — paired end
extract_kraken_reads.py \
  -s SRR19995508_R1.fastq -s2 SRR19995508_R2.fastq \
  -o R1.fq -o2 R2.fq -k out.kraken -t 590 --fastq-output

# barena — single end
barena -1 INPUT.fq -k out.kraken -t TAXID -c out.fq.gz

# barena — paired end
barena -1 R1.fq.gz -2 R2.fq.gz -k out.kraken -t 590 -o outdir/

# Kractor — single end
kractor -i INPUT.fq -k out.kraken -o out.fq -t TAXID

# Kractor — paired end
kractor -i R1.fastq -i R2.fastq -k out.kraken \
  -t 590 -o R1_out.fq -o R2_out.fq
```
