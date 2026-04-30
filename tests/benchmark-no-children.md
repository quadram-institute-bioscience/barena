| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `./barena_linux -1 tests/data/dataset.fastq.gz -o barena-out-3-nochild// -k tests/data/dataset.raw  -d _data/inspect.txt  --taxon 3` | 40.4 ± 1.0 | 39.1 | 41.5 | 1.59 ± 0.06 |
| `./barena_linux -1 tests/data/dataset.fastq.gz -c barena-out-3-nochild//barena_class.fq -k tests/data/dataset.raw  -d _data/inspect.txt  --taxon 3` | 25.4 ± 0.7 | 24.6 | 26.5 | 1.00 |
| `extract_kraken_reads.py -k tests/data/dataset.raw -t 3 -o barena-out-3-nochild//kraken.fq --fastq-output -s tests/data/dataset.fastq.gz` | 558.8 ± 33.7 | 526.2 | 614.7 | 21.99 ± 1.46 |
| `kractor -i tests/data/dataset.fastq.gz -o barena-out-3-nochild//kractor.fq -k tests/data/dataset.raw -t 3` | 45.1 ± 1.9 | 41.4 | 46.9 | 1.77 ± 0.09 |
| `kraken2 --threads 2 tests/data/dataset.fastq.gz \| ./barena_linux -1 tests/data/dataset.fastq.gz -o barena-out-3-nochild//stream  -d _data/inspect.txt --taxon 3` | 263.1 ± 6.1 | 256.2 | 270.8 | 10.35 ± 0.37 |
