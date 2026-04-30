| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `./barena_linux -1 tests/data/dataset.fastq.gz -o barena-out-3-child/ -k tests/data/dataset.raw  -d _data/inspect.txt  --taxon 3 --child` | 48.6 ± 1.5 | 46.5 | 49.8 | 1.27 ± 0.06 |
| `./barena_linux -1 tests/data/dataset.fastq.gz -c barena-out-3-child//barena_class.fq -k tests/data/dataset.raw  -d _data/inspect.txt  --taxon 3 --child` | 38.2 ± 1.2 | 36.8 | 40.5 | 1.00 |
| `extract_kraken_reads.py -k tests/data/dataset.raw -t 3 -o barena-out-3-child//kraken.fq --fastq-output -s tests/data/dataset.fastq.gz  --include-children -r _data/inspect.txt` | 609.3 ± 24.4 | 585.2 | 655.2 | 15.95 ± 0.82 |
| `kractor -i tests/data/dataset.fastq.gz -o barena-out-3-child//kractor.fq -k tests/data/dataset.raw -t 3 --children -r _data/inspect.txt.skip` | 55.6 ± 2.0 | 52.7 | 57.6 | 1.46 ± 0.07 |
| `kraken2 --threads 2 tests/data/dataset.fastq.gz \| ./barena_linux  -1 tests/data/dataset.fastq.gz -o barena-out-3-child//stream/ -k tests/data/dataset.raw  -d _data/inspect.txt  --taxon 3 --child` | 250.8 ± 5.7 | 245.2 | 258.9 | 6.57 ± 0.26 |
