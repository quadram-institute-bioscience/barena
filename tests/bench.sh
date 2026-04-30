export BINARY=./barena_linux
export KRAKEN2_DEFAULT_DB=/Users/telatina/git/amplikraken/local/db/silva/16S_SILVA138_k2db
export FASTQ=tests/data/dataset.fastq.gz
export KROUT=tests/data/dataset.raw


export OUTDIR=barena-out-3-child/

hyperfine --warmup 1 --max-runs 7 --export-csv tests/benchmark-children.csv --export-markdown tests/benchmark-children.md \
 --prepare "grep -v '^#' _data/inspect.txt > _data/inspect.txt.skip && rm -f $OUTDIR/kractor.fq" \
 "$BINARY -1 $FASTQ -o $OUTDIR -k $KROUT  -d _data/inspect.txt  --taxon 3 --child" \
 "$BINARY -1 $FASTQ -c $OUTDIR/barena_class.fq -k $KROUT  -d _data/inspect.txt  --taxon 3 --child" \
 "extract_kraken_reads.py -k $KROUT -t 3 -o $OUTDIR/kraken.fq --fastq-output -s $FASTQ  --include-children -r _data/inspect.txt" \
 "kractor -i $FASTQ -o $OUTDIR/kractor.fq -k $KROUT -t 3 --children -r _data/inspect.txt.skip" \
 "kraken2 --threads 2 $FASTQ | $BINARY  -1 $FASTQ -o $OUTDIR/stream/ -k $KROUT  -d _data/inspect.txt  --taxon 3 --child"
kractor -i $FASTQ -o $OUTDIR/kractor.fq -k $KROUT -t 3 --children -r _data/inspect.txt.skip

export OUTDIR=barena-out-3-nochild/
hyperfine --warmup 1 --max-runs 7 --export-csv tests/benchmark-no-children.csv --export-markdown tests/benchmark-no-children.md \
  --prepare "rm -f $OUTDIR/kractor.fq" \
  "$BINARY -1 $FASTQ -o $OUTDIR/ -k $KROUT  -d _data/inspect.txt  --taxon 3" \
  "$BINARY -1 $FASTQ -c $OUTDIR/barena_class.fq -k $KROUT  -d _data/inspect.txt  --taxon 3" \
  "extract_kraken_reads.py -k $KROUT -t 3 -o $OUTDIR/kraken.fq --fastq-output -s $FASTQ" \
  "kractor -i $FASTQ -o $OUTDIR/kractor.fq -k $KROUT -t 3" \
  "kraken2 --threads 2 $FASTQ | $BINARY -1 $FASTQ -o $OUTDIR/stream  -d _data/inspect.txt --taxon 3"

kractor -i $FASTQ -o $OUTDIR/kractor.fq -k $KROUT -t 3
