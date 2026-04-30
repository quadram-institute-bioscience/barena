export BINARY=./bin/barena
if [[ -d /qib/platforms/Informatics/transfer/outgoing/qib_databases/kraken2/k2_16S_Silva138_20200326 ]]; then
  export KRAKEN2_DEFAULT_DB=/qib/platforms/Informatics/transfer/outgoing/qib_databases/kraken2/k2_16S_Silva138_20200326
  echo "[DB] QIB location: $KRAKEN2_DEFAULT_DB"
elif [[ -d $USER/git/amplikraken/local/db/silva/16S_SILVA138_k2db ]]; then
  export KRAKEN2_DEFAULT_DB=$USER/git/amplikraken/local/db/silva/16S_SILVA138_k2db
  echo "[DB] Local user: $KRAKEN2_DEFAULT_DB"
else
  echo Missing kraken2 db in default locations
  exit 1
fi
export FASTQ=tests/data/dataset.fastq.gz
export KROUT=tests/data/dataset.raw
export BENCHDIR=./tests/benchmark/$(uname)

export OUTDIR=barena-out-3-child/
mkdir -p $BENCHDIR
hyperfine --warmup 1 --max-runs 7 --export-csv $BENCHDIR/benchmark-children.csv --export-markdown $BENCHDIR/benchmark-children.md \
 -n "barena (full output)" -n "barena (classified)" -n "kraken_tools" -n "kractor" -n "kraken + barena" \
 --prepare "grep -v '^#' _data/inspect.txt > _data/inspect.txt.skip && rm -f $OUTDIR/kractor.fq" \
 "$BINARY -1 $FASTQ -o $OUTDIR -k $KROUT  -d _data/inspect.txt  --taxon 3 --child" \
 "$BINARY -1 $FASTQ -c $OUTDIR/barena_class.fq -k $KROUT  -d _data/inspect.txt  --taxon 3 --child" \
 "extract_kraken_reads.py -k $KROUT -t 3 -o $OUTDIR/kraken.fq --fastq-output -s $FASTQ  --include-children -r _data/inspect.txt" \
 "kractor -i $FASTQ -o $OUTDIR/kractor.fq -k $KROUT -t 3 --children -r _data/inspect.txt.skip" \
 "kraken2 --threads 2 $FASTQ | $BINARY  -1 $FASTQ -o $OUTDIR/stream/ -k $KROUT  -d _data/inspect.txt  --taxon 3 --child"
kractor -i $FASTQ -o $OUTDIR/kractor.fq -k $KROUT -t 3 --children -r _data/inspect.txt.skip

export OUTDIR=barena-out-3-nochild/
hyperfine --warmup 1 --max-runs 7 --export-csv $BENCHDIR/benchmark-no-children.csv --export-markdown $BENCHDIR/benchmark-no-children.md \
  -n "barena (full output)" -n "barena (classified)" -n "kraken_tools" -n "kractor" -n "kraken + barena" \
  --prepare "rm -f $OUTDIR/kractor.fq" \
  "$BINARY -1 $FASTQ -o $OUTDIR/ -k $KROUT  -d _data/inspect.txt  --taxon 3" \
  "$BINARY -1 $FASTQ -c $OUTDIR/barena_class.fq -k $KROUT  -d _data/inspect.txt  --taxon 3" \
  "extract_kraken_reads.py -k $KROUT -t 3 -o $OUTDIR/kraken.fq --fastq-output -s $FASTQ" \
  "kractor -i $FASTQ -o $OUTDIR/kractor.fq -k $KROUT -t 3" \
  "kraken2 --threads 2 $FASTQ | $BINARY -1 $FASTQ -o $OUTDIR/stream  -d _data/inspect.txt --taxon 3"

kractor -i $FASTQ -o $OUTDIR/kractor.fq -k $KROUT -t 3
