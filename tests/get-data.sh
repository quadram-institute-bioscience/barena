#!/usr/bin/env bash
OUTDIR=dataset/
mkdir -p $OUTDIR/
if [[ ! -e $OUTDIR/Zymo-nanopore.fastq.gz ]]; then
  curl -L https://s3.climb.ac.uk/nanopore/Zymo-GridION-EVEN-BB-SN-PCR-R10HC-flipflop.fq.gz -o $OUTDIR/Zymo-nanopore.fastq.gz
fi

if [[ ! -e $OUTDIR/SRR19995508_R1.fastq.gz ]]; then
  curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR199/008/SRR19995508/SRR19995508_1.fastq.gz -o $OUTDIR/SRR19995508_R1.fastq.gz
fi

if [[ ! -e $OUTDIR/SRR19995508_R2.fastq.gz ]]; then
  curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR199/008/SRR19995508/SRR19995508_2.fastq.gz -o $OUTDIR/SRR19995508_R2.fastq.gz
fi

if [[ ! -d $OUTDIR/db ]]; then
  curl -L https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08_GB_20260226.tar.gz -o $OUTDIR/db.tar.gz
  cd $OUTDIR && mkdir -p tmp && cd tmp && tar xfz ../db.tar.gz && cd .. && mv tmp db && rm db.tar.gz && cd -
fi

