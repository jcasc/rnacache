#!/bin/bash

module load java/13.0

LABPATH=/home/jcascitt/si
GTFPATH=$LABPATH/gencode.v38.annotation.gtf
REFPATH=$LABPATH/gencode.v38.transcripts.fa
FLUXPATH=/home/jcascitt/opt/bin/flux-simulator
export FLUX_MEM="32G"

num=0
while [ -d data_${1}-${num} ]; do
    num=$(( $num + 1 ))
done
mkdir data_${1}-${num}
cd data_${1}-${num}

awk '$0~"transcript_id"' $GTFPATH > annotation_tid_only.gtf
# filter gtf by reference transcriptome (e.g. for protein-coding):
# $LABPATH/filter_gtf.py $REFPATH annotation_tid_only.gtf > annotation_sim.gtf
# else:
mv annotation_tid_only.gtf annotation_sim.gtf

cp ../flux_params ./flux
$FLUXPATH -x -l -s -p flux
rm flux.bed
$LABPATH/split_fastq ./flux.fastq
rm flux.fastq

num_reads=$(($(wc -l reads1.fq | awk '{print $1}')/4))
echo $num_reads >> info

