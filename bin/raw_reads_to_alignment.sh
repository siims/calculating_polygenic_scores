#!/bin/bash

set -x

NUM_THREADS=$1
BOWTIE_INDEX=$2
FASTQ1=$3
FASTQ2=$4
OUTPUT_DIR=$5
FASTQ_FILENAME=`basename $FASTQ1`
FILENAME_ROOT=${FASTQ_FILENAME::10}

echo "will start process in $NUM_THREADS threads"
echo "using bowtie2 index: $BOWTIE_INDEX"
echo "using fastq1: $FASTQ1"
echo "using fastq2: $FASTQ2"
UNSORTED_BAM="$OUTPUT_DIR/$FILENAME_ROOT""_hg38.p7.unsorted.bam"
SORTED_BAM="$OUTPUT_DIR/$FILENAME_ROOT""_hg38.p7.bam"
echo "raw results will be written into: $OUTPUT_DIR as $UNSORTED_BAM,$SORTED_BAM,$SORTED_BAM.bai"
mkdir -p $OUTPUT_DIR

time bowtie2 -p $NUM_THREADS -x $BOWTIE_INDEX -1 $FASTQ1 -2 $FASTQ2 | samtools view -bS - > $UNSORTED_BAM # -p for number of parallel threads
time samtools sort -@ $NUM_THREADS $UNSORTED_BAM -o $SORTED_BAM # decreases the size of the alignment file to around 40GB; -@ for number of parallel threads
time samtools index -@ $NUM_THREADS $SORTED_BAM # created GENOME12345_hg38.p7.bam.bai index file for fast lookups; -@ for number of parallel threads
