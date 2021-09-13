#!/usr/bin/env bash

set -x

NUM_THREADS=$1
ALIGNMENT_FILE=$2
REF_GENOME=$3
OUTPUT_DIR=$4
ALIGNMENT_FILENAME=`basename $ALIGNMENT_FILE`
FILENAME_ROOT=${ALIGNMENT_FILENAME::10}
BCF_FILE=$OUTPUT_DIR"/"$FILENAME_ROOT".GRCh38.p7.bcf"
VCF_FILE=$OUTPUT_DIR"/"$FILENAME_ROOT".GRCh38.p7.vcf"
MPILEUP_FILE=$OUTPUT_DIR"/"$FILENAME_ROOT".mpileup"
VCF_GZ_FILE=$VCF_FILE".gz"

echo "Using "$NUM_THREADS" threads"
echo "Will generate following output files "$BCF_FILE", "$VCF_FILE" and "$VCF_GZ_FILE

mkdir -p $OUTPUT_DIR
time bcftools mpileup --threads $NUM_THREADS -Ou -f $REF_GENOME $ALIGNMENT_FILE > $MPILEUP_FILE
time bcftools call --threads $NUM_THREADS -mv -Ob -o $BCF_FILE $MPILEUP_FILE
time bcftools view --threads $NUM_THREADS $BCF_FILE > $VCF_FILE

time bgzip -c $VCF_FILE > $VCF_GZ_FILE
time tabix -p vcf $VCF_GZ_FILE
