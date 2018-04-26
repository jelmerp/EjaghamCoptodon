#!/bin/bash

IND=$1
FASTQ_ID1=$2
FASTQ_ID2=$3
REFERENCE=$4
FASTQ_FOLDER=$6
BAM_FOLDER=$7

echo "Aligning fastq files to reference genome for: $IND"
echo "Fastq ID 1: $FASTQ_ID1"
echo "Fastq ID 2: $FASTQ_ID2"

bwa mem $REFERENCE -aM -t 8 -R "@RG\tID:group1\tSM:$IND\tPL:illumina\tLB:lib1" $FASTQ_FOLDER/$FASTQ_ID1.fastq.gz $FASTQ_FOLDER/$FASTQ_ID2.fastq.gz > $BAM_FOLDER/$IND.sam

echo "Done with script."





