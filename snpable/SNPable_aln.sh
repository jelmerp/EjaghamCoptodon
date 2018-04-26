#!/bin/bash

FOCALFILE=$1
SNPABLE=$2
REFERENCE=$3

FOCALFILE=$(basename $FOCALFILE)

date
echo "File: $FOCALFILE"

echo "Runnning SNPable alignment..."
bwa aln -R 1000000 -O 3 -E 3 $REFERENCE kmers/$FOCALFILE > kmers/$FOCALFILE.sai

echo "Converting .sai to .sam..."
bwa samse $REFERENCE kmers/$FOCALFILE.sai kmers/$FOCALFILE > kmer_alignments/$FOCALFILE.sam

echo "Done with script."
date