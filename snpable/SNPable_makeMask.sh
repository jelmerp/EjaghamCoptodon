#!/bin/bash

date
echo "Script: SNPable_makeMask.sh..."

SNPABLE=$1
k_value=$2
r_value=$3

printf "\n\n"
echo "Merging sam files, sorting resulting file..."
samtools merge kmer_alignments/merged.sam `ls -d $PWD/kmer_alignments/*.sam`
samtools sort -n -O sam -o kmer_alignments/merged.sorted.sam kmer_alignments/merged.sam

printf "\n\n"
echo "Generating raw mask..."
$SNPABLE/gen_raw_mask.pl kmer_alignments/merged.sorted.sam > mask_fasta/rawMask_l50_r50.fa

printf "\n\n"
echo "Generating final mask..."
$SNPABLE/gen_mask -l $k_value -r $r_value mask_fasta/rawMask_l50_r50.fa > mask_fasta/mask_l50_r50.fa

printf "\n\n"
echo "Converting to bed file..."
module load python/3.3.3
$SCR/SNPable_makeMappabilityMask.py

echo "Done with script."
date