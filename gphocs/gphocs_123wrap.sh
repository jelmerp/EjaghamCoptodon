#!/bin/bash -l

FILE_ID=$1
SCAFFOLD=$2 # SCAFFOLD=NC_022219.1
LOCUS_SIZE=$3
MIN_DISTANCE=$4
MAX_MISSING=$5
MASK_CPG=$6
MAX_SNPS=$7
CONVERT_FASTA=$8

if [ -z $CONVERT_FASTA ]; then CONVERT_FASTA=TRUE; fi

INDIR=seqdata/vcf_split

date
echo "Script: gphocs_submaster.sh"
echo "File ID: $FILE_ID"
echo "Scaffold: $SCAFFOLD"
echo "Locus size: $LOCUS_SIZE"
echo "Minimum interlocus distance: $MIN_DISTANCE"
echo "Maximum percentage of missing data: $MAX_MISSING"
echo "Maximum number of SNPS: $MAX_SNPS"
echo "Convert fasta: $CONVERT_FASTA"

if [ $CONVERT_FASTA == TRUE ]; then scripts/conversion/vcf2fasta.sh $FILE_ID $INDIR $SCAFFOLD; fi

scripts/gphocs/gphocs_2_chooseLoci.sh $SCAFFOLD # Only needs to be done once per scaffold, not for every FILE_ID set

scripts/gphocs/gphocs_3_createSeqs.sh $FILE_ID $SCAFFOLD $LOCUS_SIZE $MIN_DISTANCE $MAX_MISSING $MASK_CPG $MAX_SNPS

echo "Done with script."
date

