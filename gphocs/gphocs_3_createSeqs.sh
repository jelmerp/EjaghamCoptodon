#!/bin/bash -l

FILE_ID=$1
SCAFFOLD=$2 # SCAFFOLD=NC_022219.1
LOCUS_SIZE=$3
MIN_DISTANCE=$4
MAX_MISSING=$5
MASK_CPG=$6
MAX_SNPS=$7

date
echo "Script: gphocs_createSeqs.sh"
echo "File ID: $FILE_ID"
echo "Scaffold: $SCAFFOLD"
echo "Locus size: $LOCUS_SIZE"

Rscript scripts/gphocs/gphocs_3_createSeqs.R $FILE_ID $SCAFFOLD $PWD $LOCUS_SIZE $MIN_DISTANCE $MAX_MISSING $MASK_CPG $MAX_SNPS

echo "Done with script."
date