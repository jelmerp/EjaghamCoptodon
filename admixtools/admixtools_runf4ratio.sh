#!/bin/bash -l

## Command-line arguments:
SEQ_ID=$1
POPFILE_ID=$2

## Scripts:
ATOOLS_F4RATIO=/proj/cmarlab/users/jelmer/software/AdmixTools/bin/qpF4ratio

## Existing files:
PEDFILE=seqdata/plink/$SEQ_ID.ped 
MAPFILE=seqdata/plink/$SEQ_ID.map
INDFILE=analyses_input/admixtools/indfile_$POPFILE_ID.txt
POPFILE=analyses_input/admixtools/popfile_f4ratio_$POPFILE_ID.txt

## Files to create:
PARFILE_F4RATIO=analyses_input/admixtools/parfile_f4ratio_$POPFILE_ID.txt
F4R_OUTPUT=analyses_output/admixtools/raw/$SEQ_ID.f4ratio.out

## Report:
date
echo "Script: admixtools_runf4ratio.sh"
echo "Seq ID: $SEQ_ID"
echo "Popfile & Seqfile ID: $POPFILE_ID"
echo "Output folder: analyses/admixtools/output/raw/"

## Create parfile:
echo "Creating parfile $PARFILE_F4RATIO..."
printf "genotypename:\t$PEDFILE\n" > $PARFILE_F4RATIO
printf "snpname:\t$MAPFILE\n" >> $PARFILE_F4RATIO
printf "indivname:\t$INDFILE\n" >> $PARFILE_F4RATIO
printf "popfilename:\t$POPFILE\n" >> $PARFILE_F4RATIO
printf "printsd:\tYES\n" >> $PARFILE_F4RATIO
echo "Parfile:"
cat $PARFILE_F4RATIO

## Run f4-ratio test:
echo "Running qpf4ratio..."
$ATOOLS_F4RATIO -p $PARFILE_F4RATIO > $F4R_OUTPUT

echo "Done with script."
date