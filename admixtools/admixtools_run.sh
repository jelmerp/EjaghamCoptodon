#!/bin/bash -l

## Command-line arguments:
SEQ_ID=$1
POPFILE_ID=$2
INDFILE_ID=$3
LINE=$4

## Scripts:
ATOOLS_DSTAT=/proj/cmarlab/users/jelmer/software/AdmixTools/bin/qpDstat

## Existing files:
PEDFILE=seqdata/plink/$SEQ_ID.ped 
MAPFILE=seqdata/plink/$SEQ_ID.map
INDFILE=analyses_input/admixtools/indfile_$INDFILE_ID.txt
POPFILE=analyses_input/admixtools/popfile_$POPFILE_ID.txt

## Files to create:
PARFILE_DMODE=analyses_input/admixtools/parfile_dmode_$SEQ_ID.$POPFILE_ID.$LINE.txt
PARFILE_FMODE=analyses_input/admixtools/parfile_fmode_$SEQ_ID.$POPFILE_ID.$LINE.txt
OUTPUT_DIR=analyses_output/admixtools/raw/

## Report:
date
echo "Script: admixtools_run.sh"
echo "Seq ID: $SEQ_ID"
echo "Popfile & Indfile ID: $POPFILE_ID"
echo "Line: $LINE"
echo "Parfile (dmode): $PARFILE_DMODE"
echo "Parfile (fmode): $PARFILE_FMODE"
echo "Output folder: $OUTPUT_DIR"

## Create eigenstat parfile:
# parfile for d-mode:
printf "\n \n \n \n"
echo "Creating EIGENSTAT parfile..."
printf "genotypename:\t$PEDFILE\n" > $PARFILE_DMODE
printf "snpname:\t$MAPFILE\n" >> $PARFILE_DMODE
printf "indivname:\t$INDFILE\n" >> $PARFILE_DMODE
printf "popfilename:\t$POPFILE\n" >> $PARFILE_DMODE
printf "printsd:\tYES\n" >> $PARFILE_DMODE
echo "Parfile:"
cat $PARFILE_DMODE

# parfile for f-mode:
cp $PARFILE_DMODE $PARFILE_FMODE
printf "f4mode:\tYES\n" >> $PARFILE_FMODE
printf "\n \n \n \n"

## Run admixtools:
if [ $LINE == ALL ]
then
	echo "Running all popfile lines at once..."
	echo "Running qpDstat in D-mode..."
	$ATOOLS_DSTAT -p $PARFILE_DMODE > $OUTPUT_DIR/$SEQ_ID.dmode.out
	echo "Running qpDstat in F-mode..."
	$ATOOLS_DSTAT -p $PARFILE_FMODE > $OUTPUT_DIR/$SEQ_ID.fmode.out
else
	echo "Running one popfile line: $LINE"
	echo "Running qpDstat in D-mode..."
	$ATOOLS_DSTAT  -l $LINE -h $LINE -p $PARFILE_DMODE > $OUTPUT_DIR/$SEQ_ID.$POPFILE_ID.line${LINE}.dmode.out
	echo "Running qpDstat in F-mode..."
	$ATOOLS_DSTAT -l $LINE -h $LINE -p $PARFILE_FMODE > $OUTPUT_DIR/$SEQ_ID.$POPFILE_ID.line${LINE}.fmode.out
fi

echo "Done with script."
date