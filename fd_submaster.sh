#!/bin/bash -l

## Command line arguments:
POPFILE=$1
GENOFILE=$2
MINDEPTH=$3
MINQUAL=$4
WINSIZE=$5
STEPSIZE=$6
MINSITES=$7
MINDATA=$8
FIRSTLINE=$9
LASTLINE=${10}

## Other variables:
NCORES=4
OUTPUT_DIR=analyses/admixblocks/output
FD_SUBMISSION_SCRIPT=scripts/admixblocks/fd_run.sh

## Report:
date
echo "Script: abbababa_submaster.sh"
echo "POPFILE: $POPFILE"
echo "GENOFILE: $GENOFILE"
echo "MINDEPTH: $MINDEPTH"
echo "MINQUAL: $MINQUAL"
echo "WINSIZE: $WINSIZE"
echo "STEPSIZE: $STEPSIZE"
echo "MINSITES: $MINSITES"
echo "MINDATA: $MINDATA"

## Determine lines of "popfile" (population configurations) for which the script should be run:
if [ ! -z $FIRSTLINE ]
then
	SEQ=$(seq $FIRSTLINE $LASTLINE)
else
	FIRSTLINE=1
	LASTLINE=$(cat $POPFILE | wc -l)
fi
echo "First line: $FIRSTLINE"
echo "Last line: $LASTLINE"

## For each line of the popfile, run abbababa_run.sh:
for i in $(seq $FIRSTLINE $LASTLINE)
do
	echo "Line nr: $i"
	P1=$(head -n $i $POPFILE | tail -1 | cut -d " " -f 1)
	P2=$(head -n $i $POPFILE | tail -1 | cut -d " " -f 2)
	P3=$(head -n $i $POPFILE | tail -1 | cut -d " " -f 3)
	P4=$(head -n $i $POPFILE | tail -1 | cut -d " " -f 4)
	OUTPUT=$OUTPUT_DIR/ABBABABAoutput_EjaC.Dstat.DP$MINDEPTH.GQ$MINQUAL.win$WINSIZE.step$STEPSIZE.popfileline$i.csv
	bsub -n $NCORES -R "span[hosts=1]" -q day -o slurm.abbababa.DP$MINDEPTH.GQ$MINQUAL.win$WINSIZE.step$STEPSIZE.popfileline$i.txt $FD_SUBMISSION_SCRIPT $GENOFILE $OUTPUT $NCORES $WINSIZE $STEPSIZE $MINSITES $MINDATA $P1 $P2 $P3 $P4
done

echo "Done with script fd_submaster.sh"