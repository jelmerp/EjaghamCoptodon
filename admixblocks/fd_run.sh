#!/bin/bash -l

module load python/2.7.1 # Scripts only seem to work with this version of Python

## Command line arguments:
GENO=$1
OUTPUT=$2
NCORES=$3
WINSIZE=$4
STEPSIZE=$5
MINSITES=$6
MINDATA=$7
P1=$8
P2=$9
P3=${10}
P4=${11}

## Scripts:
FD_SCRIPT=scripts/misc/ABBABABAwindows.py

## Report:
date
echo "Script: abbababa_run.sh"
echo "Geno file: $GENO"
echo "Nr of threads: $NCORES"
echo "Window size: $WINSIZE"
echo "Step size: $STEPSIZE"
echo "Minimum nr of sites: $MINSITES"
echo "Minimum prop of good data: $MINDATA"
echo "Individuals in population 1: $P1"
echo "Individuals in population 2: $P2"
echo "Individuals in population 3: $P3"
echo "Individuals in population 4 (outgroup): $P4"

## Run Fd script:
python $FD_SCRIPT -g $GENO -o $OUTPUT -f phased --windSize $WINSIZE --stepSize $STEPSIZE --minSites $MINSITES \
	--minData $MINDATA --Threads $NCORES -P1 pop1 $P1 -P2 pop2 $P2 -P3 pop3 $P3 -O pop4 $P4

echo "Done with script."
date