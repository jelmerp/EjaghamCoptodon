#!/bin/bash -l

module load bedtools

FILE_ID=$1
SCAFFOLD=$2
START=$3 
END=$4
INDIR=$5
OUTDIR=$6

VCF_IN=$INDIR/$FILE_ID.vcf.gz
VCF_OUT=$OUTDIR/$FILE_ID.$SCAFFOLD.$START.$END.vcf

date
echo "Script: splitVCF_byCoords.sh"
echo "File ID: $FILE_ID"
echo "Scaffold: $SCAFFOLD"
echo "Start: $START"
echo "End: $END"
echo "Indir: $INDIR"
echo "Outdir: $OUTDIR"
echo "Output file: $VCF_OUT"

if [ ! -d tmpfiles ]; then mkdir tmpfiles; fi

printf "${SCAFFOLD}\t${START}\t${END}" > tmpfiles/$SCAFFOLD.$START.$END.bed

bedtools intersect -header -a $VCF_IN -b tmpfiles/$SCAFFOLD.$START.$END.bed > $VCF_OUT

printf "\n"
echo "Done with script."
date
