#!/bin/bash -l

SOURCEDIR=$1
SOURCE_ID=$2
TARGETDIR=$3
TARGET_ID=$4

DP=$5
GQ=$6
MAXMISS=$7
MAF=$8

INPUT=$SOURCEDIR/$SOURCE_ID.vcf.gz
OUTPUT=$TARGETDIR/$TARGET_ID.DP$DP.GQ$GQ.MAXMISS$MAXMISS.MAF$MAF.vcf

date
echo "Script: filterVCF.sh"
echo "Source dir: $SOURCEDIR"
echo "Source ID: $SOURCE_ID"
echo "Target dir: $TARGETDIR"
echo "Target ID: $TARGET_ID"
echo "Min DP: $DP"
echo "Min GQ: $GQ"
echo "Max missing genotypes per site: $MAXMISS"
echo "Min MAF: $MAF"

vcftools --gzvcf $INPUT --out $OUTPUT --minDP $DP --minGQ $GQ --max-missing $MAXMISS --maf $MAF --max-non-ref-af 0.99 --min-alleles 2 --recode --recode-INFO-all
mv $OUTPUT.recode.vcf $OUTPUT
gzip -f $OUTPUT

echo "Done with script."
date
