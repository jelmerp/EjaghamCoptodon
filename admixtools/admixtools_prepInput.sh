#!/bin/bash -l

## Needs loaded: bcftools

## Command-line arguments:
FILE_ID=$1
VCF=$2
INDFILE=$3
MAF=$4

## Scripts:
VCF2PLINK_SCRIPT=scripts/conversion/vcf2plink.sh
SPLITVCF_SCRIPT=scripts/conversion/splitVCF_byIndv.sh

## Report:
date
echo "Script: admixtools_prepInput.sh"
echo "FILE_ID: $FILE_ID"
if [ ! -z $INDS ]; then echo "INDS: $INDS"; echo "MAF: $MAF"; fi
echo "VCF: $VCF"

## Convert VCF to plink format:
LD_MAX=1
$VCF2PLINK_SCRIPT $FILE_ID $MAF $LD_MAX

## Create eigenstat indfile:
bcftools query -l $VCF | sed ':a;N;$!ba;s/\n/ U U\n/g' | sed 's/\([0-9]\)$/\1 U U/g' | sed 's/NOL$/NOL U U/g' > $INDFILE

echo "Done with script."
date