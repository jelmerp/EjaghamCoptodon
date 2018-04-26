#!/bin/bash -l

## Command-line arguments:
VCF_ID=$1
SCAFFOLD=$2
INDIR=$3
OUTDIR=$4

## Other variables:
VCF_IN=$INDIR/$VCF_ID.vcf.gz
VCF_OUT=$OUTDIR/$VCF_ID.$SCAFFOLD.vcf.gz

## Report:
date
echo "SCRIPT: splitVCF_byScaffold.sh"
echo "VCF ID: $VCF_ID"
echo "SCAFFOLD: $SCAFFOLD"
echo "VCF_IN: $VCF_IN"
echo "VCF_OUT: $VCF_OUT"

## Split VCF:
vcftools --gzvcf $VCF_IN --chr $SCAFFOLD --recode --recode-INFO-all --stdout | gzip -c > $VCF_OUT