#!/bin/bash

## Command-line args:
SETNAME=$1
REF=$2
GVCF_FOLDER=$3
VCF_INTERMEDIATE_FOLDER=$4
shift

count=0
while [ "$*" != "" ]
  do INDS[$count]=$1
  shift
  count=`expr $count + 1`
done

## Other variables:
INDS_COMMAND=$(for ind in ${INDS[@]}; do printf " --variant $GVCF_FOLDER/${ind}.raw_variants.g.vcf"; done)

## Software:
GATK=/proj/cmarlab/users/jelmer/software/GenomeAnalysisTK.jar

### Report settings/parameters:
date
echo "Script: SNPcalling4_jointGenotyping_allSites.sh"
echo "Set name: $SETNAME"
echo "Number of individuals: ${#INDS[@]}"
echo "Individuals: ${INDS[@]}"

### Run genotyper:
echo "Starting joint genotyping..."
java -Xmx10g -jar $GATK -T GenotypeGVCFs -R $REF $INDS_COMMAND -o $VCF_INTERMEDIATE_FOLDER/$SETNAME.output.raw.vcf

echo "Done with script."
date