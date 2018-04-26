#!/bin/bash

## Command-line args:
IND=$1
REF=$2
BAM_PROCESSED_FOLDER=$3
GVCF_FOLDER=$4

## Other variables:
INPUT=$BAM_PROCESSED_FOLDER/$IND.dedup.bam
OUTPUT=$GVCF_FOLDER/$IND.raw_variants.g.vcf

## Software:
GATK=/proj/cmarlab/users/jelmer/software/GenomeAnalysisTK.jar

## Report:
echo "Variant discovery with GATK for: $IND"
echo "Reference sequence: $REF"
echo "Input file name: $INPUT"
echo "Output file name: $OUTPUT"

## Run GATK:
echo "Starting variant discovery..."

java -Xmx14g -jar $GATK -T HaplotypeCaller -R $REF -I $INPUT --genotyping_mode DISCOVERY --emitRefConfidence GVCF -mmq 5 -nct 4 -o $OUTPUT

echo "Done with script."