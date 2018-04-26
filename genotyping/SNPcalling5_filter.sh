#!/bin/bash

## Command-line args:
FILE_ID=$1
REF=$2
VCF_INTERMEDIATE_FOLDER=$3
VCF_FINAL_FOLDER=$4

## Software:
GATK=/proj/cmarlab/users/jelmer/software/GenomeAnalysisTK.jar

## Other variables:
VCF_BASE=$VCF_INTERMEDIATE_FOLDER/$FILE_ID

### Report settings/parameters:
date
echo "Script: SNPcalling5_filter.sh"
echo "VCF input file: $VCF_BASE"

### Step A: select only SNPs:
echo "Starting Step A (select only SNPS)..." # 17 minutes
java -jar $GATK -T SelectVariants -R $REF -V $VCF_BASE.output.raw.vcf -selectType SNP -o $VCF_BASE.SNPs.raw.vcf 

### Step B: filter SNPs:
echo "Starting first part of Step B (filter SNPS with GATK)..." # 2.4h
java -jar $GATK -T VariantFiltration -R $REF -V $VCF_BASE.SNPs.raw.vcf --filterExpression "QD < 2.0" --filterName "QD_lt2" \
	--filterExpression "FS > 60.0" --filterName "FS_gt60" --filterExpression "SOR > 3" --filterName "SOR_gt3" \
	--filterExpression "MQ < 40.0" --filterName "MQ_lt40" --filterExpression "MQRankSum < -12.5" --filterName "MQRankSum_ltm12.5" \
	--filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum_ltm8" -o $VCF_BASE.SNPs.GATKfiltSoft.vcf

echo "Starting second part of Step B (filter SNPS to biallelic only with bcftools)..."

bcftools view -m2 -M2 -v snps $VCF_BASE.SNPs.GATKfilt.vcf -O z > $VCF_BASE.SNPs.GATKfiltSoft.biallelic.vcf.gz

vcftools --vcf $VCF_BASE.SNPs.GATKfiltSoft.biallelic.vcf.gz --remove-filtered-all --max-non-ref-af 0.99 --recode \
	--recode-INFO-all --stdout > $VCF_BASE.SNPs.GATKfilt.biallelic.vcf

### Step C: gzip vcf files:
echo "Starting Step E (gzip vcf files)..."
gzip -f $VCF_BASE.output.raw.vcf
gzip -f $VCF_BASE.SNPs.raw.vcf
gzip -f $VCF_BASE.SNPs.GATKfiltSoft.vcf
gzip -f $VCF_BASE.SNPs.GATKfiltSoft.biallelic.vcf
gzip -f $VCF_BASE.SNPs.GATKfilt.biallelic.vcf
mv $VCF_BASE.SNPs.GATKfilt.biallelic.vcf.gz $VCF_FINAL_FOLDER/$FILE_ID.SNPs.GATKfilt.biallelic.vcf.gz

echo "Done with script"
date