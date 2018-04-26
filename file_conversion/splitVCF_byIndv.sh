#!/bin/bash -l

## INDS should be a comma-separated list of individuals IDs that should be included. To exclude an individual, place a ^ before the ID.

SOURCEDIR=$1
SOURCE_ID=$2
INDS=$3 # INDS=Cdec088,Cdec328,Ceja262,Ceja408,Cfus085,Cfus350,Cfus503,Ckot383,Ckot499,TguiNG2,TguiNG5,TguiMA1,TguiMA2,TguiMA4
TARGET_ID=$4
TARGETDIR=$5
MAF=$6

VCF_SOURCE=$SOURCEDIR/$SOURCE_ID.vcf.gz

date
echo "Script: splitVCF_byIndv.sh"
echo "VCF source: $VCF_SOURCE"
echo "Minor allele frequency: $MAF"
echo "Individuals: $INDS"


if [ -z "$MAF" ]
then
	VCF_TARGET=/proj/cmarlab/users/jelmer/cichlids/seqdata/vcf_split/$TARGET_ID.vcf
	echo "VCF target: $VCF_TARGET.gz"
	echo "Extracting individuals..."
	bcftools view -O z -s $INDS $VCF_SOURCE > $VCF_TARGET.gz
else
	VCF_INTERMED=/proj/cmarlab/users/jelmer/cichlids/seqdata/vcf_split/tmp.$TARGET_ID.vcf.gz
	VCF_TARGET=/proj/cmarlab/users/jelmer/cichlids/seqdata/vcf_split/$TARGET_ID.MAF$MAF.vcf
	echo "VCF target: $VCF_TARGET.gz"
	
	echo "Extracting individuals and filtering by MAF $MAF..."
	bcftools view -O z -s $INDS $VCF_SOURCE > $VCF_INTERMED
	vcftools --gzvcf $VCF_INTERMED --maf $MAF --max-non-ref-af 0.99 --min-alleles 2 --recode --recode-INFO-all --stdout > $VCF_TARGET
	gzip -f $VCF_TARGET
	rm $VCF_INTERMED
fi

if [ -z "$TARGETDIR" ]
	then 
		echo "No target folder set, file will be kept in $VCF_TARGET.gz"
	else
		echo "Target folder: '$TARGETDIR'"
		mv $VCF_TARGET.gz $TARGETDIR/$TARGET_ID.MAF$MAF.vcf.gz
fi

echo "Done with script."
date
