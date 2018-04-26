################################################################################################################
#### Set-up ####
INDS=(Cdec088 Cdec328 Ceja262 Ceja408 Cfus085 Cfus350 Cfus503 Ckot383 Ckot499 SgalMA1 TguiMA1 TguiMA2 TguiMA4 TguiNG2 TguiNG5) # Individuals
SETNAME=EjaC.Ckot.Dstat

REF=seqdata/reference/Oreochromis.fna # Reference genome

FASTQ_FOLDER=seqdata/fastq_raw
BAM_ALIGNED_FOLDER=seqdata/bam_aligned
BAM_PROCESSED_FOLDER=seqdata/bam_processed
GVCF_FOLDER=seqdata/vcf_gVCF
VCF_INTERMEDIATE_FOLDER=seqdata/vcf_intermediate
VCF_FINAL_FOLDER=seqdata/vcf_master


#### Step 1: Aligning raw sequences with bwa. ####
FASTQ1_IDS=(CAM06_S6_L001_R1_001 CAM02_S2_L001_R1_001 CAM08_S8_L001_R1_001 CAM07_S7_L001_R1_001 CAM03_S3_L001_R1_001 \
	CAM05_S5_L001_R1_001 CAM04_S4_L001_R1_001 CAM10_S10_L001_R1_001 CAM09_S9_L001_R1_001 CAM43_S256_L006_R1_001 CAM37_S243_L005_R1_001 \
	CAM44_S257_L006_R1_001 CAM48_S245_L005_R1_001 CAM28_S23_L007_R1_001 CAM22_S18_L006_R1_001) 
FASTQ2_IDS=(CAM06_S6_L001_R2_001 CAM02_S2_L001_R2_001 CAM08_S8_L001_R2_001 CAM07_S7_L001_R2_001 CAM03_S3_L001_R2_001 \
	CAM05_S5_L001_R2_001 CAM04_S4_L001_R2_001 CAM10_S10_L001_R2_001 CAM09_S9_L001_R2_001 CAM43_S256_L006_R2_001 CAM37_S243_L005_R2_001 \
	CAM44_S257_L006_R2_001 CAM48_S245_L005_R2_001 CAM28_S23_L007_R2_001 CAM22_S18_L006_R2_001) 
MAX=`expr ${#INDS[@]} - 1` 
for NR in $(seq 0 $MAX)
do
	IND=${INDS[$NR]}
	FASTQ1_ID=${FASTQ1_IDS[$NR]}
	FASTQ2_ID=${FASTQ2_IDS[$NR]}
	echo "Nr: $NR IND: $IND FASTQ1_ID: $FASTQ1_ID FASTQ2_ID: $FASTQ2_ID"
	bsub -n 8 -R "span[hosts=1]" -M 16 -N -o slurm.SNPcalling1.$IND.txt -q week \
		scripts/snp_calling/SNPcalling1_alignment.sh $IND $FASTQ_ID1 $FASTQ_ID2 $REF $FASTQ_FOLDER $BAM_ALIGNED_FOLDER
done


#### Step 2: Processing (sorting, deduplication, indexing) of sam files with Picard. ####
for IND in ${INDS[@]}
do
	bsub -M 16 -N -o slurm.SNPcalling2.$IND.txt -q day \
		scripts/snp_calling/SNPcalling2_sort_dedup_index.sh $IND $BAM_ALIGNED_FOLDER $BAM_PROCESSED_FOLDER
done


#### Step 3: Variant discovery with GATK ####
for IND in ${INDS[@]}
do
	bsub -M 16 -N -o slurm.SNPcalling3.$IND.txt -n 4 -q whole_node \
		scripts/snp_calling/SNPcalling3_variantDiscovery.sh $IND $REF $BAM_PROCESSED_FOLDER $GVCF_FOLDER
done


#### Step 4: Joint genotyping of GVCF files for multiple individuals ####
bsub -M 12 -N -o slurm.SNPcalling4.$SETNAME.txt -q week \
	scripts/snp_calling/SNPcalling4_jointGenotyping.sh $SETNAME REF $GVCF_FOLDER $VCF_INTERMEDIATE_FOLDER ${INDS[@]}


#### Step 5: Filtering of VCF files ####
bsub -M 12 -N -o slurm.SNPcalling5.$SETNAME.txt -q day \
	scripts/snp_calling/SNPcalling5_filter.sh $SETNAME $REF $VCF_INTERMEDIATE_FOLDER $VCF_FINAL_FOLDER


################################################################################################################
#### Create dictionary for and index genome ####
#java -jar picard.jar CreateSequenceDictionary R=Oreochromis.fna O=Oreochromis.dict
#samtools faidx Oreochromis.fna
#bwa index Oreochromis.fna