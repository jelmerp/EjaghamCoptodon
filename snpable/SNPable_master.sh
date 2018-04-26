## Set up directories etc:
BASEDIR=analyses_output/SNPable_mask # Working directory for SNPable mask creation
SNPABLE=software/SNPable # SNPable software folder
REFERENCE=reference/Oreochromis.fna # Reference genome fasta
SCR=scripts/SNPable # Dir should have these scripts: SNPable_extractKmers.sh, SNPable_aln.sh, SNPable_makeMask.sh, SNPable_makeMappabilityMask.py

mkdir $BASEDIR $BASEDIR/kmers $BASEDIR/kmer_alignments $BASEDIR/mask_bed $BASEDIR/mask_fasta
cd $BASEDIR


## Step 1: extract all overlapping k-mer subsequences as read sequences:
bsub -q hour -N -o slurm.SNPable_extractKmers.txt $SCR/SNPable_extractKmers.sh $REFERENCE $SNPABLE


## Step 2: align all reads (produced in step 1) to the genome with BWA:        
for FOCALFILE in $(ls kmers/???)
do
	bsub -n 8 -R "span[hosts=1]" -M 16 -q hour -N -o slurm.SNPable_aln.$FOCALFILE.txt $SCR/SNPable_aln.sh $FOCALFILE $SNPABLE $REFERENCE
done


## Step 3: generate mask (first rawMask, then final mask), then convert to bed:
k_value=50
r_value=0.5
bsub -q day -N -o slurm.SNPable_makeMask.txt $SCR/SNPable_makeMask.sh $SNPABLE $k_value $r_value
rm kmer_alignments/*sam
