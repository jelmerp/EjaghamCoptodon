########################################################################################################################################################################
#### SET-UP ####
IFS=$'\n' read -d '' -a SCAFFOLDS < metadata/scaffolds.txt # List  of scaffolds
## "Controlfiles" need to be manually prepared, not done in these scripts

## Variables:
FILE_ID=EjaC.Cgal.Cgui.DP5.GQ20.MAXMISS0.5.MAF0.01
LOCUS_SIZE=1000
MIN_DISTANCE=50000
MAX_MISSING=50
MASK_CPG=FALSE
MAX_SNPS=25


########################################################################################################################################################################
#### PREP INPUT ####

## 1. Prep loci for each scaffold
for SCAFFOLD in ${SCAFFOLDS[@]:1:501}
	do bsub -q hour -o slurm.gphocs_123wrap.$FILE_ID.$SCAFFOLD.txt scripts/gphocs/gphocs_123wrap.sh $FILE_ID $SCAFFOLD $LOCUS_SIZE $MIN_DISTANCE $MAX_MISSING $MASK_CPG $MAX_SNPS
done

## 2. Combine loci into single gphocs input file
scripts/gphocs/gphocs_4_combineLoci.sh $FILE_ID


########################################################################################################################################################################
#### RUN GPHOCS ####
NPROC=12

for CFILE in $(ls analyses/gphocs/controlfiles/*ctrl)
	do bsub -n $NPROC -R "span[hosts=1]" -q week -o slurm.gphocs_run.$(basename $CFILE).txt scripts/gphocs/gphocs_5_run.sh $CFILE $NPROC
done

