## Global settings:
FILE_ID=EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01
MINQUAL=20
MINDEPTH=5

## Files to read:
VCF_FILE=seqdata/$FILE_ID.vcf.gz
POPFILE=analyses_input/admixblocks/popIDs.fd.txt

## Files to write:
GENOFILE=seqdata/variants_otherFormats/geno/$FILE_ID.geno.gz

## Convert VCF to "geno" file:
scripts/conversion/vcf2geno.sh $VCF $GENO_FILE $MINQUAL $MINDEPTH

## Run Simon Martin's ABBA-BABA by window script:
MAXMISS=0.5
WINSIZE=50000
STEPSIZE=5000
MINSITES=50
MINDATA=0.5
scripts/windowstats/fd_submaster.sh $POPFILE $GENOFILE $MINDEPTH $MINQUAL $WINSIZE $STEPSIZE $MINSITES $MINDATA
