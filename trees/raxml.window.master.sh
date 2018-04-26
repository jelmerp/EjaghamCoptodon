##### NEEDED #####

## Scripts:
scripts/trees/raxml.window.sub.sh
scripts/conversion/splitVCF_byCoords.sh
scripts/conversion/vcf2fasta.sh
scripts/conversion/vcf_tab_to_fasta_alignment.pl
scripts/trees/raxml_run.sh
/proj/cmarlab/users/jelmer/software/standard-RAxML-master/raxmlHPC-SSE3 # Change location of RaxML program on line 11 of raxml_run.sh

## Folders:
seqdata/vcf_singleScaf # Where cut-up VCFs will go (can be changed in scripts/raxml.window.sub.sh)
seqdata/variants_otherFormats/fasta/ # Where fasta files will go (can be changed in scripts/raxml.window.sub.sh)
analyses/raxml # Where stats on nr of sites per window will go (can be changed in scripts/raxml.window.sub.sh)
analyses/raxml/output/byWindow # Where raxml output will go (can be changed in this script)
seqdata/vcf_split # Where input VCF is located (can be changed in this script)
metadata # Where scaffold info is and file with windows will go (can be changed in this script)

## Files:
metadata/scaffolds_withLength.txt # Info on scaffolds; should have columns "scaffold.name" and "scaffold.length"
$INDIR/$FILE_ID.vcf OR $INDIR/$FILE_ID.vcf.gz # input VCF file, define INDIR AND FILE_ID below


##### CREATE FILE WITH ONE WINDOW PER LINE #####
DIR=metadata
WINSIZE=100000 # Window size
STEPSIZE=100000 # Step size
MIN_SCAFFOLD_SIZE=500000 # Minimum scaffold size, don't create windows for smaller scaffolds
scripts/misc/scaffoldwindows.R $DIR $WINSIZE $STEPSIZE $MIN_SCAFFOLD_SIZE  


##### RUN RAXML FOR EACH WINDOW #####
FILE_ID=EjaC.Ckot.DP5.GQ20.MAXMISS0.9.MAF0.01
INDIR=seqdata/vcf_split
OUTDIR=analyses/raxml/output/byWindow
WINDOW_FILE=$DIR/scaffoldWindows_win${WINSIZE}_step${STEPSIZE}.txt
MODEL=GTRCAT
MINSITES=250
SKIPFASTA=FALSE
OUTGROUPS="-o Ckot383,Ckot499"

for WINDOW_LINE in $(seq 1 $(cat $WINDOW_FILE | wc -l))
do
	bsub -q hour -o slurm.raxml.window.$FILE_ID.$MINSITES.sub.txt scripts/trees/raxml.window.sub.sh $FILE_ID $INDIR $OUTDIR $WINDOW_FILE $WINDOW_LINE $MODEL $MINSITES $SKIPFASTA "$OUTGROUPS"
done
