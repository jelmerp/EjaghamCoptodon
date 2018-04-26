
########################################################################################################################################	
##### RUN RAXML FOR WHOLE GENOME #####

## Global variables:
MODEL=GTRCAT
BOOTSTRAP=100
ASC_COR=FALSE
OUTDIR=analyses_output/raxml

## Ckot as outgroup:
FILE_ID=EjaC.Ckot.DP5.GQ20.MAXMISS0.5.MAF0.01
OUTGROUPS="-o Ckot383,Ckot499"
INPUT=seqdata/fasta/$FILE_ID.fasta
bsub -q day -M 8 -o slurm.raxml.$FILE_ID.$MODEL.$BOOTSTRAP scripts/raxml/raxml_run.sh $INPUT $FILE_ID $OUTDIR $MODEL $ASC_COR $BOOTSTRAP "$OUTGROUPS"

## Tgui as outgroup:
FILE_ID=EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01
OUTGROUPS="-o TguiMA1,TguiMA2,TguiMA4"
INPUT=seqdata/fasta/$FILE_ID.fasta
bsub -q day -M 12 -o slurm.raxml.$FILE_ID.$MODEL.$BOOTSTRAP scripts/raxml/raxml_run.sh $INPUT $FILE_ID $OUTDIR $MODEL $ASC_COR $BOOTSTRAP "$OUTGROUPS"

## Ckot and Tgui as outgroups:
FILE_ID=EjaC.Ckot.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01
OUTGROUPS="-o TguiMA1,TguiMA2,TguiMA4"
INPUT=seqdata/fasta/$FILE_ID.fasta
bsub -q day -M 12 -o slurm.raxml.$FILE_ID.$MODEL.$BOOTSTRAP scripts/raxml/raxml_run.sh $INPUT $FILE_ID $OUTDIR $MODEL $ASC_COR $BOOTSTRAP "$OUTGROUPS"
