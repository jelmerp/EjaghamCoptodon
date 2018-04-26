## Variables:
FILE_ID=EjaC.Dstat
OUTGROUP=Sgui

## Files to read:
WINDOW_FILE=analyses_output/admixblocks/fdblocks_forDFOIL.txt # Produced by admixblocks pipeline

## Prep input:
for WINDOW_LINE in $(seq 1 $(cat $WINDOW_FILE | wc -l))
do
	echo $WINDOW_LINE
	bsub -q hour -o slurm.dfoil_inputprep.txt scripts/dfoil/dfoil_inputprep.sh $FILE_ID $WINDOW_FILE $WINDOW_LINE $OUTGROUP
done

## Run DFOIL:
bsub -q day -o slurm.dfoil_run.txt scripts/dfoil/dfoil_run.sh $FILE_ID
