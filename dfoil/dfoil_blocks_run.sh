##### SET-UP #####
module load python/2.7.1

FILE_ID=$1

date
echo "Script: dfoil_run.sh"
echo "File ID: $FILE_ID"

## Other variables:
DFOIL_INFILES=$(ls analyses_input/dfoil/$FILE_ID.*dfoil.in*)
DFOIL_OUTFILES=$(echo ${DFOIL_INFILES[@]} | sed 's/input/output/g' | sed 's/dfoil.in/doil.out/g')

## Dfoil scripts:
DFOIL=software/dfoil/dfoil.py
DFOIL_ANALYZE=software/dfoil/dfoil_analyze.py

## Run DFOIL:
$DFOIL --infile ${DFOIL_INFILES[@]} --out ${DFOIL_OUTFILES[@]} --mode dfoilalt --mincount 5 --mintotal 25

echo "Done with script."
date