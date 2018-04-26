## Global variables:
FILE_ID=EjaC.Cgal.Cgui.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb_allGeneTrees.tre

## Files to read:
INPUT=analyses_output/raxml/$FILE_ID
INDFILE=analyses_input/astral/indfile_EjaC.Cgal.Cgui.txt

## Files to write:
OUTPUT=analyses_output/astral/output/$FILE_ID

## Software:
ASTRAL=/proj/cmarlab/users/jelmer/software/ASTRAL-multiind/Astral/astral.5.4.3.jar

## Run Astral:
java -jar $ASTRAL -i $INPUT -o $OUTPUT -a $INDFILE