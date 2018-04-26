#!/bin/bash -l

module load python/2.7.1

FILE_ID=$1
WINDOW_FILE=$2
WINDOW_LINE=$3
OUTGROUP=$4

VCF_IN=seqdata/vcf_split/$FILE_ID

date
echo "Script: dfoil_inputprep.sh"
echo "File ID: $FILE_ID"
echo "Window file: $WINDOW_FILE"
echo "Window line: $WINDOW_LINE"
echo "VCF input file: $VCF_IN"


##### STEP 1 - SPLIT VCF AND CONVERT TO FASTA #####
SCAFFOLD_M=$(cat $WINDOW_FILE | sed -ne "$WINDOW_LINE,${WINDOW_LINE}p;${WINDOW_LINE}q" | cut -f 1 | sed 's/^[ \t]*//;s/[ \t]*$//')
START_M=$(cat $WINDOW_FILE | sed -ne "$WINDOW_LINE,${WINDOW_LINE}p;${WINDOW_LINE}q" | cut -f 2 | sed 's/^[ \t]*//;s/[ \t]*$//')
END_M=$(cat $WINDOW_FILE | sed -ne "$WINDOW_LINE,${WINDOW_LINE}p;${WINDOW_LINE}q" | cut -f 3 | sed 's/^[ \t]*//;s/[ \t]*$//')
INDIR=seqdata/vcf_split

echo "Cutting up vcfs in windows... (Scaffold: $SCAFFOLD_M, Start: $START_M, End: $END_M)"
scripts/conversion/splitVCF_byCoords.sh $FILE_ID $SCAFFOLD_M $START_M $END_M $INDIR
	
FILE_ID_FULL=$FILE_ID.$SCAFFOLD_M.$START_M.$END_M
INDIR_SPLIT=seqdata/vcf_singleScaf
echo "File ID - split VCF: $FILE_ID_FULL"
	
echo "Converting vcf to fasta..."
scripts/conversion/vcf2fasta.sh $FILE_ID_FULL $INDIR_SPLIT ALL

FASTA=seqdata/fasta/$FILE_ID_FULL.fasta


##### STEP 2 - EXCLUDE ONE IND IN THREE CONFIGS #####
echo "Creating three fasta configurations..."

## A - configuration without Cfus:
cat $FASTA | sed 's/Cdec[0-9][0-9][0-9]/Cdec/' | sed 's/Ceja[0-9][0-9][0-9]/Ceja/' | sed 's/Cfus[0-9][0-9][0-9]/Cfus/' | sed 's/SgalMA1/Cmam/' \
| sed 's/TguiNG[0-9]/Cgui/' | sed 's/TguiMA[0-9]/Sgui/' | sed 's/Ckot[0-9][0-9][0-9]/Ckot/' | awk '/^>/ {P=index($0,"Cfus")==0} {if(P) print} ' > $FASTA.noCfus.dfoil

## B - configuration without Cdec:
cat $FASTA | sed 's/Cdec[0-9][0-9][0-9]/Cdec/' | sed 's/Ceja[0-9][0-9][0-9]/Ceja/' | sed 's/Cfus[0-9][0-9][0-9]/Cfus/' | sed 's/SgalMA1/Cmam/' \
| sed 's/TguiNG[0-9]/Cgui/' | sed 's/TguiMA[0-9]/Sgui/' | sed 's/Ckot[0-9][0-9][0-9]/Ckot/' | awk '/^>/ {P=index($0,"Cdec")==0} {if(P) print} ' > $FASTA.noCdec.dfoil

## C - configuration without Ceja:
cat $FASTA | sed 's/Cdec[0-9][0-9][0-9]/Cdec/' | sed 's/Ceja[0-9][0-9][0-9]/Ceja/' | sed 's/Cfus[0-9][0-9][0-9]/Cfus/' | sed 's/SgalMA1/Cmam/' \
| sed 's/TguiNG[0-9]/Cgui/' | sed 's/TguiMA[0-9]/Sgui/' | sed 's/Ckot[0-9][0-9][0-9]/Ckot/' | awk '/^>/ {P=index($0,"Ceja")==0} {if(P) print} ' > $FASTA.noCeja.dfoil


##### STEP 3 - PREP DFOIL INPUT FROM FASTA #####
DFOIL_INFILE=analyses_input/dfoil/$FILE_ID_FULL.dfoil.in
FASTA2DFOIL=software/dfoil/fasta2dfoil.py

echo "Converting fasta to DFOIL..."
$FASTA2DFOIL $FASTA.noCdec.dfoil --out $DFOIL_INFILE.noCdec --names Ceja Cfus Cgui Cmam $OUTGROUP
$FASTA2DFOIL $FASTA.noCeja.dfoil --out $DFOIL_INFILE.noCeja --names Cdec Cfus Cgui Cmam $OUTGROUP
$FASTA2DFOIL $FASTA.noCfus.dfoil --out $DFOIL_INFILE.noCfus --names Cdec Ceja Cgui Cmam $OUTGROUP


echo "Done with script."
date