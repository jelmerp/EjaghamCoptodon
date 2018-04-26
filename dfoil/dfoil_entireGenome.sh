##### SET-UP #####
module load python/2.7.12

## Key variables:
FILE_ID=EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01; OUTGROUP=Sgui

## Other variables:
FASTA=seqdata/fasta/$FILE_ID.fasta
DFOIL_INFILE=analyses_input/dfoil/$FILE_ID.dfoil.in
DFOIL_OUTFILE=analyses_output/dfoil/$FILE_ID.dfoil.out

## Dfoil scripts:
FASTA2DFOIL=/proj/cmarlab/users/jelmer/software/dfoil/fasta2dfoil.py
DFOIL=/proj/cmarlab/users/jelmer/software/dfoil/dfoil.py
DFOIL_ANALYZE=/proj/cmarlab/users/jelmer/software/dfoil/dfoil_analyze.py


##### STEP 1 - PREP FASTA #####
## A - configuration without Cfus:
cat $FASTA | sed 's/Cdec[0-9][0-9][0-9]/Cdec/' | sed 's/Ceja[0-9][0-9][0-9]/Ceja/' | sed 's/Cfus[0-9][0-9][0-9]/Cfus/' | sed 's/SgalMA1/Cmam/' \
| sed 's/TguiNG[0-9]/Cgui/' | sed 's/TguiMA[0-9]/Sgui/' | sed 's/Ckot[0-9][0-9][0-9]/Ckot/' | awk '/^>/ {P=index($0,"Cfus")==0} {if(P) print} ' > $FASTA.noCfus.dfoil

## B - configuration without Cdec:
cat $FASTA | sed 's/Cdec[0-9][0-9][0-9]/Cdec/' | sed 's/Ceja[0-9][0-9][0-9]/Ceja/' | sed 's/Cfus[0-9][0-9][0-9]/Cfus/' | sed 's/SgalMA1/Cmam/' \
| sed 's/TguiNG[0-9]/Cgui/' | sed 's/TguiMA[0-9]/Sgui/' | sed 's/Ckot[0-9][0-9][0-9]/Ckot/' | awk '/^>/ {P=index($0,"Cdec")==0} {if(P) print} ' > $FASTA.noCdec.dfoil

## C - configuration without Ceja:
cat $FASTA | sed 's/Cdec[0-9][0-9][0-9]/Cdec/' | sed 's/Ceja[0-9][0-9][0-9]/Ceja/' | sed 's/Cfus[0-9][0-9][0-9]/Cfus/' | sed 's/SgalMA1/Cmam/' \
| sed 's/TguiNG[0-9]/Cgui/' | sed 's/TguiMA[0-9]/Sgui/' | sed 's/Ckot[0-9][0-9][0-9]/Ckot/' | awk '/^>/ {P=index($0,"Ceja")==0} {if(P) print} ' > $FASTA.noCeja.dfoil


##### STEP 2 - PREP DFOIL INPUT FROM FASTA #####
$FASTA2DFOIL $FASTA.noCdec.dfoil --out $DFOIL_INFILE.noCdec --names Ceja Cfus Cgui Cmam $OUTGROUP
$FASTA2DFOIL $FASTA.noCeja.dfoil --out $DFOIL_INFILE.noCeja --names Cdec Cfus Cgui Cmam $OUTGROUP
$FASTA2DFOIL $FASTA.noCfus.dfoil --out $DFOIL_INFILE.noCfus --names Cdec Ceja Cgui Cmam $OUTGROUP

##### STEP 3 - RUN DFOIL #####
$DFOIL --infile $DFOIL_INFILE.noCfus --out $DFOIL_OUTFILE.noCfus --mode dfoilalt
$DFOIL --infile $DFOIL_INFILE.noCdec --out $DFOIL_OUTFILE.noCdec --mode dfoilalt
$DFOIL --infile $DFOIL_INFILE.noCeja --out $DFOIL_OUTFILE.noCeja --mode dfoilalt