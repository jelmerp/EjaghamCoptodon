## Software:
MPEST=software/mpest_1.5/src/mpest

## Files to read:
FILE_ID=raxml_EjaC.Ckot.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb
INFILE_GENETREES=analyses_output/raxml/${FILE_ID}_allGeneTrees.tre

CONTROLFILE=analyses_input/mpest/mpest.controlfile_$FILE_ID

## Files to write:
OUTFILE=analyses/mpest/output/${FILE_ID}_mpest.tre

## Run MP-EST:
cp $INFILE_GENETREES trees.tre
$MPEST $CONTROLFILE
mv trees.tre.tre $OUTFILE


## Controlfile format:
# trees.tre # the name of the gene tree file
# 0  # 1: calculate triple distance among trees. 0: donot calculate
# 6950387  # seed
# 5  # number of independent runs
# 8740 6 # number of genes and number of species
# Cdec 2 Cdec088 Cdec328
# Ceja 2 Ceja262 Ceja408
# Cfus 3 Cfus085 Cfus350 Cfus503
# Tgui 2 TguiNG2 TguiNG5
# Cgal 1 SgalMA1
# Ckot 2 Ckot383 Ckot499
# 0 # usertree provided
