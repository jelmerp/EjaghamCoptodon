###### RUN FOR PRE-DEFINED BLOCKS #####
module load r
FILE_ID=EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01
BLOCK_P=0.001
Rscript scripts/windowstats/fd_process3_runHC.R $FILE_ID $BLOCK_P > hybridcheck_knownblocks.txt