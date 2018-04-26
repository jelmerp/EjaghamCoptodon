##### GLOBAL VARIABLES #####
SEQ_ID=EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01

##### PREP INPUT #####
VCF=seqdata/$FILE_ID.vcf.gz
INDFILE=analyses_input/admixtools/indfile_$FILE_ID.txt
MAF=0.01
bsub -q day -o slurm.admixtools.prepinput.$SEQ_ID scripts/admixtools/admixtools_prepInput.sh $SEQ_ID $VCF $INDFILE $MAF


##### RUN ADMIXTOOLS QPDSTAT #####
POPFILE_ID=EjaC.Dstat
INDFILE_ID=EjaC.Dstat
POPFILE=analyses/admixtools/input/popfile_$POPFILE_ID.txt

for LINE in $(seq 1 $(cat $POPFILE | wc -l)
do
	sbatch -t 2:00:00 -n 1 --mem-per-cpu=24G -o slurm.admixtools.$SEQ_ID.$LINE scripts/admixstats/admixtools_run.sh $SEQ_ID $POPFILE_ID $INDFILE_ID $LINE
done

grep -h "result" analyses/admixtools/output/raw/*$SEQ_ID.$POPFILE_ID.line*dmode* > analyses_output/admixstats/$SEQ_ID.$POPFILE_ID.dmode.out
grep -h "result" analyses/admixtools/output/raw/*$SEQ_ID.$POPFILE_ID.line*fmode* > analyses_output/admixstats/$SEQ_ID.$POPFILE_ID.fmode.out


##### RUN ADMIXTOOLS QPF4RATIO #####
SEQ_ID=EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01
POPFILE_ID=EjaC.Dstat
sbatch -t 1:00:00 --mem=24G -e error.f4ratio.$SEQ_ID -o slurm.f4ratio.$SEQ_ID -J atools.f4ratio.$SEQ_ID scripts/admixstats/admixtools_runf4ratio.sh $SEQ_ID $POPFILE_ID

grep -h "result" analyses/admixtools/output/raw/*$SEQ_ID*f4ratio*out > analyses_output/admixstats/$SEQ_ID.f4ratio.out