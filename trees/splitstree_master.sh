## Nexus files are first created from fasta using PGDSspider, see PGDSspider_master.sh. For vcf-to-fasta conversion, see vcf2fasta.sh 

## Run splitstree:
FILE_ID=EjaC.Ckot.DP5.GQ20.MAXMISS0.5.MAF0.01

NEXUS_IN=seqdata/nexus/$FILE_ID.nexus
NEXUS_OUT=analyses_output/splitstree/splitstree_$FILE_ID.nexus

bsub -q day -M 12 -o slurm.splitstree.$FILE_ID scripts/trees/splitstree_run.sh $NEXUS_IN $NEXUS_OUT
