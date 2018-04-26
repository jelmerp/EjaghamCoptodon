## Prepare Saguaro input from VCF:
./VCF2HMMFeature -i EjaC.Cgal.Cgui.DP5.GQ20.MAXMISS0.5.MAF0.01.vcf -o cichlids_eja.feature

## Run Saguaro:
./Saguaro -f cichlids_eja.feature -o cichlids_eja -iter 30