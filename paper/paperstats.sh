
##### NR OF SNPs, COVERAGE #####
VCF=seqdata/EjaC.Ckot.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01.vcf.gz
zcat $VCF | grep -v "##" | wc -l # number of SNPs
vcftools --gzvcf $VCF --depth # mean depth


##### NR OF GENE TREES #####
cat analyses_output/raxml_EjaC.Dstat.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb_allGeneTrees.tre | wc -l # 1766
cat analyses_output/raxml_EjaC.Ckot.Dstat.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb_allGeneTrees.tre | wc -l # 1532
cat analyses_output/raxml_EjaC.Ckot.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb_allGeneTrees.tre | wc -l # 2559


##### GPHOCS #####
less analyses_input/EjaC.Cgal.Cgui.DP5.GQ20.MAXMISS0.5.MAF0.01.gphocsInput.txt # Check nr of loci
