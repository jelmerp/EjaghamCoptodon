##### INFER NETWORK WITH ML AND BOOTSTRAPS #####
## Ckot as outgroup:
FILE_ID=EjaC.Ckot.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb
TREEFILE=analyses_output/raxml/${FILE_ID}_allGeneTrees.tre
TAXONMAP="Tgui:TguiNG2,TguiNG5; Cgal:SgalMA1; Ckot:Ckot383,Ckot499; Cfus:Cfus085,Cfus350,Cfus503; Cdec:Cdec088,Cdec328; Ceja:Ceja262,Ceja408"
NMIG=0
NEXUS_IN=analyses_input/phylonet/$FILE_ID.NMIG$NMIG.MLbootstrap.nexus # To be created
NEXUS_OUT=analyses_output/phylonet/$FILE_ID.NMIG$NMIG.MLbootstrap
PHYLONET_CMD="InferNetwork_ML_Bootstrap (all) $NMIG -pl 8 -di"
bsub -n 8 -R "span[hosts=1]" -q week -o slurm.phylonet.MLbootstrap.$FILE_ID.$NMIG scripts/trees/phylonet_run.sh "$TREEFILE" "$NEXUS_IN" "$NEXUS_OUT" "$PHYLONET_CMD" "$TAXONMAP"

## Tgui as outgroup:
FILE_ID=EjaC.Dstat.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb
TREEFILE=analyses_output/raxml/${FILE_ID}_allGeneTrees.tre
TAXONMAP="Tgui:TguiNG2,TguiNG5; Cgal:SgalMA1; Sgui:TguiMA1,TguiMA2,TguiMA4; Cfus:Cfus085,Cfus350,Cfus503; Cdec:Cdec088,Cdec328; Ceja:Ceja262,Ceja408"
NMIG=0
NEXUS_IN=analyses_input/phylonet/$FILE_ID.NMIG$NMIG.MLbootstrap.nexus # To be created
NEXUS_OUT=analyses_output/phylonet/$FILE_ID.NMIG$NMIG.MLbootstrap
PHYLONET_CMD="InferNetwork_ML_Bootstrap (all) $NMIG -pl 8 -di"
bsub -n 8 -R "span[hosts=1]" -q week -o slurm.phylonet.MLbootstrap.$FILE_ID.$NMIG scripts/trees/phylonet_run.sh "$TREEFILE" "$NEXUS_IN" "$NEXUS_OUT" "$PHYLONET_CMD" "$TAXONMAP"


##### INFER SPECIES TREE #####
## Ckot as outgroup:
FILE_ID=EjaC.Dstat.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb
TREEFILE=analyses_output/raxml/${FILE_ID}_allGeneTrees.tre
TAXONMAP="Tgui:TguiNG2,TguiNG5; Cgal:SgalMA1; Sgui:TguiMA1,TguiMA2,TguiMA4; Cfus:Cfus085,Cfus350,Cfus503; Cdec:Cdec088,Cdec328; Ceja:Ceja262,Ceja408"
NEXUS_IN=analyses_input/phylonet/$FILE_ID.ST_MDC.nexus # To be created
NEXUS_OUT=analyses_output/phylonet/$FILE_ID.ST_MDC
PHYLONET_CMD="Infer_ST_MDC (all)"
scripts/trees/phylonet_run.sh "$TREEFILE" "$NEXUS_IN" "$NEXUS_OUT" "$PHYLONET_CMD" "$TAXONMAP"

## Tgui as outgroup:
FILE_ID=EjaC.Ckot.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb
TREEFILE=analyses_output/raxml/${FILE_ID}_allGeneTrees.tre
TAXONMAP="Tgui:TguiNG2,TguiNG5; Cgal:SgalMA1; Ckot:Ckot383,Ckot499; Cfus:Cfus085,Cfus350,Cfus503; Cdec:Cdec088,Cdec328; Ceja:Ceja262,Ceja408"
NEXUS_IN=analyses_input/phylonet/$FILE_ID.ST_MDC.nexus # To be created
NEXUS_OUT=analyses_output/phylonet/$FILE_ID.ST_MDC
PHYLONET_CMD="Infer_ST_MDC (all)"
scripts/trees/phylonet_run.sh "$TREEFILE" "$NEXUS_IN" "$NEXUS_OUT" "$PHYLONET_CMD" "$TAXONMAP"



