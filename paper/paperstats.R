##### SET-UP #####

## Other scripts to load:
source('scripts/gphocs/gphocs_analyze_fun.R')
source('scripts/admixblocks/admixblocks_plotfun.R')

## Files to read:
dp.file  <- 'analyses_output/coverage.idepth'

saguaro_summary.file <- 'analyses/output/saguaro_output_summary.txt'

gphocs_sumr.tt.file <- 'analyses_output/gphocs_summary_theta_tau.txt'
gphocs_sumr.m.file <- 'analyses_output/gphocs_summary_m.txt'

fdblocks.file <- 'analyses_output/admixblocks_fd.only_EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01.txt'
admixblocks.highconfidence.file <- 'analyses_output/admixblocks.highconfidence.txt'

dfoil_admixblocks.likely.file <- 'analyses/dfoil/admixblocks.dfoil.likely.txt'
dfoil_admixblocks.highconfidence.file <- 'analyses/dfoil/admixblocks.dfoil.highconfidence.txt'


##### METADATA: SEQUENCING DEPTH #####
dp <- read.table(dp.file, header = TRUE)
mean(dp$MEAN_DEPTH)
arrange(dp, MEAN_DEPTH)


##### TREES: SAGUARO #####
cac <- read.table(saguaro_summary.file, header = TRUE, stringsAsFactors = FALSE)
sum(cac$percentage) # 100

sum(filter(cac, EjaC.clade.mono == 1)$percentage) # 90.01968 --> EjaC monophyletic as clade
sum(filter(cac, EjaC.species.mono == 1)$percentage) # 87.61261 --> EjaC species monophyletic
sum(filter(cac, EjaC.species.para == 1)$percentage) # 9.2794 --> EjaC species paraphyletic
sum(filter(cac, Mam.Gui.mono == 1)$percentage) # 90.01968 --> Mam+Gui mono


##### GPHOCS #####

## Summaries for migration from Mam:
sumr.tt <- dplyr::filter(Log, var %in% c('tau', 'theta'), migfromRun %in% c('none', 'Mam')) %>%
  dplyr::group_by(migfromRun, migtoType3, var, pop) %>%
  dplyr::summarise(mean = mean(cval), min = hpd.min(cval), max = hpd.max(cval))
write.table(sumr.tt, gphocs_sumr.tt.file, sep  = '\t', quote = FALSE, row.names = FALSE)

sumr.m1 <- dplyr::filter(Log, var %in% c('m.prop', '2Nm'), migfromRun == 'Mam') %>%
  dplyr::group_by(migfromRun, migtoType3, migto, var) %>%
  dplyr::summarise(mean = mean(val), min = hpd.min(val), max = hpd.max(val))
sumr.m2 <- dplyr::filter(Log, var == 'm', migfromRun == 'Mam') %>%
  dplyr::group_by(migfromRun, migtoType3, migto, var) %>%
  dplyr::summarise(mean = mean(cval), min = hpd.min(cval), max = hpd.max(cval))
sumr.m <- rbind(sumr.m1, sumr.m2) %>% arrange(migfromRun, migtoType3, migto, var)
write.table(sumr.m, gphocs_sumr.m1.file, sep  = '\t', quote = FALSE, row.names = FALSE)

tau <- sumr.tt %>% dplyr::filter(var == 'tau')
(root.split.point <- dplyr::filter(tau, pop == 'root', migtoType3 == 'all')$mean) # 9760.314
(root.split.min <- min(dplyr::filter(tau, pop == 'root')$min)) # 6587.2

## Summaries for migration from Gui:
m.prop.Gui <- dplyr::filter(Log, migfrom == 'Gui', var == 'm.prop') %>%
  dplyr::group_by(migtoType3, migto, var) %>%
  dplyr::summarise(mean = mean(val), min = hpd.min(val), max = hpd.max(val))
m.2Nm.Gui <- dplyr::filter(Log, migfrom == 'Gui', var == '2Nm') %>%
  dplyr::group_by(migtoType3, migto, var) %>%
  dplyr::summarise(mean = mean(val), min = hpd.min(val), max = hpd.max(val))
m.tot.Gui <- dplyr::filter(Log, migfrom == 'Gui', var == 'm') %>%
  dplyr::group_by(migtoType3, migto, var) %>%
  dplyr::summarise(mean = mean(cval), min = hpd.min(cval), max = hpd.max(cval))

## Summaries for migration AU -> DEF:
m.2Nm.AU <- dplyr::filter(Log, migfrom == 'AU', var == '2Nm') %>%
  dplyr::group_by(migtoType3, migto, var) %>%
  dplyr::summarise(mean = mean(val), min = hpd.min(val), max = hpd.max(val))

## Summaries for intra-radiation gene flow:
m.2Nm.wi <- filter(Log, var == '2Nm' & migfromRun %in% c('DE', 'Dec', 'Eja', 'Fus') & migto != 'Gui' & migto != 'Mam') %>%
  dplyr::group_by(migtoType3, migfrom, migto, var) %>%
  dplyr::summarise(mean = mean(val), min = hpd.min(val), max = hpd.max(val))

## Migration to Cfus compared to Ceja:
m.tot <- sumr.m2
Fus.migAll <- filter(m.tot, migto == 'Fus', migtoType3 == 'all')$mean
Eja.migAll <- filter(m.tot, migto == 'DE', migtoType3 == 'all')$mean + filter(m.tot, migto == 'Eja', migtoType3 == 'all')$mean
Dec.migAll <- filter(m.tot, migto == 'DE', migtoType3 == 'all')$mean + filter(m.tot, migto == 'Dec', migtoType3 == 'all')$mean
(Eja.migAll - Fus.migAll) / Eja.migAll
(Dec.migAll - Fus.migAll) / Dec.migAll

## Difference since split Fus - DE:
Fus.migSingle <- filter(m.tot, migto == 'Fus', migtoType3 == 'single')$mean
Eja.migSingle <- filter(m.tot, migto == 'DE', migtoType3 == 'single')$mean + filter(m.tot, migto == 'Eja', migtoType3 == 'single')$mean
Dec.migSingle <- filter(m.tot, migto == 'DE', migtoType3 == 'single')$mean + filter(m.tot, migto == 'Dec', migtoType3 == 'single')$mean
(Eja.migSingle - Fus.migSingle) / Eja.migSingle # 0.4321928 less migration to Cfus than Ceja since split of Cfus from DE
(Dec.migSingle - Fus.migSingle) / Dec.migSingle # 0.4062653 less migration to Cfus than Cdec since split of Cfus from DE

## Difference since split root - DEF:
DEF.mig <- filter(m.tot, migto == 'DEF', migtoType3 == 'single')$mean
Fus.migTot <- Fus.migSingle + DEF.mig
Eja.migTot <- Eja.migSingle + DEF.mig
Dec.migTot <- Dec.migSingle + DEF.mig
(Eja.migTot - Fus.migTot) / Eja.migTot
(Eja.migTot - Dec.migTot) / Eja.migTot
(Dec.migTot - Fus.migTot) / Fus.migTot


##### ADMIXTUREBLOCKS ######
fdblocks <- read.table(fdbocks.file, header = TRUE)
blocks.hic <- read.table(blocks.highconfidence.file, header = TRUE, as.is = TRUE)

## "Best" blocks:
grep('NC_022214.1_6220001_FEA|NC_022214.1_6215001_FDA', blocks.hic$ID)
arrange(blocks.hic, desc(fd.sum))$ID %>% grep(pattern = 'NC_022214.1_6220001_FEA|NC_022214.1_6215001_FDA') # 1,2
arrange(blocks.hic, p.min)$ID %>% grep(pattern = 'NC_022214.1_6220001_FEA|NC_022214.1_6215001_FDA') # 33,34
arrange(blocks.hic, p.hc)$ID %>% grep(pattern = 'NC_022214.1_6220001_FEA|NC_022214.1_6215001_FDA') # 2, 5
arrange(blocks.hic, desc(nwin))$ID %>% grep(pattern = 'NC_022214.1_6220001_FEA|NC_022214.1_6215001_FDA') # 1, 2

## All blocks prior to overlap with Hybridcheck:
fdblocks <- filter(fdblocks, nwin > 4, consistent != 0)
nrow(fdblocks) # 1,138

## All blocks after overlap with Hybridcheck:
## Likely blocks:
nrow(blocks.lik) # 259
table(blocks.lik$riverpop)
#Cgui Cmam
#82  177

## High-confidence blocks:
nrow(blocks.hic) # 146
table(blocks.hic$riverpop)
#Cgui Cmam
#29  117

## Mean age for blocks:
mean(filter(blocks.lik, riverpop == 'Cgui')$age.mean) # 4555
mean(filter(blocks.lik, riverpop == 'Cmam')$age.mean) # 2938

mean(filter(blocks.hic, riverpop == 'Cgui')$age.mean) # 1972.517
mean(filter(blocks.hic, riverpop == 'Cmam')$age.mean) # 1345.216

table(blocks.lik$unq)
table(blocks.hic$unq)

## Dfoil blocks:
dfoil.lik <- read.delim(dfoil_admixblocks.likely.file, header = TRUE, as.is = TRUE)
dfoil.hic <- read.delim(dfoil_admixblocks.highconfidence.file, header = TRUE, as.is = TRUE)

table(dfoil.lik$lakepoptype)
table(dfoil.hic$lakepoptype)

table(dfoil.hic$riverpop)

table(dfoil.hic$ancDir)
table(dfoil.hic$migDir)
table(dfoil.hic$migpattern)