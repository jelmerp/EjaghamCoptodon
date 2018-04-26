## Set-up:
library(plyr)

## Files to read:
dd.filename <- 'analyses/dfoil/intermediateFiles/fdblocks.dfoil.overlap.txt'
tau.filename <- 'analyses/gphocs/output_summaries/theta_tau.txt'

## Files to write:
dfoil.lik.file <- 'analyses_output/admixblocks/admixblocks.dfoil.likely.txt'
dfoil.hic.file <- 'analyses_output/admixblock/admixblocks.dfoil.highconfidence.txt'

## Read files:
dfoil.unfiltered <- read.table(dd.filename, header = TRUE, as.is = TRUE)
dfoil.unfiltered$lakepoptype <- mapvalues(dfoil.unfiltered$lakepop,
                                          from = c('D', 'E', 'F', 'DE', 'DEF'),
                                          to = c(rep('extant species', 3), 'ancestor DE', 'ancestor DEF'))
tau <- read.table(tau.filename, header = TRUE, as.is = TRUE)

## Get ages from gphocs analyses:
tau <- tau %>% dplyr::filter(var == 'tau')
root.split.point <- dplyr::filter(tau, pop == 'root', migtoType3 == 'all')$mean
root.split.min <- min(dplyr::filter(tau, pop == 'root')$min)

## First filter by significance level, etc:
dfoil.all <- filter(dfoil.unfiltered, consistent.fd != 0, singleHC == 0, p.hc < 0.001)

## Then filter by ages:
dfoil.lik <- filter(dfoil.all, age.mean < root.split.point)
dfoil.hic <- filter(dfoil.all, age.max < root.split.min)

## Write tables:
write.table(dfoil.lik, dfoil.lik.file, sep = '\t', quote = FALSE, row.names = FALSE)
write.table(dfoil.hic, dfoil.hic.file, sep = '\t', quote = FALSE, row.names = FALSE)
