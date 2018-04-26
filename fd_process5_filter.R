##### SET-UP #####

## Files to read:
fdblocks.HCtest.file <- 'analyses_output/admixblocks/fdblocks.HCtest.txt'
gphocs_output.summary <- 'analyses_output/gphocs/gphocs_summary_theta_tau.txt'

## Files to write:
admixblocks.all.file <- 'analyses_output/admixblocks/admixblocks.all.txt'
admixblocks.lik.file <- 'analyses_output/admixblocks/admixblocks.likely.txt'
admixblocks.hic.file <- 'analyses_output/admixblocks/admixblocks.highconfidence.txt'

## Libraries:
library(tidyverse)


##### READ TABLES #####
blocks.unfiltered <- read.table(fdblocks.HCtest.file, header = TRUE, as.is = TRUE)
tau <- read.table(gphocs_output.summary, header = TRUE, as.is = TRUE)


##### PROCESS #####
## Get ages from gphocs analyses:
tau <- tau %>% dplyr::filter(var == 'tau')
root.split.point <- dplyr::filter(tau, pop == 'root', migtoType3 == 'all')$mean
root.split.min <- min(dplyr::filter(tau, pop == 'root')$min)

## First filter by nr of windows etc:
blocks.all <- filter(blocks.unfiltered, nwin > 4, consistent != 0, singleHC == 0, p.hc < 0.001)

## Then filter by ages:
blocks.lik <- filter(blocks.all, age.mean < root.split.point)
blocks.hic <- filter(blocks.all, age.max < root.split.min)


##### WRITE TABLES #####
write.table(blocks.all, admixblocks.all.file, sep = '\t', quote = FALSE, row.names = FALSE)
write.table(blocks.lik, admixblocks.lik.file, sep = '\t', quote = FALSE, row.names = FALSE)
write.table(blocks.hic, admixblocks.hic.file, sep = '\t', quote = FALSE, row.names = FALSE)
