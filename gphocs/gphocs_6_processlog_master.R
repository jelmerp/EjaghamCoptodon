## Other scripts:
source('scripts/gphocs/gphocs_analyze_fun.R')
processScript <- 'scripts/gphocs/gphocs_processlog.R'

## Variables:
gentime <- 1
mutrate.gen <- 7.5e-9
m.scale <- 1000
t.scale <- 0.0001

## Process G-PhoCS output:
commandArgs <- function() c('process', 'analyses_output/gphocs/raw/all/'); source(processScript)

Log <- readRDS('analyses_output/gphocs/output/all/mergedlog.RDS')
write.table(Log, 'analyses_output/gphocs/output/all/mergedlog.RDS', sep = '\t', quote = FALSE, row.names = FALSE)
