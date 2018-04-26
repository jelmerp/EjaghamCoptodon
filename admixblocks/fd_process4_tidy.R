##### SET-UP ####

## Files to read:
fdblocks.file.in <- 'analyses_output/admixblocks/fdblocks.txt'
blocks.hc.file.in <- 'analyses_output/admixblocks/fdblocks.HCtest_intermed.txt'

## Files to write:
fdblocks.file.out <- 'analyses_output/admixblocks/fdblocks.HCtest.txt'
fdblocks.forDfoil.file <- 'analyses_output/admixblocks/fdblocks.forDFOIL.txt'

## Libraries:
library(dplyr)
file.id <- 'EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01'


##### READ DATA #####
fdblocks <- read.table(fdblocks.file.in, header = TRUE, as.is = TRUE)
fdblocks <- dplyr::arrange(fdblocks, desc(fd.sum))
fd.hc <- read.table(blocks.hc.file.in, header = TRUE, as.is = TRUE)


##### POST-PROCESSING #####
## Add column specifying whether HC identified block in one or both individuals:
frequencies <- table(fd.hc$ID)
unq <- names(frequencies[frequencies == 1])
fd.hc$singleHC <- 0
fd.hc$singleHC[match(unq, fd.hc$ID)] <- 1

## Get mean of stats for HC blocks for each individual:
smr.tmp <- fd.hc %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(p.hc = mean(p.hc), age.min = round(mean(age.min)),
                   age.mean = round(mean(age.mean)), age.max = round(mean(age.max)))
fd.hc <- merge(smr.tmp, fd.hc[, -which(names(fd.hc) %in% c('p.hc', 'age.min', 'age.mean', 'age.max'))],
               by = 'ID', all.y = FALSE)

## Edit some columns:
fd.hc <- fd.hc[-which(duplicated(fd.hc$ID)), ]
fd.hc <- merge(fdblocks, fd.hc[, c('ID', 'duo.hc', 'p.hc', 'age.min', 'age.mean', 'age.max', 'singleHC')], by = 'ID')
fd.hc$St <- round(fd.hc$start / 1000000, 3)
fd.hc$En <- round(fd.hc$end / 1000000, 3)

fd.hc$triplet <- gsub('Cdec', 'Cdec', fd.hc$triplet)
fd.hc$triplet <- gsub('Ceja', 'Ceja', fd.hc$triplet)
fd.hc$triplet <- gsub('Cfus', 'Cfus', fd.hc$triplet)
fd.hc$triplet <- gsub('Cgui', 'Cgui', fd.hc$triplet)
fd.hc$triplet <- gsub('Cmam', 'Cmam', fd.hc$triplet)
fd.hc$riverpop <- ifelse(fd.hc$popA %in% c('Cmam', 'Cgui'), fd.hc$popB, fd.hc$popC)
fd.hc$lakepop <- ifelse(fd.hc$popA %in% c('Cdec', 'Ceja', 'Cfus'), fd.hc$popB, fd.hc$popC)
fd.hc$type <- ifelse(fd.hc$popA %in% c('Cdec', 'Ceja', 'Cfus'), 'lake.diff', 'river.diff')

## Manually edit some blocks:
fd.hc$unq <- factor(fd.hc$unq, levels = c('unq', 'shared1', 'shared2'))
fd.hc$unq[fd.hc$ID == 'NC_022219.1_5450001_DFA'] <- 'shared2'
fd.hc$unq[fd.hc$ID == 'NC_022213.1_18320001_AUE'] <- 'shared2'
fd.hc$unq[fd.hc$ID == 'NC_022214.1_5330001_UAE'] <- 'shared2'
fd.hc$unq[fd.hc$ID == 'NC_022217.1_21905001_FEA'] <- 'shared2'
fd.hc$unq[fd.hc$ID == 'NT_167529.1_605001_AUD'] <- 'shared1'
to.remove.bad <- c('NC_022217.1_24340001_DEU')
to.remove.dup <- c('NC_022205.1_10001_DEU', 'NC_022218.1_17730001_FDU', 'NC_022219.1_5455001_UAF')
fd.hc <- fd.hc[-which(fd.hc$ID %in% c(to.remove.bad, to.remove.dup)), ]

##### WRITE DATA #####
write.table(fd.hc, blocks.file, sep = '\t', quote = FALSE, row.names = FALSE)

## Select only a few columns for DFOIL:
fd.hc.dfoil <- select(fd.hc, scaffold, start, end)
write.table(fd.hc.dfoil, fdblocks.forDfoil.file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)