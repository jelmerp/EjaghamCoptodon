##### SET-UP #####

## Files to read:
popduos.file <-'analyses_input/sumstats/popduos_sumstats.txt'
poptriplets.file <- 'analyses_input/admixblocks/poptrios_fd.txt'

fd.file.in <- 'analyses_output/admixblocks/fd_intermediate.txt'
fdblocks.file.in <- 'analyses_output/admixblocks/fdblocks_intermediate.txt'

## Files to write:
fd.file.out <- 'analyses_output/admixblocks/fd.txt'
fdblocks.file.out <- 'analyses_output/admixblocks/fdblocks.txt'

## Libraries:
library(data.table)
library(plyr)
library(dplyr)
library(GenomicRanges)

## Read files:
duos <- read.table(popduos.file, as.is = TRUE)
triplets <- read.table(poptriplets.file, header= TRUE, as.is = TRUE)

## Variables:
popnames.long <- c('Cdec', 'Ceja', 'Cfus', 'Cgui', 'Cmam')
popnames.short <- c('D', 'E', 'F', 'U', 'A')

## Triplets:
triplets.long <- paste0(triplets$popA, '.', triplets$popB, '.', triplets$popC)
triplets.short <- gsub('\\.', '', triplets.long)
triplets.short <- gsub('Cfus', 'F', triplets.short)
triplets.short <- gsub('Cdec', 'D', triplets.short)
triplets.short <- gsub('Ceja', 'E', triplets.short)
triplets.short <- gsub('Cmam', 'A', triplets.short)
triplets.short <- gsub('Cgui', 'U', triplets.short)


##### FUNCTIONS #####
overlap.get <- function(rownr, focalblocks, refblocks, min.overlap = 0.5) {

  foc <- focalblocks[rownr, ]
  foc.gr <- as(paste0(foc$scaffold, ':', foc$start, '-', foc$end), "GRanges")
  ref.gr <- as(paste0(refblocks$scaffold, ':', refblocks$start, '-', refblocks$end), "GRanges")

  if(min.overlap < 1) min.overlap <- foc$nwin * min.overlap
  hits <- findOverlaps(foc.gr, ref.gr, minoverlap = min.overlap)
  hits <- refblocks[hits@to, ]
  hits.unique <- unique(c(foc$triplet.short, hits$triplet.short))
  hits.unique <- hits.unique[-which(hits.unique == foc$triplet.short)]

  cat('Row number:', rownr, '   Number of hits:', nrow(hits),
      'Number of unique hits:', length(hits.unique),'\n')

  popmatches <- as.integer(gsub(TRUE, 1, gsub(FALSE, 0, triplets.short %in% hits.unique)))

  return(c(foc$ID, foc$triplet.short, nrow(hits), length(hits.unique), popmatches))
}

overlap.wrap <- function(focalblocks, refblocks, min.overlap = 0.5) {

  overlaps.list <- lapply(1:nrow(focalblocks), overlap.get,
                          focalblocks = focalblocks, refblocks = refblocks, min.overlap = min.overlap)
  overlaps <- rbind.fill(lapply(overlaps.list, function(x){as.data.frame(t(x), stringsAsFactors = FALSE)}))

  colnames(overlaps) <- c('ID', 'triplet.short', 'ov.fd', 'ov.fd.un', triplets.short)
  intcols <- 5:ncol(overlaps)
  overlaps[, intcols] = apply(overlaps[, intcols], 2, function(x) as.integer(x))
  return(overlaps)
}


##### APPLY #####
fd <- as.data.frame(fread(fd.file.in, stringsAsFactors = FALSE))
blocks <- as.data.frame(fread(fdblocks.file.in, stringsAsFactors = FALSE))

## Add block info:
if(any(is.na(blocks$pop))) blocks <- blocks[-which(is.na(blocks$pop)), ]
if(any(blocks$win.first > blocks$win.last)) blocks <- blocks[-which(blocks$win.first > blocks$win.last), ] # REMOVE BLOCKS CROSSING SCAFFOLD BORDERS
blocks <- blocks[-which(blocks$nwin < 3), ]
blocks <- blocks[-which(blocks$fd.max < 0.5), ]

## Formatting:
blocks <- blocks %>% dplyr::rename(start = win.first, end = win.last, triplet = pop)
blocks$recip <- unlist(lapply(strsplit(blocks$triplet, split = '\\.'), '[[', 2))
blocks$donor <- unlist(lapply(strsplit(blocks$triplet, split = '\\.'), '[[', 3))
a <- mapvalues(unlist(lapply(strsplit(blocks$triplet, split = '\\.'), '[[', 1)),
               from = popnames.long, to = popnames.short, warn_missing = FALSE)
b <- mapvalues(blocks$recip, from = popnames.long, to = popnames.short, warn_missing = FALSE)
c <- mapvalues(blocks$donor, from = popnames.long, to = popnames.short, warn_missing = FALSE)
blocks$triplet.short <- paste0(a, b, c)
blocks$ID <- paste0(blocks$scaffold, '_', blocks$start, '_', blocks$triplet.short)

## Remove blocks with Cgui as popA and crater-lake species as popB:
blocks <- dplyr::filter(blocks, triplet.short != 'UDA', triplet.short != 'UEA', triplet.short != 'UFA')

## Remove blocks with start and end window = 1 (why are they there??!!)
blocks <- blocks[-which(blocks$start == 1 & blocks$end == 1), ]

## Add column to fd whether window is in a block:
getBlockWins <- function(blocksrow) {
  # blocksrow <- 1
  block <- blocks[blocksrow, ]
  wins <- seq(block$start, block$end, by = 5000)
  win.IDs <- paste0(block$scaffold, '_', block$triplet, '_', wins)
}
blockWins <- unlist(sapply(1:nrow(blocks), getBlockWins))
winID <- paste0(fd$scaffold, '_', fd$pop, '_', fd$start)
fd$is.block <- FALSE
fd$is.block[match(blockWins, winID)] <- TRUE

## Edit popnames:
fd$pop <- gsub('Cdec', 'Dec', fd$pop)
fd$pop <- gsub('Ceja', 'Eja', fd$pop)
fd$pop <- gsub('Cfus', 'Fus', fd$pop)
fd$pop <- gsub('Cgui', 'Gui', fd$pop)
fd$pop <- gsub('Cmam', 'Mam', fd$pop)

##### GET OVERLAP #####
blocks.ovl <- overlap.wrap(focalblocks = blocks, refblocks = blocks, min.overlap = 0.25)
blocks <- merge(blocks, blocks.ovl, by = c('triplet.short', 'ID'))

blocks$popA <- sapply(strsplit(blocks$triplet, split = '\\.'), '[', 1)
blocks$popB <- sapply(strsplit(blocks$triplet, split = '\\.'), '[', 2)
blocks$popC <- sapply(strsplit(blocks$triplet, split = '\\.'), '[', 3)
blocks$ov.fd.un <- as.integer(blocks$ov.fd.un)

blocks$St <- round(blocks$start / 1000000, 3)
blocks$En <- round(blocks$end / 1000000, 3)
blocks <- arrange(blocks, scaffold, start)


##### ADD UNIQUE/SHARED STATUS #####

## Crater lake pops 3rd, Cmam 2nd:
crater3.cmam <- blocks[blocks$popB == 'Cmam', ]
crater3.cmam$unq <- NA
crater3.cmam$unq[crater3.cmam$triplet == 'Cgui.Cmam.Cdec' & crater3.cmam$UAE == 0 & crater3.cmam$UAF == 0] <- 'unq'
crater3.cmam$unq[crater3.cmam$triplet == 'Cgui.Cmam.Ceja' & crater3.cmam$UAD == 0 & crater3.cmam$UAF == 0] <- 'unq'
crater3.cmam$unq[crater3.cmam$triplet == 'Cgui.Cmam.Cfus' & crater3.cmam$UAD == 0 & crater3.cmam$UAE == 0] <- 'unq'
crater3.cmam$unq[crater3.cmam$UAE >= 1 | crater3.cmam$UAD >= 1 | crater3.cmam$UAF >= 1] <- 'shared1' # Block shared with at least one other crater lake sp
crater3.cmam$unq[crater3.cmam$triplet == 'Cgui.Cmam.Cdec' & crater3.cmam$UAE >= 1 & crater3.cmam$UAF >= 1] <- 'shared2' # Block shared with both other crater lake sp
crater3.cmam$unq[crater3.cmam$triplet == 'Cgui.Cmam.Ceja' & crater3.cmam$UAD >= 1 & crater3.cmam$UAF >= 1] <- 'shared2'
crater3.cmam$unq[crater3.cmam$triplet == 'Cgui.Cmam.Cfus' & crater3.cmam$UAD >= 1 & crater3.cmam$UAE >= 1] <- 'shared2'
crater3.cmam$consistent <- 2

## Crater lake pops 3rd, Cgui 2nd:
crater3.cgui <- blocks[blocks$popB == 'Cgui', ]
crater3.cgui$unq <- NA
crater3.cgui$unq[crater3.cgui$triplet == 'Cmam.Cgui.Cdec' & crater3.cgui$AUE == 0 & crater3.cgui$AUF == 0] <- 'unq'
crater3.cgui$unq[crater3.cgui$triplet == 'Cmam.Cgui.Ceja' & crater3.cgui$AUD == 0 & crater3.cgui$AUF == 0] <- 'unq'
crater3.cgui$unq[crater3.cgui$triplet == 'Cmam.Cgui.Cfus' & crater3.cgui$AUD == 0 & crater3.cgui$AUE == 0] <- 'unq'
crater3.cgui$unq[crater3.cgui$AUE >= 1 | crater3.cgui$AUD >= 1 | crater3.cgui$AUF >= 1] <- 'shared1' # Block shared with at least one other crater lake sp
crater3.cgui$unq[crater3.cgui$triplet == 'Cmam.Cgui.Cdec' & crater3.cgui$AUE >= 1 & crater3.cgui$AUF >= 1] <- 'shared2' # Block shared with both other crater lake sp
crater3.cgui$unq[crater3.cgui$triplet == 'Cmam.Cgui.Ceja' & crater3.cgui$AUD >= 1 & crater3.cgui$AUF >= 1] <- 'shared2'
crater3.cgui$unq[crater3.cgui$triplet == 'Cmam.Cgui.Cfus' & crater3.cgui$AUD >= 1 & crater3.cgui$AUE >= 1] <- 'shared2'
crater3.cgui$consistent <- 2

## Cmam 3rd:
cmam3 <- blocks[blocks$popC == 'Cmam', ]
cmam3$unq <- 'shared1'
cmam3$unq[cmam3$popB == 'Cdec' & cmam3$FEA == 0 & cmam3$DEA == 0  & cmam3$EFA == 0 & cmam3$DFA == 0] <- 'unq'
cmam3$unq[cmam3$popB == 'Ceja' & cmam3$EDA == 0 & cmam3$FDA == 0  & cmam3$EFA == 0 & cmam3$DFA == 0] <- 'unq'
cmam3$unq[cmam3$popB == 'Cfus' & cmam3$FEA == 0 & cmam3$DEA == 0  & cmam3$EDA == 0 & cmam3$FDA == 0] <- 'unq'

cmam3$consistent <- 1
cmam3$consistent[cmam3$popB == 'Cdec' & cmam3$UAD == 0] <- 0
cmam3$consistent[cmam3$popB == 'Ceja' & cmam3$UAE == 0] <- 0
cmam3$consistent[cmam3$popB == 'Cfus' & cmam3$UAF == 0] <- 0
table(cmam3$consistent)

## Cgui 3rd:
cgui3 <- blocks[blocks$popC == 'Cgui', ]
cgui3$unq <- 'shared1'
cgui3$unq[cgui3$popB == 'Cdec' & cgui3$FEU == 0 & cgui3$DEU == 0  & cgui3$EFU == 0 & cgui3$DFU == 0] <- 'unq'
cgui3$unq[cgui3$popB == 'Ceja' & cgui3$EDU == 0 & cgui3$FDU == 0  & cgui3$EFU == 0 & cgui3$DFU == 0] <- 'unq'
cgui3$unq[cgui3$popB == 'Cfus' & cgui3$FEU == 0 & cgui3$DEU == 0  & cgui3$EDU == 0 & cgui3$FDU == 0] <- 'unq'

cgui3$consistent <- 1
cgui3$consistent[cgui3$popB == 'Cdec' & cgui3$AUD == 0] <- 0
cgui3$consistent[cgui3$popB == 'Ceja' & cgui3$AUE == 0] <- 0
cgui3$consistent[cgui3$popB == 'Cfus' & cgui3$AUF == 0] <- 0
table(cgui3$consistent)

blocks <- rbind(crater3.cmam, crater3.cgui, cgui3, cmam3) %>% arrange(scaffold, start)


##### WRITE FILES #####
cat('Writing df:', fd.file.out, '\n')
write.table(fd, fd.file.out, sep = '\t', quote = FALSE, row.names = FALSE)

cat('Writing df:', fdblocks.file.out, '\n')
write.table(blocks, fdblocks.file.out, sep = '\t', quote = FALSE, row.names = FALSE)