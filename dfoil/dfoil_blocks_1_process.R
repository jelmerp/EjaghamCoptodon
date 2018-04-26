##### SET-UP #####

## Files to read:
dfoil.output.files <- list.files('analyses_output/dfoil/byWindow', full.names = TRUE)
fdblocks.hc.file <- 'analyses_output/admixblocks/fdblocks.HCtest.txt'

## Files to write:
file.dfoil.all <- 'analyses_output/dfoil/dfoil.all.txt'
dd.filename <- 'analyses/dfoil/intermediateFiles/fdblocks.dfoil.overlap.txt'

## Libraries:
library(reshape2)
library(tidyverse)
library(GenomicRanges)
library(plyr)

## Variables:
pops <- c('noCdec', 'noCeja', 'noCfus')
migs <- c('A2D', 'A2E', 'A2F', 'A2DE', 'A2DEF', 'U2D', 'U2E', 'U2F', 'U2DE', 'U2DEF',
          'D2A', 'E2A', 'F2A', 'D2U', 'E2U', 'F2U')
curpatterns <- c('A2D', 'A2E', 'A2F', 'U2D', 'U2E', 'U2F', 'D2A', 'E2A', 'F2A', 'D2U', 'E2U', 'F2U')


##### FUNCTIONS #####
overlap.get <- function(rownr, focalblocks, refblocks, min.overlap = 0.5) {

  foc <- focalblocks[rownr, ]
  ref <- refblocks[(rownr + 1):nrow(refblocks), ]

  foc.gr <- as(paste0(foc$scaffold, ':', foc$start, '-', foc$end), "GRanges")
  ref.gr <- as(paste0(ref$scaffold, ':', ref$start, '-', ref$end), "GRanges")

  if(min.overlap < 1) min.overlap <- (foc$end - foc$start) * min.overlap

  hits <- findOverlaps(foc.gr, ref.gr, minoverlap = min.overlap)
  hits <- ref[hits@to, ]

  hits.net <- length(unique(hits[which(hits$intr2 != foc$intr2), ]$hits$intr2))
  migpatterns <- unique(c(foc$intr2, hits$intr2))

  popmatches <- as.integer(gsub(TRUE, 1, gsub(FALSE, 0, migs %in% migpatterns)))
  match.IDs <- hits$ID

  cat('Row number:', rownr, '   Number of hits:', nrow(hits), '\n')

  cat('Focal ID:', foc$ID, '         First match ID:', match.IDs[1], '\n')
  return(list(c(foc$ID, foc$popconfig, foc$intr, foc$intr2, nrow(hits), hits.net, popmatches, paste(migpatterns, collapse = '_')),
              match.IDs))
}

overlap.wrap <- function(focalblocks, refblocks, min.overlap = 0.25, lastblock = NULL) {

  if(is.null(lastblock)) lastblock <- nrow(focalblocks) - 1 # lastblock <- 25
  overlaps.list <- lapply(1:lastblock, overlap.get,
                          focalblocks = focalblocks, refblocks = refblocks, min.overlap = min.overlap)

  match.IDs <- unique(unlist(lapply(overlaps.list, '[', 2)))

  overlaps <- lapply(overlaps.list, '[', 1)
  overlaps <- rbind.fill(lapply(overlaps, function(x){as.data.frame(t(unlist(x)), stringsAsFactors = FALSE)}))
  colnames(overlaps) <- c('ID', 'popconfig', 'intr', 'intr2', 'ov.dfoil', 'ov2.dfoil', migs, 'migpatterns')

  overlaps <- overlaps[! overlaps$ID %in% match.IDs, ]

  return(overlaps)
}

##### PREP  #####
dfoil.all <- do.call(rbind, lapply(dfoil.output.files, read.delim, header = TRUE, as.is = TRUE))

dfoil.all$X.chrom <- gsub('\\.dfoil', '', gsub('fasta\\.', '', gsub(".*(N[T|C]_.*)", '\\1', dfoil.all$X.chrom)))
dfoil.all$popconfig <- gsub('.*(noC.*)', '\\1', dfoil.all$X.chrom)
dfoil.all$X.chrom <- gsub('(.*)\\.noC.*', '\\1', dfoil.all$X.chrom)
dfoil.all$popA.dfoil <- ifelse(dfoil.all$popconfig == 'noCdec', 'Ceja', 'Cdec')
dfoil.all$popB.dfoil <- ifelse(dfoil.all$popconfig == 'noCfus', 'Ceja', 'Cfus')

dfoil.all <- dfoil.all %>% dplyr::rename(ID = X.chrom, intr.nrs = introgression,
                                 DFO.l = DFO_left, DFO.r = DFO_right, DFO.t = DFO_total, DFO.stat = DFO_stat, DFO.p = DFO_Pvalue,
                                 DIL.l = DIL_left, DIL.r = DIL_right, DIL.t = DIL_total, DIL.stat = DIL_stat, DIL.p = DIL_Pvalue,
                                 DFI.l = DFI_left, DFI.r = DFI_right, DFI.t = DFI_total, DFI.stat = DFI_stat, DFI.p = DFI_Pvalue,
                                 DOL.l = DOL_left, DOL.r = DOL_right, DOL.t = DOL_total, DOL.stat = DOL_stat, DOL.p = DOL_Pvalue)

dfoil.all <- dfoil.all %>% mutate(DFO.stat = round(DFO.stat, 2), DFO.p = round(DFO.p, 4), DFO.diff = abs(DFO.l - DFO.r),
                          DIL.stat = round(DIL.stat, 2), DIL.p = round(DIL.p, 4), DIL.diff = abs(DIL.l - DIL.r),
                          DFI.stat = round(DFI.stat, 2), DFI.p = round(DFI.p, 4), DFI.diff = abs(DFI.l - DFI.r),
                          DOL.stat = round(DOL.stat, 2), DOL.p = round(DOL.p, 4), DOL.diff = abs(DOL.l - DOL.r))

dfoil.all$DFO <- ifelse(dfoil.all$DFO.diff < 10, '0', ifelse(dfoil.all$DFO.stat < 0, '-', '+'))
dfoil.all$DIL <- ifelse(dfoil.all$DIL.diff < 10, '0', ifelse(dfoil.all$DIL.stat < 0, '-', '+'))
dfoil.all$DFI <- ifelse(dfoil.all$DFI.diff < 10, '0', ifelse(dfoil.all$DFI.stat < 0, '-', '+'))
dfoil.all$DOL <- ifelse(dfoil.all$DOL.diff < 10, '0', ifelse(dfoil.all$DOL.stat < 0, '-', '+'))

dfoil.all$intr2.nrs <- 0
dfoil.all$intr2.nrs[dfoil.all$DFO == '+' & dfoil.all$DIL == '+' & dfoil.all$DFI == '+' & dfoil.all$DOL == 0] <- 13
dfoil.all$intr2.nrs[dfoil.all$DFO == '+' & dfoil.all$DIL == 0 & dfoil.all$DFI == '+' & dfoil.all$DOL == '+'] <- 31
dfoil.all$intr2.nrs[dfoil.all$DFO == '-' & dfoil.all$DIL == '-' & dfoil.all$DFI == 0 & dfoil.all$DOL == '+'] <- 14
dfoil.all$intr2.nrs[dfoil.all$DFO == '-' & dfoil.all$DIL == 0 & dfoil.all$DFI == '+' & dfoil.all$DOL == '+'] <- 41
dfoil.all$intr2.nrs[dfoil.all$DFO == '+' & dfoil.all$DIL == '+' & dfoil.all$DFI == '-' & dfoil.all$DOL == 0] <- 23
dfoil.all$intr2.nrs[dfoil.all$DFO == 0 & dfoil.all$DIL == '+' & dfoil.all$DFI == '-' & dfoil.all$DOL == '-'] <- 32
dfoil.all$intr2.nrs[dfoil.all$DFO == '-' & dfoil.all$DIL == '-' & dfoil.all$DFI == 0 & dfoil.all$DOL == '-'] <- 24
dfoil.all$intr2.nrs[dfoil.all$DFO == 0 & dfoil.all$DIL == '-' & dfoil.all$DFI == '-' & dfoil.all$DOL == '-'] <- 42
dfoil.all$intr2.nrs[dfoil.all$DFO == '+' & dfoil.all$DIL == '+' & dfoil.all$DFI == 0 & dfoil.all$DOL == 0] <- 123
dfoil.all$intr2.nrs[dfoil.all$DFO == '-' & dfoil.all$DIL == '-' & dfoil.all$DFI == 0 & dfoil.all$DOL == 0] <- 124

dfoil.all$intr <- 0
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr.nrs == 13] <- 'D2U'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr.nrs == 31] <- 'U2D'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr.nrs == 14] <- 'D2A'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr.nrs == 41] <- 'A2D'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr.nrs == 23] <- 'E2U'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr.nrs == 32] <- 'U2E'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr.nrs == 24] <- 'E2A'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr.nrs == 42] <- 'A2E'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr.nrs == 123] <- 'U2DE'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr.nrs == 124] <- 'A2DE'

dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 13] <- 'D2U'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 31] <- 'U2D'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 14] <- 'D2A'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 41] <- 'A2D'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 23] <- 'F2U'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 32] <- 'U2F'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 24] <- 'F2A'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 42] <- 'A2F'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 123] <- 'U2DEF'
dfoil.all$intr[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 124] <- 'A2DEF'

dfoil.all$intr[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 13] <- 'E2U'
dfoil.all$intr[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 31] <- 'U2E'
dfoil.all$intr[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 14] <- 'E2A'
dfoil.all$intr[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 41] <- 'A2E'
dfoil.all$intr[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 23] <- 'F2U'
dfoil.all$intr[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 32] <- 'U2F'
dfoil.all$intr[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 24] <- 'F2A'
dfoil.all$intr[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 42] <- 'A2F'
dfoil.all$intr[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 123] <- 'U2DEF'
dfoil.all$intr[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr.nrs == 124] <- 'A2DEF'

dfoil.all$intr2 <- 0
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr2.nrs == 13] <- 'D2U'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr2.nrs == 31] <- 'U2D'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr2.nrs == 14] <- 'D2A'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr2.nrs == 41] <- 'A2D'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr2.nrs == 23] <- 'E2U'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr2.nrs == 32] <- 'U2E'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr2.nrs == 24] <- 'E2A'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr2.nrs == 42] <- 'A2E'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr2.nrs == 123] <- 'U2DE'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Ceja' & dfoil.all$intr2.nrs == 124] <- 'A2DE'

dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 13] <- 'D2U'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 31] <- 'U2D'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 14] <- 'D2A'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 41] <- 'A2D'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 23] <- 'F2U'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 32] <- 'U2F'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 24] <- 'F2A'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 42] <- 'A2F'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 123] <- 'U2DEF'
dfoil.all$intr2[dfoil.all$popA == 'Cdec' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 124] <- 'A2DEF'

dfoil.all$intr2[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 13] <- 'E2U'
dfoil.all$intr2[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 31] <- 'U2E'
dfoil.all$intr2[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 14] <- 'E2A'
dfoil.all$intr2[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 41] <- 'A2E'
dfoil.all$intr2[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 23] <- 'F2U'
dfoil.all$intr2[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 32] <- 'U2F'
dfoil.all$intr2[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 24] <- 'F2A'
dfoil.all$intr2[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 42] <- 'A2F'
dfoil.all$intr2[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 123] <- 'U2DEF'
dfoil.all$intr2[dfoil.all$popA == 'Ceja' & dfoil.all$popB == 'Cfus' & dfoil.all$intr2.nrs == 124] <- 'A2DEF'

dfoil.all$intr2[dfoil.all$intr != 'none' & dfoil.all$intr != 'na' & dfoil.all$intr != 0] <- dfoil.all$intr[dfoil.all$intr != 'none' & dfoil.all$intr != 'na' & dfoil.all$intr != 0]

dfoil.all$scaffold <- paste0(sapply(strsplit(dfoil.all$ID, split = '\\.'), '[', 1), '.1')
dfoil.all$start <- as.integer(sapply(strsplit(dfoil.all$ID, split = '\\.'), '[', 3))
dfoil.all$end <- as.integer(sapply(strsplit(dfoil.all$ID, split = '\\.'), '[', 4))
dfoil.all$ID <- paste0(dfoil.all$ID, '_', dfoil.all$popconfig)

dfoil.all <- dfoil.all %>% filter(intr.nrs %in% c(123, 124, 14, 41, 13, 31, 24, 42, 23, 32) | intr2.nrs %in% c(123, 124, 14, 41, 13, 31, 24, 42, 23, 32))

dfoil.all <- dfoil.all %>% select(ID, scaffold, start, end, popconfig, popA.dfoil, popB.dfoil,
                                  total, dtotal,
                                  DFO.l, DFO.r, DFO.t, DFO.stat, DFO.diff, DFO.p,
                                  DIL.l, DIL.r, DIL.t, DIL.stat, DIL.diff, DIL.p,
                                  DFI.l, DFI.r, DFI.t, DFI.stat, DFI.diff, DFI.p,
                                  DOL.l, DOL.r, DOL.t, DOL.stat, DOL.diff, DOL.p,
                                  DFO, DIL, DFI, DOL, intr, intr2)

write.table(dfoil.all, file.dfoil.all, sep = '\t', quote = FALSE, row.names = FALSE)

##### OVERLAP #####
do <- overlap.wrap(focalblocks = dfoil.all, refblocks = dfoil.all, min.overlap = 0.25)
do$ID2 <- paste0(sapply(strsplit(do$ID, split = '_'), '[', 1), '_',
                 sapply(strsplit(do$ID, split = '_'), '[', 2))
do <- dplyr::rename(do, ID.dfoil = ID)
backup <- do


##### COMBINE WITH FD.HC RESULTS #####
## Get rid of duplicate blocks in fd.hc (same coordinates, different triplet):
fd.hc <- read.table(fdblocks.hc.file, header = TRUE, as.is = TRUE)
fd.hc <- dplyr::rename(fd.hc, tri = triplet.short)
fd.hc$ID2 <- paste0(fd.hc$scaffold, '.', fd.hc$start, '.', fd.hc$end)
fd.hc <- dplyr::arrange(fd.hc, scaffold, start, end, desc(p.hc))
fd.hc$tri2 <- NA
fd.hc$dup <- 0
for(i in 1:(nrow(fd.hc) - 1)) {
  if(fd.hc[i, ]$ID2 == fd.hc[i + 1, ]$ID2) {
    fd.hc[i, ]$dup <- 1
    if(is.na(fd.hc[i, ]$tri2)) fd.hc[i + 1, ]$tri2 <- paste0(fd.hc[i, ]$tri, '_', fd.hc[i + 1, ]$tri)
    if(!is.na(fd.hc[i, ]$tri2)) fd.hc[i + 1, ]$tri2 <- paste0(fd.hc[i, ]$tri2, '_', fd.hc[i + 1, ]$tri)
  }
}
fd.hc <- dplyr::filter(fd.hc, dup == 0) %>% dplyr::select(-dup)

## Merge:
dd <- merge(fd.hc, do, by = 'ID2')


##### GET MIGRATION PATTERN #####
dd$pattern <- 0
dd$ancDir <- 0
dd$pattern[dd$ov2.dfoil == 0] <- dd$intr2[do$ov2.dfoil == 0]
dd$pattern[dd$A2DE == 1] <- 'A2DE'
dd$pattern[dd$A2DEF == 1] <- 'A2DEF'
dd$pattern[dd$U2DE == 1] <- 'U2DE'
dd$pattern[dd$U2DEF == 1] <- 'U2DEF'
dd$pattern[dd$U2D == 1 & dd$U2F == 1] <- 'U2DEF'

dfoil.popA <- sapply(strsplit(dd$pattern, split = '2'), '[', 1)
dfoil.popB <- sapply(strsplit(dd$pattern, split = '2'), '[', 2)

match1 <- sapply(1:nrow(dd), FUN = function(i) grepl(pattern = dfoil.popA[i], x = dd$duo.hc[i]))
match2 <- sapply(1:nrow(dd), FUN = function(i) any(unlist(strsplit(dfoil.popB[2], split = '')) %in%
                                                     unlist(strsplit(dd$duo.hc[2], split = ''))))
dd$dfoil.consistent <- ifelse(match1 + match2 < 2, 0, 1)

## Manually modify some inconsistent inferences:
#filter(dd, dfoil.consistent != 0, consistent == 1)
dd$pattern[dd$ID == 'NC_022199.1_25900001_EFU'] <- 'U2F'; dd$dfoil.consistent[dd$ID == 'NC_022199.1_25900001_EFU'] <- 1

#filter(dd, intr != 0, intr != pattern)
dd$pattern[dd$ID == 'NT_167899.1_130001_FEU'] <- dd$intr[dd$ID == 'NT_167899.1_130001_FEU']
dd$pattern[dd$ID == 'NT_167671.1_400001_UAE'] <- dd$intr[dd$ID == 'NT_167671.1_400001_UAE']
dd$pattern[dd$ID == 'NC_022200.1_18710001_UAD'] <- dd$intr[dd$ID == 'NC_022200.1_18710001_UAD']

## Filter:
dd <- filter(dd, dfoil.consistent != 0)

## Specify migration patterns:
dd$ancDir[dd$pattern == 'A2DEF' & (dd$A2D == 1 | dd$A2E == 1 | dd$A2F == 1)] <- 'into'
dd$ancDir[dd$pattern == 'A2DE' & (dd$A2D == 1 | dd$A2E == 1 | dd$A2F == 1)] <- 'into'
dd$ancDir[dd$pattern == 'A2DEF' & (dd$D2A == 1 | dd$E2A == 1 | dd$F2A == 1)] <- 'from'
dd$ancDir[dd$pattern == 'A2DE' & (dd$D2A == 1 | dd$E2A == 1 | dd$F2A == 1)] <- 'from'
dd$ancDir[dd$pattern == 'U2DEF' & (dd$U2D == 1 | dd$U2E == 1 | dd$U2F == 1)] <- 'into'
dd$ancDir[dd$pattern == 'U2DE' & (dd$U2D == 1 | dd$U2E == 1 | dd$U2F == 1)] <- 'into'
dd$ancDir[dd$pattern == 'U2DEF' & (dd$D2U == 1 | dd$E2U == 1 | dd$F2U == 1)] <- 'from'
dd$ancDir[dd$pattern == 'U2DE' & (dd$D2U == 1 | dd$E2U == 1 | dd$F2U == 1)] <- 'from'
dd$ancDir[dd$U2D == 1 & dd$U2F == 1] <- 'into'

dd$AncCur <- ifelse(grepl('DE|DEF', dd$pattern), 'anc', ifelse(dd$pattern %in% curpatterns, 'cur', 0))
dd$migDir <- ifelse(dd$pattern %in% c('A2D', 'A2E', 'A2F', 'U2D', 'U2E', 'U2F'), 'into',
                          ifelse(dd$pattern %in% c('D2A', 'E2A', 'F2A', 'D2U', 'E2U', 'F2U'), 'from',
                          ifelse(dd$ancDir == 'into', 'into', ifelse(dd$ancDir == 'from', 'from', 0))))

dd$lakepop <- sapply(strsplit(dd$pattern, split = 2), '[', 2)
dd$lakepop[dd$lakepop %in% c('A', 'U')] <- sapply(strsplit(dd$pattern[dd$lakepop %in% c('A', 'U')], split = 2), '[', 1)


##### MERGE OVERLAPS BACK WITH DFOIL #####
dd <- dplyr::select(dd, -ID2, -popconfig, -intr2, -ov.dfoil, -ov2.dfoil, -recip, -donor, -triplet, -p.min,
                    -A2D, -A2E, -A2F, -A2DE, -A2DEF, -U2D, -U2E, -U2F, -U2DE, -U2DEF, -D2A, -E2A, -F2A, -D2U,
                    -E2U, -F2U, -FDA, -FDU, -FEA, -FEU, -DEA, -DEU, -DFA, -DFU, -EFA, -EFU, -EDA, -EDU, -AUD,
                    -AUE, -AUF, -UAD , -UAE, -UAF, -UDA, -UEA, -UFA)
dd <- dplyr::rename(dd, consistent.fd = consistent, consistent.dfoil = dfoil.consistent,
                    migpatterns.dfoil = migpatterns, migpattern = pattern)
dd <- dplyr::select(dd, ID, ID.dfoil, tri, tri2, scaffold, nwin, start, end, fd.sum, fd.mean, fd.max, p.mean,
                    fst.AC, fst.BC, dxy.AC, dxy.BC, pi.A, pi.B, pi.C, tajD.A, tajD.B, tajD.C, ov.fd, ov.fd.un,
                    popA, popB, popC, St, En, duo.hc, p.hc, age.min, age.mean, age.max, ancDir, AncCur,
                    migDir, riverpop, lakepop, singleHC, consistent.fd, consistent.dfoil, unq,
                    migpatterns.dfoil, migpattern, intr)
dfoil.all <- select(dfoil.all, -intr, -intr2, -scaffold, -start, -end)
dd <- merge(dd, dfoil.all, by.x = 'ID.dfoil', by.y = 'ID')


##### FINAL EDITS #####
dd <- dplyr::filter(dd, ID != 'NC_022204.1_4860001_UAD')
dd$lakepop <- as.character(dd$lakepop)
dd$lakepoptype <- mapvalues(dd$lakepop,
                            from = c('D', 'E', 'F', 'DE', 'DEF'),
                            to = c(rep('extant species', 3), 'ancestor DE', 'ancestor DEF'))
dd$lakepoptype <- factor(dd$lakepoptype, levels = c('extant species', 'ancestor DE', 'ancestor DEF'))


##### WRITE TABLES #####
write.table(dd, file = dd.filename, sep = '\t', quote = FALSE, row.names = FALSE)
