##### SET-UP #####

## Files to read:
pgfile <- 'analyses_output/sumstats/pg.collected_EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01.win50k.step5k.txt'
tajD.file <- 'analyses_output/sumstats/tajD.combined.EjaC.Cgal.Cgui.DP5.GQ20.MAXMISS0.5.MAF0.01.win10000.txt'

popIDs.file <- 'analyses_input/admixblocks/popIDs_fd.txt'
poptriplets.file <- 'analyses_input/admixblocks/poptrios_fd.txt'

## Files to write:
fd.file <- 'analyses_output/admixblocks/fd_intermediate.txt'
blocks.file <- 'analyses_output/admixblocks/fdblocks_intermediate.txt'

## Libraries:
library(data.table)
library(plyr)
library(dplyr)

## Variables:
popnames.long <- c('Cdec', 'Ceja', 'Cfus', 'Cgui', 'Cmam')
popnames.short <- c('D', 'E', 'F', 'U', 'A')

## Read data files:
pg <- as.data.frame(fread(pgfile, stringsAsFactors = FALSE))
pg$pop <- paste0(pg$popA, '.', pg$popB)
tajD <- as.data.frame(fread(tajD.file))

## Triplets:
triplets <- read.table(poptriplets.file, as.is = TRUE)
triplets.long <- paste0(triplets$popA, '.', triplets$popB, '.', triplets$popC)
triplets.short <- gsub('\\.', '', triplets.long)
triplets.short <- gsub('Cfus', 'F', triplets.short)
triplets.short <- gsub('Cdec', 'D', triplets.short)
triplets.short <- gsub('Ceja', 'E', triplets.short)
triplets.short <- gsub('Cmam', 'A', triplets.short)
triplets.short <- gsub('Cgui', 'U', triplets.short)

## Populations:
inds <- read.table(popIDs.file, as.is = TRUE)
get.pop <- function(colnr) substr(sapply(strsplit(inds[, colnr], split = ','), '[', 1), 1, 4)
pops <- data.frame(sapply(1:ncol(inds), get.pop))
colnames(pops) <- c('popA', 'popB', 'popC', 'popD')
pops$popA <- gsub('Sgal', 'Cmam', pops$popA)
pops$popA <- gsub('Tgui', 'Cgui', pops$popA)
pops$popB <- gsub('Sgal', 'Cmam', pops$popB)
pops$popB <- gsub('Tgui', 'Cgui', pops$popB)
pops$popC <- gsub('Sgal', 'Cmam', pops$popC)
pops$popC <- gsub('Tgui', 'Cgui', pops$popC)


##### FUNCTIONS #####
combine.fd <- function(file.id, poplines = 'all') {

  if(poplines[1] == 'all') {
    cat('Running all files... \n')
    fd.and.blocks <- lapply(1:nrow(pops), read.fd, file.id = file.id)
  } else {
    cat('Running files:', poplines, '\n')
    fd.and.blocks <- lapply(poplines[1]:poplines[2], read.fd, file.id = file.id)
  }

  fd <- do.call(rbind, lapply(fd.and.blocks, '[[', 1))
  blocks <- do.call(rbind, lapply(fd.and.blocks, '[[', 2))

  return(list(fd, blocks))
}

read.fd <- function(linenr, file.id = 'EjaC.Dstat.DP5.GQ20.win50000.step5000') {

  filename <- paste0('analyses/admixblocks/fd/output/ABBABABAoutput_',
                     file.id, '.popfileline', linenr, '.csv')

  if(file.exists(filename)) {
    cat('Reading', filename, '\n')
    out <- read.csv(filename, as.is = TRUE)

    out$fd[which(is.na(out$fd))] <- 0
    out$fd[which(is.infinite(out$fd))] <- 0
    out$fd[which(out$fd < 0)] <- 0
    out$fd[which(out$fd > 1)] <- 1

    out$popA <- pops[linenr, 1]
    out$popB <- pops[linenr, 2]
    out$popC <- pops[linenr, 3]
    out$pop <- paste0(out$popA, '.', out$popB, '.', out$popC)
    cat(out$pop[1], '\n')

    out$ABBA <- round(out$ABBA)
    out$BABA <- round(out$BABA)

    out$p.fd <- round(qp.p(out, var = 'fd'), 5)

    if(any(out$p.fd < 0.05)) {
      blocks <- qpBlockStats(blocks = qpBlock(out), out = out)
    } else {
      cnames.blocks <- c('pop', 'scaffold', 'nwin', 'win.first', 'win.last', 'fd.sum', 'fd.mean', 'fd.max', 'p.mean', 'p.min',
                         'fst.AC', 'fst.BC', 'dxy.AC', 'dxy.BC', 'pi.A', 'pi.B', 'pi.C', 'tajD.A', 'tajD.B', 'tajD.C')
      blocks <- data.frame(matrix(nrow = 1, ncol = length(cnames.blocks), dimnames = list(1, cnames.blocks)))
    }

    out <- select(out, popA, popB, popC, scaffold, start, end, mid, sites, sitesUsed, ABBA, BABA, D, fd, p.fd)

  } else {
    cat("Can't find file:", filename, '\n')
    cnames.out <- c('popA', 'popB', 'popC', 'scaffold', 'start', 'end', 'mid',
                    'sites', 'sitesUsed', 'ABBA', 'BABA', 'D', 'fd', 'p.fd')
    cnames.blocks <- c('pop', 'scaffold', 'nwin', 'win.first', 'win.last',
                       'fd.sum', 'fd.mean', 'fd.max', 'p.mean', 'p.min',
                       'fst.AC', 'fst.BC', 'dxy.AC', 'dxy.BC', 'pi.A', 'pi.B', 'pi.C',
                       'tajD.A', 'tajD.B', 'tajD.C')
    out <- data.frame(matrix(nrow = 1, ncol = length(cnames.out), dimnames = list(1, cnames.out)))
    blocks <- data.frame(matrix(nrow = 1, ncol = length(cnames.blocks), dimnames = list(1, cnames.blocks)))
  }

  return(list(out, blocks))
}



qp.p <- function(out, var = 'fd', p.adjust = TRUE) {

  varvals <- out[, var]

  Z <- scale(varvals, center = TRUE, scale = TRUE)
  p <- 2 * pnorm(-abs(Z))

  if(out$popA[1] == 'Cmam' & out$popB[1] == 'Cgui') {
    p.adjust <- FALSE
    p <- p * 10
  }

  if(p.adjust == TRUE) p <- p.adjust(p, method = 'BH')

  if(any(p < 0.05)) cat('Proportion significant:', round(length(which(p < 0.05)) / length(p), 3), '\n')
  if(all(p >= 0.05)) cat('None significant!\n')

  return(p)
}

qpBlock <- function(out) {
  sig <- which(out$p.fd < 0.05)
  blocks <- split(sig, cumsum(c(1, diff(sig) > 5)))
  blocks <- data.frame(cbind(first = sapply(blocks, '[', 1), last = sapply(blocks, tail, 1)))

  sigtest <- sig[1:25]
  blocks <- split(sigtest, cumsum(c(1, diff(sigtest) > 5)))
  blocks <- data.frame(cbind(first = sapply(blocks, '[', 1), last = sapply(blocks, tail, 1)))


  return(blocks)
}

qpBlockStats <- function(blocks = NULL, out = NULL) {

  qpBlockStats.one <- function(blocksrow = NULL, blocks = NULL, out = NULL) {
    cat(blocksrow, '\n')
    block <- blocks[blocksrow, ]

    scaffold <- out$scaffold[block$first]
    pop <- out$pop[block$first]

    nwin <- (block$last - block$first) + 1
    win.first <- out$start[block$first]
    win.last <- out$start[block$last]

    fd.sum <- sum(out[block$first:block$last, ]$fd)
    fd.mean <- mean(out[block$first:block$last, ]$fd)
    fd.max <- max(out[block$first:block$last, ]$fd)
    p.mean <- mean(out[block$first:block$last, ]$p.fd)
    p.min <- min(out[block$first:block$last, ]$p.fd)

    popA <- unlist(strsplit(pop, split = '\\.'))[1]
    popB <- unlist(strsplit(pop, split = '\\.'))[2]
    popC <- unlist(strsplit(pop, split = '\\.'))[3]

    pg.AC.block <- pg[pg$scaffold == scaffold & pg$pop == paste0(popA, '.', popC), ]
    if(nrow(pg.AC.block) > 0) {
      if(any(pg.AC.block$start >= win.first)) pg.AC.blockstart <- which(pg.AC.block$start >= win.first)[1] else pg.AC.blockstart <- NULL
      if(any(pg.AC.block$start <= win.last)) pg.AC.blockend <- tail(which(pg.AC.block$start <= win.last), 1) else pg.AC.blockend <- NULL

      if(!is.null(pg.AC.blockstart) & !is.null(pg.AC.blockend)) {
        fst.AC <- mean(pg.AC.block[pg.AC.blockstart:pg.AC.blockend, ]$fst, na.rm = TRUE)
        dxy.AC <- mean(pg.AC.block[pg.AC.blockstart:pg.AC.blockend, ]$dxy, na.rm = TRUE)
        pi.A <- mean(pg.AC.block[pg.AC.blockstart:pg.AC.blockend, ]$pi_popA, na.rm = TRUE)
      } else {
        fst.AC <- NA
        dxy.AC <- NA
        pi.A <- NA
      }
    } else {
      fst.AC <- NA
      dxy.AC <- NA
      pi.A <- NA
    }

    pg.BC.block <- pg[pg$scaffold == scaffold & pg$pop == paste0(popB, '.', popC), ]
    if(nrow(pg.BC.block) > 0) {
      if(any(pg.BC.block$start >= win.first)) pg.BC.blockstart <- which(pg.BC.block$start >= win.first)[1] else pg.BC.blockstart <- NULL
      if(any(pg.BC.block$start <= win.last)) pg.BC.blockend <- tail(which(pg.BC.block$start <= win.last), 1) else pg.BC.blockend <- NULL

      if(!is.null(pg.BC.blockstart) & !is.null(pg.BC.blockend))  {
        fst.BC <- mean(pg.BC.block[pg.BC.blockstart:pg.BC.blockend, ]$fst, na.rm = TRUE)
        dxy.BC <- mean(pg.BC.block[pg.BC.blockstart:pg.BC.blockend, ]$dxy, na.rm = TRUE)
        pi.B <- mean(pg.BC.block[pg.BC.blockstart:pg.BC.blockend, ]$pi_popA, na.rm = TRUE)
        pi.C <- mean(pg.BC.block[pg.BC.blockstart:pg.BC.blockend, ]$pi_popB, na.rm = TRUE)
      } else {
        fst.BC <- NA
        dxy.BC <- NA
        pi.B <- NA
        pi.C <- NA
      }
    } else {
      fst.BC <- NA
      dxy.BC <- NA
      pi.B <- NA
      pi.C <- NA
    }

    tajD.A.block <- tajD[tajD$scaffold == scaffold & tajD$pop == popA, ]
    if(nrow(tajD.A.block) > 0) {
      if(any(tajD.A.block$start >= win.first)) tajD.A.blockstart <- which(tajD.A.block$start >= win.first)[1] else tajD.A.blockstart <- NULL
      if(any(tajD.A.block$start <= win.last)) tajD.A.blockend <- tail(which(tajD.A.block$start <= win.last), 1) else tajD.A.blockend <- NULL
      if(!is.null(tajD.A.blockstart) & !is.null(tajD.A.blockend)) {
        tajD.A <- mean(tajD.A.block[tajD.A.blockstart:tajD.A.blockend, ]$tajD, na.rm = TRUE)
      } else tajD.A <- NA
    } else tajD.A <- NA

    tajD.B.block <- tajD[tajD$scaffold == scaffold & tajD$pop == popB, ]
    if(nrow(tajD.B.block) > 0) {
      if(any(tajD.B.block$start >= win.first)) tajD.B.blockstart <- which(tajD.B.block$start >= win.first)[1] else tajD.B.blockstart <- NULL
      if(any(tajD.B.block$start <= win.last)) tajD.B.blockend <- tail(which(tajD.B.block$start <= win.last), 1) else tajD.B.blockend <- NULL
      if(!is.null(tajD.B.blockstart) & !is.null(tajD.B.blockend)) {
        tajD.B <- mean(tajD.B.block[tajD.B.blockstart:tajD.B.blockend, ]$tajD, na.rm = TRUE)
      } else tajD.B <- NA
    } else tajD.B <- NA

    tajD.C.block <- tajD[tajD$scaffold == scaffold & tajD$pop == popC, ]
    if(nrow(tajD.C.block) > 0) {
      if(any(tajD.C.block$start >= win.first)) tajD.C.blockstart <- which(tajD.C.block$start >= win.first)[1] else tajD.C.blockstart <- NULL
      if(any(tajD.C.block$start <= win.last)) tajD.C.blockend <- tail(which(tajD.C.block$start <= win.last), 1) else tajD.C.blockend <- NULL
      if(!is.null(tajD.C.blockstart) & !is.null(tajD.C.blockend)) {
        tajD.C <- mean(tajD.C.block[tajD.C.blockstart:tajD.C.blockend, ]$tajD, na.rm = TRUE)
        } else tajD.C <- NA
    } else tajD.C <- NA

    return(c(pop, scaffold, nwin, win.first, win.last, fd.sum, fd.mean, fd.max, p.mean, p.min,
             fst.AC, fst.BC, dxy.AC, dxy.BC, pi.A, pi.B, pi.C, tajD.A, tajD.B, tajD.C))
  }

  blockstats <- do.call(rbind, lapply(1:nrow(blocks), qpBlockStats.one, blocks = blocks, out = out))
  blockstats <- data.frame(blockstats, stringsAsFactors = FALSE)
  colnames(blockstats) <- c('pop', 'scaffold', 'nwin', 'win.first', 'win.last',
                            'fd.sum', 'fd.mean', 'fd.max', 'p.mean', 'p.min',
                            'fst.AC', 'fst.BC', 'dxy.AC', 'dxy.BC', 'pi.A', 'pi.B', 'pi.C',
                            'tajD.A', 'tajD.B', 'tajD.C')

  blockstats$nwin <- as.integer(blockstats$nwin)
  blockstats$win.first <- as.integer(blockstats$win.first)
  blockstats$win.last <- as.integer(blockstats$win.last)
  blockstats$fd.sum <- round(as.numeric(blockstats$fd.sum), 2)
  blockstats$fd.mean <- round(as.numeric(blockstats$fd.mean), 2)
  blockstats$fd.max <- round(as.numeric(blockstats$fd.max), 2)
  blockstats$p.mean <- round(as.numeric(blockstats$p.mean), 5)
  blockstats$p.min <- round(as.numeric(blockstats$p.min), 5)
  blockstats$fst.AC <- round(as.numeric(blockstats$fst.AC), 4)
  blockstats$fst.BC <- round(as.numeric(blockstats$fst.BC), 4)
  blockstats$dxy.AC <- round(as.numeric(blockstats$dxy.AC), 4)
  blockstats$dxy.BC <- round(as.numeric(blockstats$dxy.BC), 4)
  blockstats$pi.A <- round(as.numeric(blockstats$pi.A), 4)
  blockstats$pi.B <- round(as.numeric(blockstats$pi.B), 4)
  blockstats$pi.C <- round(as.numeric(blockstats$pi.C), 4)
  blockstats$tajD.A <- round(as.numeric(blockstats$tajD.A), 2)
  blockstats$tajD.B <- round(as.numeric(blockstats$tajD.B), 2)
  blockstats$tajD.C <- round(as.numeric(blockstats$tajD.C), 2)

  return(blockstats)
}


##### APPLY #####
fd.and.blocks <- combine.fd(file.id = file.id)

fd <- fd.and.blocks[[1]]
fd$pop <- paste0(fd$popA, '.', fd$popB, '.', fd$popC)
write.table(fd, fd.file, sep = '\t', quote = FALSE, row.names = FALSE)

blocks <- fd.and.blocks[[2]]
if(any(is.na(blocks$pop))) blocks <- blocks[-which(is.na(blocks$pop)), ]
write.table(blocks, blocks.file, sep = '\t', quote = FALSE, row.names = FALSE)

