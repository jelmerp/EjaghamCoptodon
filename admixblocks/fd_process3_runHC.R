##### SET-UP #####
args <- commandArgs(trailingOnly = TRUE)
file.id <- args[1]
my.block.p <- as.numeric(args[2])
my.mutrate <- 7.5e-9

## Files to read:
fdblocks.file <- 'analyses_output/admixblocks/fdblocks.txt'
varposfile <- paste0('seqdata/fasta/', file.id, '.', scaffold, '.varpos')
fasta.file <- paste0('seqdata/fasta/', file.id, '.', scaffold, '.fasta')
poptriplets.file <- 'analyses_input/admixblocks/poptrios_fd.txt'

## Files to write:
blocks.hc.file <- 'analyses_output/admixblocks/fdblocks.HCtest_intermed.txt'

## Load libraries:
library(HybridCheck, lib.loc = '/nas02/home/j/e/jelmerp/R/x86_64-pc-linux-gnu-library/3.3')
library(data.table)
library(plyr)
library(tidyr)
library(dplyr)

## Define populations:
populations <- list(
  Cdec = c("Cdec088", "Cdec328"),
  Ceja = c("Ceja262", "Ceja408"),
  Cfus = c("Cfus085", "Cfus350", "Cfus503"),
  Cmam = c("SgalMA1"),
  Cgui = c("TguiNG2", "TguiNG5"),
  Sgui = c('TguiMA1', 'TguiMA2', 'TguiMA4')
)
inds <- c('Cdec088', 'Cdec328', 'Ceja262', 'Ceja408', 'Cfus085', 'Cfus350', 'Cfus503', 'TguiNG2', 'TguiNG5', 'SgalMA1')
popnames.long <- c('Cdec', 'Cdec', 'Ceja', 'Ceja', 'Cfus', 'Cfus', 'Cfus', 'Cgui', 'Cgui', 'Cmam')
popnames.short <- c('D', 'D', 'E', 'E', 'F', 'F', 'F', 'U', 'U', 'A')

## Triplets:
triplets <- read.table(triplets.file, header = TRUE, as.is = TRUE)
triplets.long <- paste0(triplets$popA, '.', triplets$popB, '.', triplets$popC)
triplets.short <- gsub('\\.', '', triplets.long)
triplets.short <- gsub('Cfus', 'F', triplets.short)
triplets.short <- gsub('Cdec', 'D', triplets.short)
triplets.short <- gsub('Ceja', 'E', triplets.short)
triplets.short <- gsub('Cmam', 'A', triplets.short)
triplets.short <- gsub('Cgui', 'U', triplets.short)


##### FUNCTIONS #####
abs2rel.pos <- function(scaffold, start.abs, end.abs) {
  varpos <- as.data.frame(fread(varposfile))$V2

  start.rel <- which(varpos > start.abs)[1]
  end.rel <- which(varpos < end.abs)[length(which(varpos < end.abs))]

  return(c(start.rel, end.rel))
}

rel2abs.pos <- function(hcblocks) {
  scaffold <- hcblocks$scaffold[1]
  varpos <- as.data.frame(fread(varposfile))$V2

  hcblocks$start <- varpos[hcblocks$start]
  hcblocks$end <- varpos[hcblocks$end]

  return(hcblocks)
}

testblock <- function(fdblock, mutrate, block.p) {
  scaffold <- fdblock$scaffold

  if(!file.exists(fasta.file)) {

    writeLines(text = scaffold, con = scaffold)
    cat("...........FASTA NOT FOUND\n")

  } else {

    inds.A <- populations[[fdblock$recip]]
    inds.B <- populations[[fdblock$donor]]
    cat('............inds.A:', inds.A, 'inds.B:', inds.B, '\n')

    if(any(grepl('SgalMA1|TguiNG', inds.A))) {
      donor <- inds.A
      recip <- inds.B
    } else {
      donor <- inds.B
      recip <- inds.A
    }
    comb1 <- c(donor[1], recip[1])
    comb2 <- c(donor[1], recip[2])

    cat('............fd block:', fdblock$ID, 'comb1:', comb1, 'comb2:', comb2, '\n')

    ## Get relative positions (sequence with variable sites only):
    startpos.abs <- fdblock$start
    endpos.abs <- fdblock$end + 50000
    pos.rel <- abs2rel.pos(scaffold = scaffold, start.abs = startpos.abs, end.abs = endpos.abs)
    startpos.rel <- pos.rel[1]
    endpos.rel <- pos.rel[2]

    ## Get length difference between total and varsites-only seq, for adjustment of block ages:
    length.rel <- endpos.rel - startpos.rel
    length.abs <- endpos.abs - startpos.abs
    divide.by <- length.abs / length.rel
    cat('Difference between abs and rel length:', divide.by, '\n')

    ## Run Hybridcheck:
    anl <- HC$new(fasta.file)
    anl$setPopulations(populations)
    scaffold.length <- anl$DNA$getFullLength()
    anl$setParameters("BlockDating", MutationRate = my.mutrate, PValue = my.block.p, BonfCorrection = TRUE)
    anl$addUserBlock(comb1, firstbp = startpos.rel, lastbp = endpos.rel)
    anl$addUserBlock(comb2, firstbp = startpos.rel, lastbp = endpos.rel)
    anl$tabulateUserBlocks()
    anl$dateUserBlocks()
    hcb <- anl$tabulateUserBlocks()

    ## Post-processing:
    if(nrow(hcb) >= 1) {
      cat('............', nrow(hcb), 'significant block(s) found\n')

      hcb <- dplyr::rename(hcb, duo = Sequence_Pair, length = ApproxBpLength,
                           start = FirstBP, end = LastBP,
                           nrSNP = SNPs, nrSNPcor = CorrectedSNPs, p.hc = P_Value,
                           age.max = fiveAge, age.mean = fiftyAge, age.min = ninetyFiveAge)

      hcb <- mutate(hcb, age.min = round(age.min / divide.by),
                    age.mean = round(age.mean / divide.by),
                    age.max = round(age.max / divide.by))
      hcb$scaffold <- fdblock$scaffold
      hcb <- rel2abs.pos(hcb)

      hcb <- hcb %>% separate(duo, into = c("f1", "f2"), sep = ":")
      hcb$popA <- mapvalues(hcb$f1, from = inds, to = popnames.short, warn_missing = FALSE)
      hcb$popB <- mapvalues(hcb$f2, from = inds, to = popnames.short, warn_missing = FALSE)
      hcb$duo.hc <- paste0(hcb$popA, hcb$popB)

      hcb$ID <- fdblock$ID
      hcb <- hcb %>% select(ID, start, end, duo.hc, p.hc, age.min, age.mean, age.max)

      return(hcb)
    } else {
      cat("............No significant blocks found.\n")
    }
  }
}

##### APPLY #####

## Read fd blocks:
fdblocks <- read.table(fdblocks.file, header = TRUE, as.is = TRUE)
fdblocks <- dplyr::arrange(fdblocks, desc(fd.sum))

## Run Hybridcheck:
lastblock <- nrow(fdblocks)
cat('Number of blocks:', lastblock, '\n')
fdblocks.hc <- lapply(1:lastblock, FUN = function(x) testblock(fdblocks[x, ], block.p = my.block.p, mutrate = my.mutrate))
fdblocks.hc <- do.call(rbind, fdblocks.hc)

##### WRITE TABLES #####
write.table(fdblocks.hc, blocks.file, sep = '\t', quote = FALSE, row.names = FALSE)
