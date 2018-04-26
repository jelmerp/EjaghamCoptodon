#!/usr/bin/env Rscript

## sytem("cd /home/jelmer/Dropbox/sc_fish/cichlids")
## system("Rscript /home/jelmer/Dropbox/sc_fish/cichlids/scripts/gphocs/chooseLoci.R NT_167980.1 $PWD 500)
## system("Rscript /home/jelmer/Dropbox/sc_fish/cichlids/scripts/gphocs/createSeqs.R EjaC.wOut NT_167980.1 $PWD 1000 50000 25 TRUE")

#scaffold firstbase lastbase nvar.all nvar.noCpG missing.percent
#NC_022208.1   1514798  1515797      140        127            1.71
#varseq.df[match(intersect(1514798:1515797, varpos), varpos), ]

#### SET-UP ####
cat("\n\nScript: createSeqs.R\n")
args <- commandArgs(trailingOnly = TRUE)

file.id <- args[1]
scaffold <- args[2]
basedir <- args[3]
if(length(args) >= 4) locus.length <- as.integer(args[4]) else locus.length <- 1000
if(length(args) >= 5) min.distance <- as.integer(args[5]) else min.distance <- 50000
if(length(args) >= 6) max.missing <- as.integer(args[6]) else max.missing <- 25
if(length(args) >= 7) mask.cpg <- args[7] else mask.cpg <- FALSE
if(length(args) >= 8) max.snps <- as.integer(args[8]) else max.snps <- 25

# file.id <- 'EjaC.Cgal.Cgui'; locus.length=1000; min.distance=50000; max.missing=25; mask.cpg=FALSE; max.snps=50; scaffold='NT_168049.1'; basedir='/home/jelmer/Dropbox/sc_fish/cichlids'

library(plyr)
setwd(basedir)

cat("File.id:", file.id, '\n')
cat("Scaffold:", scaffold, '\n')
cat("Locus length:", locus.length, '\n')
cat("Minimum distance between loci:", min.distance, '\n')
cat("Max % of missing data:", max.missing, '\n')
cat("Mask CpG:", mask.cpg, '\n')
cat("Max nr of SNPs in one locus:", max.snps, '\n')

#### FUNCTIONS ####

## Sliding window mean:
slideFunct <- function(data, window, stepsize){
  total <- length(data)
  spots <- seq(from = 1, to = (total - window), by = stepsize)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i] + window)])
  }
  return(result)
}

getSeqForInd <- function(ind.nr, Mask.cpg = mask.cpg) {
  varseq.ind <- as.character(varseq.df[, ind.nr])
  seq.ind <- ref
  seq.ind[varpos] <- varseq.ind
  if(Mask.cpg == TRUE) seq.ind[var.cpg] <- 'x' # REPLACE CpG by NN
  seq.ind <- gsub('c', 'C', seq.ind)
  seq.ind <- gsub('g', 'G', seq.ind)
  return(seq.ind)
}

getSeqForLoc <- function(seq.ind, locus.nr)
  return(paste(seq.ind[bed$firstbase[locus.nr]:bed$lastbase[locus.nr]], collapse = ''))

getLocus <- function(ind.nr, locus.nr)
  getSeqForLoc(get(paste0('seq.', inds[ind.nr])), locus.nr)

## Determine 1000bp seq with least Ns:
createLocus <- function(locus.nr) {
  loc.split <- strsplit(unlist(sapply(1:length(inds), getLocus, locus.nr)), split = '')
  loc.posindex <- bed$firstbase[locus.nr]:bed$lastbase[locus.nr]

  Ns <- lapply(loc.split, grepl, pattern = 'N')
  Ncount <- sapply(1:length(Ns[[1]]), FUN = function(x) sum(unlist(lapply(Ns, '[[', x))))
  slideMean <- slideFunct(data = Ncount, window = locus.length, stepsize = 1)

  firstbase.rel <- which(slideMean == min(slideMean))[1]
  lastbase.rel <- firstbase.rel + (locus.length - 1)
  firstbase.abs <- loc.posindex[firstbase.rel]
  lastbase.abs <- firstbase.abs + (locus.length - 1)

  varpos.locus.all <- intersect(seq(firstbase.abs, lastbase.abs), varpos)
  nvar.locus.all <- length(varpos.locus.all)
  var.cpg.locus <- intersect(seq(firstbase.abs, lastbase.abs), var.cpg)
  varpos.locus.noCpG <- setdiff(varpos.locus.all, var.cpg.locus)
  nvar.locus.noCpG <- length(varpos.locus.noCpG)
  if(mask.cpg == TRUE) nvar <- nvar.locus.noCpG
  if(mask.cpg == FALSE) nvar <- nvar.locus.all
  missing <- ifelse(nvar >= 1, round(sum(Ncount[firstbase.rel:lastbase.rel]) / ((length(loc.split) * nvar)) * 100, 2), 0)

  locus.stats[locus.nr, ] <<- c(scaffold, firstbase.abs, lastbase.abs, nvar.locus.all, nvar.locus.noCpG, missing)
  cat("Locus", locus.nr, "... startpos:", firstbase.abs, "/ all varpos:", nvar.locus.all,
      "/ varpos (no CpG):", nvar.locus.noCpG,  "/ % N in varpos:", missing, '\n')

  loc.split <- lapply(loc.split, gsub, pattern = 'x', replacement = 'N')

  loc.final <- unlist(lapply(lapply(loc.split, '[', firstbase.rel:lastbase.rel), paste, collapse = ''))
  loc.final <- paste(inds, loc.final)
  firstline <- paste(paste0(file.id, '_', scaffold, '_', firstbase.abs, '_', lastbase.abs), length(inds), locus.length)
  loc.final <- c(firstline, loc.final)

  return(loc.final)
}

locToExcl <- function(i, TooClose) {
  if(locus.stats$missing[TooClose[i]] <= locus.stats$missing[TooClose[i] + 1])
    return(TooClose[i] + 1) else return(TooClose[i])
}

getToExcl <- function(locus.stats) {
  dists <- unlist(sapply(1:nrow(locus.stats) - 1, FUN = function(x) locus.stats$firstbase[x + 1] - locus.stats$lastbase[x]))
  TooClose <- which(dists < 50000)
  if(length(TooClose) >= 1) return(unique(sapply(1:length(TooClose), locToExcl, TooClose)))
  if(length(TooClose) == 1) return(NULL)
}

subtractLoci <- function(locus.stats) {
  toExclude <- getToExcl(locus.stats)
  return(locus.stats[-toExclude, ])
}

count.alleles <- function(i) {
  seq.pos <- varseq.df[i, ]

  if(any(seq.pos == 'N')) seq.pos <- seq.pos[-grep('N', seq.pos)]
  seq.pos <- mapvalues(from = c('M', 'R', 'W', 'S', 'Y', 'K'), to = c('AC', 'AG', 'AT', 'CG', 'CT', 'GT'),
                       seq.pos, warn_missing = FALSE)
  seq.pos <- unlist(strsplit(seq.pos, split = ''))

  nr.alleles <- length(unique(seq.pos))

  return(nr.alleles)
}

count.minorAllele <- function(i) {
  seq.pos <- varseq.df[i, ]

  if(any(seq.pos == 'N')) seq.pos <- seq.pos[-grep('N', seq.pos)]
  seq.pos <- mapvalues(from = c('A', 'C', 'G', 'T'), to = c('AA', 'CC', 'GG', 'TT'),
                       seq.pos, warn_missing = FALSE)
  seq.pos <- mapvalues(from = c('M', 'R', 'W', 'S', 'Y', 'K'), to = c('AC', 'AG', 'AT', 'CG', 'CT', 'GT'),
                       seq.pos, warn_missing = FALSE)
  seq.pos <- unlist(strsplit(seq.pos, split = ''))

  nr.alleles <- length(unique(seq.pos))

  if(nr.alleles > 1) {
    minor.alle.count <- min(table(seq.pos))
    return(minor.alle.count)
  } else return(0)
}

get.minorAlleleInds <- function(i) {
  seq.pos <- varseq.df[i, ]

  seq.pos <- mapvalues(from = c('A', 'C', 'G', 'T'), to = c('AA', 'CC', 'GG', 'TT'),
                       seq.pos, warn_missing = FALSE)
  seq.pos <- mapvalues(from = c('M', 'R', 'W', 'S', 'Y', 'K'), to = c('AC', 'AG', 'AT', 'CG', 'CT', 'GT'),
                       seq.pos, warn_missing = FALSE)
  seq.pos <- unlist(strsplit(seq.pos, split = ''))

  major.allele <- names(table(seq.pos))[which(table(seq.pos) == max(table(seq.pos)))][1]

  if(major.allele != 'N') {
    seq.pos[grep('N', seq.pos)] <- major.allele
    #cat(major.allele, '\n')
    if(any(seq.pos != major.allele)) {
    minor.allele.indices <- which(seq.pos != major.allele)
    return(minor.allele.indices)
    }
  }
}



#### MAIN ####

## Create sequence:
bedfile.name <- paste0('analyses/gphocs/input_prep/', scaffold, '.bed')

if(file.exists(bedfile.name)) {
  bed <- read.table(bedfile.name, col.names = c('scaffold', 'firstbase', 'lastbase'))
  ref <- readRDS(paste0('analyses/gphocs/input_prep/', scaffold, '.hardmasked.rds')) # Created in chooseLoci.R
  varpos <- as.integer(read.table(paste0('seqdata/variants_otherFormats/fasta/', file.id, '.', scaffold, '.varpos'))$V2)
  cat("Head of varpos:", head(varpos), '\n')
  varseq <- readLines(paste0('seqdata/variants_otherFormats/fasta/', file.id, '.', scaffold, '.varscaffold'))
  inds <- sapply(strsplit(varseq, split = ' '), '[[', 1)
  varseq <- sapply(strsplit(varseq, split = ' '), '[[', 2)
  varseq.df <- matrix(unlist(sapply(1:length(varseq), FUN = function(x) strsplit(varseq[x], ''))), nrow = nchar(varseq[[1]]), byrow = FALSE)

  ## Remove invariable positions (only differ from ref):
  nr.alleles <- sapply(1:nrow(varseq.df), count.alleles)
  print(table(nr.alleles))
  nonvar <- which(nr.alleles < 2)
  if(length(nonvar) >= 1) {
    varpos <- varpos[-nonvar]
    varseq.df <- varseq.df[-nonvar, ]
  }

  maf.table <- table(unlist(sapply(1:nrow(varseq.df), get.minorAlleleInds)))
  maf.table <- maf.table[seq(1, 20, 2)] + maf.table[seq(2, 20, 2)]
  names(maf.table) <- inds
  cat('Minor allele occurences per ind:\n')
  print(maf.table)

  maf.table2 <- table(unlist(sapply(1:nrow(varseq.df), count.minorAllele)))
  cat('Minor allele freq per site:\n')
  print(maf.table2)

  if(length(unique(nchar(varseq))) != 1) cat("OOPS! Not all sequence strings of same length:", nchar(varseq), '\n')
  if(length(varpos) !=  nrow(varseq.df)) print("OOPS! Indices of variable positions do not have the same length as the sequences!\n")
  if(length(varpos) ==  nrow(varseq.df)) cat('Sequence length (variable positions):', length(varpos), '\n')

  var.cpg <- intersect(varpos, grep('c|g', ref)) # CpGs have been turned into lowercase c and g in previous script
  cat("Number of variable CpG sites:", length(var.cpg), '\n')

  ## Locus stats df:
  locus.stats <- data.frame(matrix(nrow = nrow(bed), ncol = 6))
  colnames(locus.stats) <- c('scaffold', 'firstbase', 'lastbase', 'nvar.all', 'nvar.noCpG', 'missing.percent')
  locus.stats$scaffold <- scaffold

  ## Get loci:
  for(ind.nr in 1:length(inds)) assign(paste0('seq.', inds[ind.nr]), getSeqForInd(ind.nr))
  loci <- lapply(1:nrow(bed), createLocus)

  ## Write all loci to file:
  write.table(locus.stats, paste0('analyses/gphocs/input_stats/locus_stats/allLoci/', scaffold, '.locus_stats.txt'), sep = '\t', quote = FALSE, row.names = FALSE)

  locus.stats$firstbase <- as.integer(locus.stats$firstbase)
  locus.stats$lastbase <- as.integer(locus.stats$lastbase)
  locus.stats$missing.percent <- as.numeric(locus.stats$missing.percent)
  locus.stats$missing.percent[is.na(locus.stats$missing.percent)] <- 0
  locus.stats$nvar.noCpG <- as.integer(locus.stats$nvar.noCpG)
  locus.stats$nvar.all <- as.integer(locus.stats$nvar.all)

  ## Loci to exclude due to too much missing data:
  toExclude.missing <- which(locus.stats$missing.percent >= max.missing)

  ## Loci to exclude due to too much many SNPs:
  if(mask.cpg == TRUE) toExclude.var <- which(locus.stats$nvar.noCpG > max.snps)
  if(mask.cpg == FALSE) toExclude.var <- which(locus.stats$nvar.all > max.snps)

  ## Exclude loci 1:
  toExclude1 <- unique(c(toExclude.var, toExclude.missing))
  if(length(toExclude1) >= 1) {
    locus.stats <- locus.stats[-toExclude1, ]
    rownames(locus.stats) <- 1:nrow(locus.stats)
    loci <- loci[-toExclude1]
  }

  ## Loci to exclude due to interlocus distance:
  if(!is.null(getToExcl(locus.stats))) locus.stats <- subtractLoci(locus.stats)
  if(!is.null(getToExcl(locus.stats))) locus.stats <- subtractLoci(locus.stats)
  if(!is.null(getToExcl(locus.stats))) locus.stats <- subtractLoci(locus.stats)
  if(!is.null(getToExcl(locus.stats))) locus.stats <- subtractLoci(locus.stats)
  if(!is.null(getToExcl(locus.stats))) locus.stats <- subtractLoci(locus.stats)
  if(!is.null(getToExcl(locus.stats))) locus.stats <- subtractLoci(locus.stats)
  toExclude.distance <- setdiff(1:length(loci), rownames(locus.stats))

  ## Exclude loci 2:
  if(length(toExclude.distance) >= 1) {
    loci <- loci[-toExclude.distance]
  }

  cat("NUMBER OF EXCLUDED LOCI DUE TO MISSING DATA:", length(toExclude.missing), '\n')
  cat("NUMBER OF EXCLUDED LOCI DUE TO TOO MANY SNPS:", length(toExclude.var), '\n')
  cat("MAX NR OF SNPS IN A LOCUS WITH/WITHOUT CPG:", max(locus.stats$nvar.all), max(locus.stats$nvar.noCpG), '\n')
  cat("NUMBER OF EXCLUDED LOCI DUE TO INTERLOCUS DISTANCE:", length(toExclude.distance), '\n')
  cat("FINAL NUMBER OF LOCI:", nrow(locus.stats), length(loci), '\n')

  loci <- unlist(loci)

  ## Write files:
  if(nrow(locus.stats) >= 1) {
    if(nrow(locus.stats) > 1) {
      closest <- min(unlist(sapply(1:nrow(locus.stats) - 1, FUN = function(x) locus.stats$firstbase[x + 1] - locus.stats$lastbase[x])))
      cat("SHORTEST DISTANCE BETWEEN TWO LOCI:", closest, '\n')
    }

    write.table(locus.stats, paste0('analyses/gphocs/input_stats/locus_stats/', file.id, '_', scaffold, '.locus_stats.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
    input.file <- paste0('analyses/gphocs/input_prep/', file.id, '_', scaffold, '.input')
    cat("Creating GPHOCS single-scaffold input file:", input.file, '\n')
    write(loci, input.file)
  }
} else {
  cat("No bed file", bedfile.name, "; scaffold probably contains no loci.\n")
  scafstats <- read.table(paste0('analyses/gphocs/scaffold_stats/', scaffold, '.maskStats.txt'), header = TRUE)
  print(scafstats)
}
