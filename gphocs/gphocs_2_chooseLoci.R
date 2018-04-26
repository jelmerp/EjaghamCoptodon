#!/usr/bin/env Rscript

## Usage: Rscript /home/jelmer/Dropbox/sc_fish/cichlids/scripts/gphocs/chooseLoci.R NC_022219.1 /home/jelmer/Dropbox/sc_fish/cichlids/

## MASK INDELS??
## MASK FILTERED SNPS??

## Set-up:
cat("\n\nScript: chooseLoci.R\n")
args <- commandArgs(trailingOnly = TRUE) 
scaffold <- args[1]
basedir <- args[2]
if(length(args) >= 3) exon.buffer <- as.integer(args[3]) else exon.buffer <- 500 # Number of bases on each side of exon that should be masked

# rm(list = ls()); scaffold <- 'NT_167980.1'; basedir <- '/home/jelmer/Dropbox/sc_fish/cichlids/'; exon.buffer <- 500

cat("Scaffold:", scaffold, '\n')
cat("Exon buffer:", exon.buffer, '\n')
setwd(basedir)

## Fasta reference sequence:
ref <- tail(readLines(paste0('seqdata/reference/byScaffold/', scaffold, '.fa')), -1)
ref <- unlist(strsplit(paste(ref, collapse = ''), split = ''))
scaffold.length <- length(ref)

## Repeats in reference sequence (lowercase):
lowercase.ref <- grep('[a-z]', ref) # lowercase
N.ref <- grep('N', ref)

## SNPable poor regions:
SNPable <- tail(readLines(paste0('analyses/SNPable/byScaffold/', scaffold, '.mask')), -1)
SNPable <- unlist(strsplit(paste(SNPable, collapse = ''), split = ''))
SNPable.excl  <- grep('0|1', SNPable)

## Exons:
exons.df <- read.delim(gzfile(paste0('seqdata/annot/byScaffold/', scaffold, '.gff.gz')), skip = 1, header = F)[, 3:5]
colnames(exons.df) <- c('element', 'firstbase', 'lastbase')
exons.df <- exons.df[exons.df$element == 'exon', ]
exons <- unlist(sapply(1:nrow(exons.df), FUN = function(x) seq(from = exons.df$firstbase[x] - exon.buffer, to = exons.df$lastbase[x] + exon.buffer)))
if(any(exons < 0)) exons <- exons[-which(exons < 0)]
if(any(exons > scaffold.length)) exons <- exons[-which(exons > scaffold.length)]
exons <- sort(unique(exons))

## Combine bad bases and turn into Ns:
exclude.hard <- sort(unique(c(lowercase.ref, N.ref, SNPable.excl, exons)))
ref[exclude.hard] <- 'N'

## Report:
prop.lowercase <- round(length(lowercase.ref) / scaffold.length, 4)
prop.N <- round(length(N.ref) / scaffold.length, 4)
prop.exon <- round(length(exons) / scaffold.length, 4)
prop.SNPable <- round(length(SNPable.excl) / scaffold.length, 4)
prop.totalmask <- round(length(exclude.hard) / scaffold.length, 4)

cat("Scaffold length:", scaffold.length, '\n')
cat("Proportion soft masked bases (lowercase) in reference genome scaffold:", prop.lowercase, "\n")
cat("Proportion hard masked bases (Ns) in reference genome scaffold:", prop.N, "\n")
cat("Proportion of bases in exons:", prop.exon, "\n")
cat("Proportion of bases excluded by SNPable:", prop.SNPable, "\n")
cat("TOTAL PROPORTION OF EXCLUDED BASES:", prop.totalmask, "\n\n")

## Find GPHOCS loci:
cat("Finding suitable regions for GPHOCS loci..\n")
N.locs <- grep('N', ref)
N.intervals <- diff(N.locs) # interval lengths between Ns
if(max(N.intervals > 1001)) {
  potloc <- which(N.intervals > 1001) # intervals between Ns longer than 1kb
  firstbase <- N.locs[potloc] + 1
  lastbase <- N.locs[potloc + 1] - 1
  
  bed <- data.frame(scaffold, firstbase, lastbase)
  bedfile <- paste0('analyses/gphocs/input_prep/', scaffold, '.bed')
  cat("NUMBER OF LOCI:", nrow(bed), "\n")
  cat("Writing bedfile:", bedfile, "\n")
  write.table(bed, bedfile, sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
  nloci <- nrow(bed)
} else {
  cat('\n\nNO SUITABLE LOCI FOUND! Largest interval with no Ns:', max(N.intervals), '\n\n')
  nloci <- 0
}
  
## CpG sites:
ref <- unlist(strsplit(gsub('CG', 'cg', paste(ref, collapse = '')), split = ''))
CpG.count <- length(grep('c', ref))
cat("CpG count:", CpG.count, "\n")

## Save scaffold sequence file:
RDS.file <- paste0('analyses/gphocs/input_prep/', scaffold, '.hardmasked.rds')
cat("Saving RDS file with final masked sequences:", RDS.file, '\n\n')
saveRDS(ref, file = RDS.file)

## Write stats to table:
id.column <- c('prop.lowercase', 'prop.N', 'prop.exon', 'prop.SNPable', 'prop.totalmask', 'scaffold.length', 'nloci_unfiltered')
value.column <- c(prop.lowercase, prop.N, prop.exon, prop.SNPable, prop.totalmask, scaffold.length, nloci)
df <- data.frame(scaffold = scaffold, id.column, value.column)
write.table(df, paste0('analyses/gphocs/input_stats/scaffold_stats/', scaffold, '.maskStats.txt'), sep = '\t', quote = FALSE, row.names = FALSE)

cat("Done with script.\n")

###################################################################################################################

## Split mask by scaffold:
#cd /home/jelmer/Dropbox/sc_fish/cichlids/analyses/SNPable/byScaffold
#/home/jelmer/Dropbox/sc_fish/cichlids/scripts/conversion/fastaSplitter.pl /home/jelmer/Dropbox/sc_fish/cichlids/analyses/SNPable/mask_l50_r50.fa

## Rename mask files:
# names.old <- list.files('analyses/SNPable/byScaffold', full.names = TRUE)
# names.new <- paste0(sapply(strsplit(names.old, split = ' '), '[[', 1), '.mask')
# file.rename(names.old, names.new)
