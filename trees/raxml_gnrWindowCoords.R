#!/usr/bin/env Rscript

# Run using:
# Rscript scripts/conversion/gnrWindowCoords.R $WINDOW.SIZE 2> R-error-output.txt

###################################################################################################
##### SET-UP #####

## Set variables:
args <- commandArgs(trailingOnly = TRUE)
window.size <- as.integer(as.character(args[1])) # window.size <- 100

cat("Window size:", window.size, '\n')

#rm(list = ls()); gc()
#setwd('/home/jelmer/Dropbox/sc_fish/cichlids/')

###################################################################################################
##### FUNCTIONS #####

qpScaf <- function(line) {
  cat('Line number:', line, 'Scaffold name:', as.character(sc$scaffold.name[line]))
  scaflength <- sc$scaffold.length[line]

  if(scaflength > 100000) {
    start <- seq(0, scaflength, by = window.size * 1000)
    start <- start[-length(start)] # get rid of last starting point to avoid getting outside of scaffold boundary
    end <- start + (window.size * 1000) - 1
    sc.df <- data.frame(scaffold = sc$scaffold.name[line], start, end)
    cat(' Number of windows:', nrow(sc.df), '\n')
    return(sc.df)
  } else {
    cat(' Scaffold length smaller than window length \n')
  }
}

###################################################################################################
##### CALCULATE WINDOWS #####
sc <- read.table('metadata/scaffolds_withLength.txt', header = TRUE)
sc.df <- do.call(rbind, lapply(1:nrow(sc), qpScaf))
sc.df$start <- format(as.integer(sc.df$start), scientific = FALSE)

filename <- paste0('metadata/scaffoldWindows_', window.size, 'kb.txt')
write.table(sc.df, filename, sep ='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
