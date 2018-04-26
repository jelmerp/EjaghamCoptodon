#!/usr/bin/env Rscript

# Run using:
# Rscript scripts/trees/phylonet_prepnexus.R $TREEFILE $NEXUSFILE $PHYLONET.CMD $NR.RETIC $TAXONMAP 2> R-error-output.txt

###################################################################################################
##### SET-UP #####

## Set variables:
args <- commandArgs(trailingOnly = TRUE)
treefile <- args[1]
nexus.in <- args[2]
nexus.out <- args[3]
phylonet.cmd <- args[4]
taxonmap <- args[5]
focalnetworkname <- args[6]
focalnetwork <- args[7]

cat('Number of characters:', nchar(focalnetwork), "\n")

if(nchar(focalnetwork) == 0) {
  cat("No focal network defined...\n")
  focalnetworkname <- 'none'
  focalnetwork <- 'none'
}

cat('\n\n\nScript: phylonet_prepnexus.R\n')
cat('Tree file:', treefile, '\n')
cat('Nexus in:', nexus.in, '\n')
cat('Nexus out:', nexus.out, '\n')
cat('Phylonet command:', phylonet.cmd, '\n')
cat('Taxon map:', taxonmap, '\n')
cat('Focal network (for CalGTProb):', focalnetwork, '\n')
cat('Focal network name (for CalGTProb):', focalnetworkname, '\n')

## Get trees:
trees <- readLines(treefile)
nr.trees <- length(trees)
trees <- paste0('Tree genetree', 1:nr.trees, '=', trees)

## Define command line:
txt.out <- paste0(nexus.out, '.phylonet.output.txt')
nexus.out.full <- paste0(nexus.out, '.phylonet.output.nexus')
nexus.out.cmd <- paste0('Nexus_Out "', nexus.out.full, '";')

full.cmd <- paste0(phylonet.cmd, ' -a <', taxonmap, '> ', '"', txt.out, '"', ';')
cat('Full phylonet command:', full.cmd, '\n\n')

## Write nexus file:
cat('Writing nexus file:', nexus.in, '\n')

write('#NEXUS', nexus.in)
write('', nexus.in, append = TRUE)
if(focalnetwork != 'none') {
  write('BEGIN NETWORKS;', nexus.in, append = TRUE)
  write(paste0('Network ', focalnetworkname, '=', focalnetwork, ';'), nexus.in, append = TRUE)
  write('END;', nexus.in, append = TRUE)
  write('', nexus.in, append = TRUE)
}
write('BEGIN TREES;', nexus.in, append = TRUE)
write('', nexus.in, append = TRUE)
write(trees, nexus.in, append = TRUE)
write('', nexus.in, append = TRUE)
write('END;', nexus.in, append = TRUE)
write('', nexus.in, append = TRUE)
write('BEGIN PHYLONET;', nexus.in, append = TRUE)
write('', nexus.in, append = TRUE)
write(nexus.out.cmd, nexus.in, append = TRUE)
write(full.cmd, nexus.in, append = TRUE)
write('', nexus.in, append = TRUE)
write('END;', nexus.in, append = TRUE)
write('', nexus.in, append = TRUE)

cat('Done with phylonet_prepnexus.R\n')
