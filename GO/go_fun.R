##### SET-UP #####

## Files to read:
# (See scripts/go/go_makeAnnotTables.R for creation of these tables)
go.df.file <- 'analyses_input/GO_gene2goCat.txt'
annot.ncbi.file <- 'analyses_input/GO_annot.ncbi.df.txt'

## Libraries:
library(refGenome)
library(IRanges)
library(tidyverse)
library(goseq)

## Get GO and annotation tables:
go.df <- read.delim('', header = TRUE, stringsAsFactors = FALSE)
annot.ncbi <- read.delim(annot.ncbi.file, header = TRUE, stringsAsFactors = FALSE)


##### FUNCTIONS #####
## intersect.bed function: intersect bed file with reference (annotation) file.
## to.return = 'genes' [return the genes found in the ref] or 'lookup' [merge bed and ref for overlapping areas]
## ref = annotation or other geneset to use as reference
intersect.bed <- function(bed, ref = annot.ncbi, ref.cols = c('entrezgene', 'gene.name'),
                          to.return = 'genes', ignore.strand = TRUE) {

  cat('Intersecting', nrow(bed), 'regions...')

  genes.list <- apply(bed, 1, get.genes, ref = ref,
                      to.return = to.return, ref.cols = ref.cols, ignore.strand = ignore.strand)
  if(is.null(genes.list)) cat("No overlap found!\n")

  if(!is.null(genes.list)) {
    genes <- do.call(rbind, genes.list)

    if(to.return == 'genes') {
      genes <- unique(as.character(genes[, 1]))
      cat(length(genes), 'genes found\n')
    }
    if(to.return == 'lookup') {
      colnames(genes)[1:3] <- c('scaffold', 'start', 'end')
      genes <- arrange(genes, scaffold, start)
    }
    return(genes)
  }

}

## get.genes function: performs matching of one row of the bed file with the reference. Called in intersect.bed().
get.genes <- function(bed.row, ref, ref.cols, to.return, ignore.strand) {
  # bed.row <- bed[5, ]; ref = ann.bils.g; ignore.strand = TRUE

  bed.row <- data.frame(bed.row)
  if(nrow(bed.row) > 1) bed.row <- t(bed.row)

  ref <- ref[as.character(ref$scaffold) == as.character(bed.row[1, 1]), ] # subset ref df to correct scaffold
  if(ignore.strand == FALSE) ref <- ref[as.character(ref$strand) == as.character(bed.row[1, 4]), ] # subset ref df to correct strand

  ref.range <- IRanges(start = ref$start, end = ref$end)
  lift.range <- IRanges(start = as.integer(bed.row[1, 2]), end = as.integer(bed.row[1, 3]))

  ov <- as.data.frame(findOverlaps(ref.range, lift.range))[, 1] # indices of ref rows that match

  if(length(ov) >= 1) {
    if(to.return == 'genes') one.row <- data.frame(ref[ov, ref.cols], row.names = NULL)
    if(to.return == 'lookup') one.row <- data.frame(cbind(bed.row, ref[ov, ref.cols], row.names = NULL))
    return(one.row)
  }

}

## GO_geneset: GO or KEGG test on any gene selection.
runGO <- function(focgenes, refgenes = annot.ncbi$entrezgene, go.map = go.df,
                  write.table = FALSE, filename = NULL, setname = NA) {
  # focgenes <- gfoc.ncbi$entrezgene; refgenes <- annot.ncbi$entrezgene; go.map <- go.df; setname = 'allblocks'
  cat('Test', length(focgenes), 'specified genes...', focgenes[1], 'etc \n')

  ## GO map:
  go.map <- go.map %>% dplyr::filter(!(is.na(go.id))) %>% select(entrezgene, go.id)

  ## Focal genes:
  focgenes <- unique(as.character(focgenes)[!is.na(focgenes)])

  ## Reference genes:
  refgenes <- unique(as.character(refgenes[!is.na(refgenes)]))
  refgenes <- refgenes[refgenes %in% go.df$entrezgene] # Only keep genes with GO mapping
  refgenes <- refgenes[-which(refgenes %in% focgenes)] # Get rid of focal genes

  ## Probability weighting function:
  pwf <- data.frame(DEgenes = c(rep(1, length(focgenes)), rep(0, length(refgenes))),
                    bias.data = 1, pwf = 1)
  rownames(pwf) <- c(focgenes, refgenes)

  ## Perform test:
  go.res <- goseq(pwf = pwf, gene2cat = go.map, method = "Wallenius")

  ## Prep table:
  colnames(go.res) <- c("category", "p.over", "p.under", 'numFocInCat', 'numInCat',
                        'term', 'ontology')
  go.res$FDR.over <- p.adjust(go.res$p.over, method = "BH")
  go.res$FDR.under <- p.adjust(go.res$p.under, method = "BH")
  go.res <- go.res %>% select(category, p.over, FDR.over, p.under, FDR.under,
                              numFocInCat, numInCat, term, ontology)
  go.res$setname <- setname

  ## Significant results:
  sig <- go.res[go.res$FDR.over < 0.05, ]

  if(nrow(sig) == 0) {
    cat("\t\t\tNO SIGNIFICANT GO CATEGORIES\n")
  } else {
    cat(nrow(sig), 'significant categories\n')
    print(sig$term)
  }

  if(write.table == TRUE) {
    if(is.null(filename)) filename <- setname
    filename.full <- paste0('analyses/GO/output/go.res_', filename, '.txt')
    cat('Writing table:', filename, '\t')
    write.table(sig, quote = FALSE, sep = "\t", row.names = FALSE)
  }

  return(sig)
}
