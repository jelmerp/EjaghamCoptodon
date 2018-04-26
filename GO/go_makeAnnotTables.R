##### SET-UP #####

## Files to read:
gff.ncbi.file <- 'seqdata/annot/GCF_000188235.2_Orenil1.1_genomic.gff'


## Files to write:
annot.ncbi.file <- 'analyses_input/GO_annot.ncbi.df.txt'
gene2goCat.file <- 'analyses_input/GO_gene2goCat.txt'

## Libraries:
library(biomaRt)


##### FUNCTIONS #####

## Function to turn GFF files into dataframes:
# collect: fields from the 9th column to extract, e.g. geneID, gene.name, etc.
# col.names: column names to give to the extracted fields, by default the same as the names for "collect".
# type: the annotation elements to select from column 3, e.g. "exon", "gene", "mRNA", etc.
ggf2df <- function(gff.file, collect = c('ID', 'geneID', 'gene.name'), col.names = NULL, type = 'all') {
  # type = 'gene'; collect = c('ID', 'Dbxref', 'Name', 'gene')

  ## If no col.names are provided, they should be the same as "collect" (which fields to collect from the gff):
  if(is.null(col.names)) col.names <- collect

  ## Read file:
  annot <- read.delim(gff.file, head = FALSE, skip = 7)
  if(ncol(annot) != 9) cat("WARNING: FILE DOES NOT HAVE 9 COLUMNS")

  ## Select type of element (all, or e.g. only 'gene' or 'exon' or 'mRNA' from column 3):
  if(type != 'all') annot <- annot[annot[, 3] == type, ]

  ## Collect desired date from column 9, which has fields separated by semicolons:
  annot9 <- strsplit(as.character(annot[, 9]), ';')

  collect_function_1 <- function(collect_item) as.character(sapply(annot9, collect_function_2, collect_item))

  collect_function_2 <- function(annot9_row, collect_item) {
    element <- grep(paste0('^', collect_item, '='), unlist(annot9_row), value = TRUE)
    if (length(element) == 0) element <- c('NA')
    unlist(strsplit(element, '='))[length(unlist(strsplit(element, '=')))]
  }

  annotAttr <- lapply(collect, collect_function_1)
  annot9 <- do.call(cbind, annotAttr)
  colnames(annot9) <- collect

  ## Get columns 1-8 and combine with data extracted from column 9:
  annot1to8 <- annot[, 1:8]
  annot <- cbind(annot1to8, annot9)

  ## Column names:
  colnames(annot)[1:8] <- c('scaffold', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase')
  colnames(annot)[9:ncol(annot)] <- col.names

  return(annot)
}


##### PROCESS NCBI GFF #####
annot.ncbi <- ggf2df(gff.ncbi.file, type = 'gene', collect = c('ID', 'Dbxref', 'Name'),
                col.names = c('gene.id', 'entrezgene', 'gene.name'))
write.table(annot.ncbi, annot.ncbi.file, sep = '\t', quote = FALSE, row.names = FALSE)


##### GET BIOMART ANNOTATIONS #####
mart.oni <- useMart('ensembl', dataset = 'oniloticus_gene_ensembl')
go <- getBM(mart = mart.oni, attributes = c('ensembl_gene_id', 'go_id', 'name_1006',
                                            'namespace_1003', 'external_gene_name',
                                            'hgnc_symbol', 'hgnc_id', 'entrezgene'))
go$entrezgene <- paste0('GeneID:', go$entrezgene)
colnames(go) <- c('ensembl.gene.id', 'go.id', 'go.name', 'go.ontology',
                  'gene.name', 'hgnc.symbol', 'hgnc.id', 'entrezgene')
go$go.id[grep('^$', go$go.id)] <- NA
go$go.name[grep('^$', go$go.name)] <- NA
go$gene.name[grep('^$', go$gene.name)] <- NA
go$go.ontology[grep('^$', go$go.ontology)] <- NA
go$hgnc.symbol[grep('^$', go$hgnc.symbol)] <- NA
go$hgnc.id[grep('^$', go$hgnc.id)] <- NA
write.table(go, gene2goCat.file, sep = '\t', quote = FALSE, row.names = FALSE)
