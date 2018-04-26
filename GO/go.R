##### SET-UP #####

## Other scripts:
source('scripts/GO/go_fun.R')

## Files to read:
admixblocks.hiConfidence.file <- 'analyses_output/admixblocks.highconfidence.txt'
admixblocks.likely.file <- 'analyses_output/admixblocks.likely.txt'

## Get focal genes:
blocks.hic <- read.table(admixblocks.hiConfidence.file, header = TRUE, as.is = TRUE)
dfoil.lik <- read.delim(admixblocks.likely.file, header = TRUE, as.is = TRUE)


##### APPLY #####
blocks.lik.go <- blocks.lik %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.go <- blocks.hic %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()

## Blocks by species:
blocks.hic.dec.go <- blocks.hic %>% filter(lakepop == 'Cdec') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.eja.go <- blocks.hic %>% filter(lakepop == 'Ceja') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.fus.go <- blocks.hic %>% filter(lakepop == 'Cfus') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()

nrow(blocks.hic.eja.go)
blocks.hic.eja.go$category %in% blocks.hic.go$category
nrow(blocks.hic.dec.go)
blocks.hic.dec.go$category %in% blocks.hic.go$category

## Among unique/shared blocks:
blocks.hic.unq.go <- blocks.hic %>% filter(unq == 'unq') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.shared1.go <- blocks.hic %>% filter(unq == 'shared1') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.shared2.go <- blocks.hic %>% filter(unq == 'shared2') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.shared1.noCfus.go <- blocks.hic %>% filter(unq == 'shared1', lakepop != 'Cfus') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()

nrow(blocks.hic.shared1.go)
blocks.hic.shared1.go$category %in% blocks.hic.go$category

nrow(blocks.hic.shared2.go)
blocks.hic.shared2.go$category %in% blocks.hic.go$category

## No significant categories among "unique" admix blocks by species:
blocks.hic.unq.Dec.go <- blocks.hic %>% filter(unq == 'unq', lakepop == 'Cdec') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.unq.Eja.go <- blocks.hic %>% filter(unq == 'unq', lakepop == 'Ceja') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()

## Which genes:
blocks.hic.g <- blocks.hic %>% select(scaffold, start, end) %>% intersect.bed()
g.lookup <- intersect.bed(select(blocks.hic, scaffold, start, end, ID), to.return = 'lookup')
