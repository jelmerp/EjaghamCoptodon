##### SET-UP #####
options(scipen = 1)

## Other scripts:
source('scripts/GO/go_fun.R')

## Files to write:
table2.file <- 'tables/Table2_GO.txt'


##### TABLE 2: GO #####
table2 <- blocks.hic %>%
  select(scaffold, start, end) %>%
  intersect.bed() %>%
  runGO() %>%
  select(ontology, category, term, FDR.over, numFocInCat) %>%
  dplyr::rename(FDR = FDR.over, ngenes = numFocInCat)

blocks.hic.dec.go <- blocks.hic %>% filter(lakepop == 'Cdec') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.eja.go <- blocks.hic %>% filter(lakepop == 'Ceja') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.fus.go <- blocks.hic %>% filter(lakepop == 'Cfus') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.unq.go <- blocks.hic %>% filter(unq == 'unq') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.shared1.go <- blocks.hic %>% filter(unq == 'shared1') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()
blocks.hic.shared2.go <- blocks.hic %>% filter(unq == 'shared2') %>% select(scaffold, start, end) %>% intersect.bed() %>% runGO()

table2$C.deckerti <- as.integer(table2$category %in% blocks.hic.dec.go$category)
table2$C.ejagham <- as.integer(table2$category %in% blocks.hic.eja.go$category)
table2$C.fusiforme <- as.integer(table2$category %in% blocks.hic.fus.go$category)
table2$unique <- as.integer(table2$category %in% blocks.hic.unq.go$category)
table2$shared1 <- as.integer(table2$category %in% blocks.hic.shared1.go$category)
table2$shared2 <- as.integer(table2$category %in% blocks.hic.shared2.go$category)

write.table(table2, table2.file, sep = '\t', quote = FALSE, row.names = FALSE)
