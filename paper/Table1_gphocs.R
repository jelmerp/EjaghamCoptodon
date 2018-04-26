##### SET-UP #####
options(scipen = 1)

## Other scripts:
source('scripts/gphocs/gphocs_analyze_fun.R')

## Files to write:
table1a.file <- 'tables/Table1a_gphocs.txt'
table1b.file <- 'tables/Table1b_gphocs.txt'
table1c.file <- 'tables/Table1c_gphocs.txt'


##### TABLE 1: G-PHOCS #####
Log$migtoType3 <- gsub( 'ancestral/extant', 'anc.ext', Log$migtoType3)

## Summary for migration from Mam -- tau and theta:
sumr.tt <- dplyr::filter(Log, var %in% c('tau', 'theta'), migfromRun %in% c('none', 'Mam')) %>%
  dplyr::group_by(migtoType3, var, pop) %>%
  dplyr::summarise(mean = round(mean(cval)), min = round(hpd.min(cval)), max = round(hpd.max(cval))) %>%
  dplyr::filter(!(var %in% c('theta') & pop %in% c('root', 'AU', 'Gui', 'Mam'))) %>%
  data.table() %>%
  data.table::dcast(var + pop ~ migtoType3, value.var = c('mean', 'min', 'max')) %>%
  mutate(range_none = paste0(min_none, '-', max_none),
         range_single = paste0(min_single, '-', max_single),
         range_anc.ext = paste0(min_anc.ext, '-', max_anc.ext),
         range_all = paste0(min_all, '-', max_all)) %>%
  select(var, pop, mean_single, mean_anc.ext, mean_all, mean_none,
         range_single, range_anc.ext, range_all, range_none)
write.table(sumr.tt, table1a.file, sep  = '\t', quote = FALSE, row.names = FALSE)

## Summary for migration from Mam -- 2Nm and total migrate:
sumr.2Nm <- dplyr::filter(Log, var == '2Nm', migfromRun == 'Mam') %>%
  dplyr::group_by(migtoType3, migto, var) %>%
  dplyr::summarise(mean = round(mean(val), 2), min = round(hpd.min(val), 2), max = round(hpd.max(val), 2))
sumr.mtot <- dplyr::filter(Log, var == 'm', migfromRun == 'Mam') %>%
  dplyr::group_by(migtoType3, migto, var) %>%
  dplyr::summarise(mean = round(mean(cval), 2), min = round(hpd.min(cval), 2), max = round(hpd.max(cval), 2))
sumr.m1 <- rbind(sumr.2Nm, sumr.mtot) %>%
  dplyr::filter(migto != 'Gui') %>%
  dplyr::rename(pop = migto) %>%
  data.table() %>%
  data.table::dcast(var + pop ~ migtoType3, value.var = c('mean', 'min', 'max')) %>%
  mutate(range_single = paste0(min_single, '-', max_single),
         range_anc.ext = paste0(min_anc.ext, '-', max_anc.ext),
         range_all = paste0(min_all, '-', max_all)) %>%
  select(var, pop, mean_single, mean_anc.ext, mean_all,
         range_single, range_anc.ext, range_all)
write.table(sumr.m1, table1b.file, sep  = '\t', quote = FALSE, row.names = FALSE)

## Prop migrants:
sumr.mprop <- dplyr::filter(Log, var == 'm.prop', migfromRun == 'Mam') %>%
  dplyr::group_by(migtoType3, migto, var) %>%
  dplyr::summarise(mean = formatC(mean(cval), format = "e", digits = 2),
                   min = formatC(hpd.min(cval), format = "e", digits = 2),
                   max = formatC(hpd.max(cval), format = "e", digits = 2)) %>%
  dplyr::filter(migto != 'Gui') %>%
  dplyr::rename(pop = migto) %>%
  data.table() %>%
  data.table::dcast(var + pop ~ migtoType3, value.var = c('mean', 'min', 'max')) %>%
  mutate(range_single = paste0(min_single, '-', max_single),
         range_anc.ext = paste0(min_anc.ext, '-', max_anc.ext),
         range_all = paste0(min_all, '-', max_all)) %>%
  select(var, pop, mean_single, mean_anc.ext, mean_all,
         range_single, range_anc.ext, range_all)
write.table(sumr.mprop, table1c.file, sep  = '\t', quote = FALSE, row.names = FALSE)

