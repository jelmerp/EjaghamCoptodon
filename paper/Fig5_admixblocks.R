##### SET-UP #####

## Other scripts to load:
source('scripts/admixblocks/admixblocks_plotfun.R')

## Files to read:
admixblocks.hiConfidence.file <- 'analyses_output/admixblocks.highconfidence.txt'
admixblocks.likely.file <- 'analyses_output/admixblocks.likely.txt'
outputSummary.file <- 'analyses_output/gphocs_summary_theta_tau.txt'

## Files to write:
fig5.file <- 'figures/Fig5_admixblocks.png'

## Libraries:
library(ggpubr)
library(cowplot)

## Get blocks:
blocks.hic <- read.table(admixblocks.hiConfidence.file, header = TRUE, as.is = TRUE)
dfoil.lik <- read.delim(admixblocks.likely.file, header = TRUE, as.is = TRUE)

## Get ages:
tau <- read.table(outputSummary.file, header = TRUE, as.is = TRUE)
tau <- dplyr::filter(tau, var == 'tau')
tau2 <- filter(tau, migtoType3 == 'all')


##### FIG 5A - age of admixblocks #####
p <- ggplot()
p <- p + geom_violin(data = blocks.hic, aes(riverpop, age.mean / 1000, color = lakepop),
                     fill = "white", position = position_dodge(width = 0.9))
p <- p + geom_point(data = blocks.hic, aes(riverpop, age.mean / 1000, color = lakepop),
                    position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1))
p <- p + geom_rect(data = tau2, alpha = 0.35,
                   aes(xmin = -Inf, xmax = Inf, ymin = min / 1000, ymax = max / 1000, fill = pop))
p <- p + geom_hline(yintercept = DEF.split.point / 1000, colour = 'darkgreen')
p <- p + geom_hline(yintercept = DE.split.point / 1000, colour = '#0072EE')
p <- p + scale_x_discrete(labels = c(expression(italic('C. guineensis')),
                                     expression(italic('C. sp. Mamfé'))))
p <- p + scale_color_manual(values = sapply(c('Cdec', 'Ceja', 'Cfus'), getpopcol),
                            name = 'Lake species ',
                            labels = c(expression(italic('C. deckerti ')),
                                       expression(italic('C. ejagham ')),
                                       expression(italic('C. fusiforme'))))
p <- p + scale_fill_manual(values = sapply(as.character(tau2$pop), getpopcol, fcol = 'pop.gphocsfile2'))
p <- p + scale_y_continuous(expand = c(0, 0))
p <- p + labs(x = "Riverine source", y = "Age of block (ky)")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 16),
               axis.text.y = element_text(size = 16),
               axis.title.x = element_text(size = 18),
               axis.title.y = element_text(size = 18),
               legend.position = 'top',
               legend.title = element_text(size = 14, face = 'bold'),
               legend.text = element_text(size = 14),
               legend.margin = margin(0, 0, 0, 0),
               legend.box.margin = margin(-5, -5, -5, -5),
               legend.key.height = unit(0.5, "cm"),
               legend.key.width = unit(0.5, "cm"),
               plot.margin = margin(0.5, 1.25, 0.5, 0.75, "cm"))
p <- p + guides(fill = FALSE)
p <- p + coord_cartesian(ylim = c(0, 4.5))
p <- ggarrange(p, ncol = 1, nrow = 1)
p <- p + draw_plot_label(label = c('DE', 'DEF'), size = 16, hjust = 0,
                         colour = sapply(c('DE', 'DEF'), getpopcol, fcol = 'pop.gphocsfile2'),
                         x = c(0.925, 0.925), y = c(0.275, 0.395))
p.ages <- p


##### FIG 5B: Block sharing histogram #####
p <- ggplot()
p <- p + geom_histogram(data = blocks.hic, aes(lakepop, fill = unq),
                        stat = 'count', colour = 'grey20')
p <- p + scale_fill_discrete(name = 'Block sharing ',
                             labels = c('unique ', 'shared: 2 sp. ', 'shared: 3 sp.'))
p <- p + scale_x_discrete(labels = c(expression(italic('C. deckerti')),
                                     expression(italic('C. ejagham')),
                                     expression(italic('C. fusiforme'))))
p <- p + labs(x = "Lake Ejagham species", y = "Number of blocks")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
p <- p + theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))
p <- p + theme(legend.position = 'top')
p <- p + theme(legend.title = element_text(size = 14, face = 'bold'),
               legend.text = element_text(size = 14),
               legend.margin = margin(0, 0, 0, 0),
               legend.box.margin = margin(-5, -5, -5, -5),
               legend.key.height = unit(0.5, "cm"),
               legend.key.width = unit(0.5, "cm"),
               plot.margin = margin(0.5, 0.5, 0.5, 0.75, "cm"))
p.sharing <- p


##### FIG 5C: dfoil blocks #####
dfoil.lik$lakepoptype <- factor(dfoil.lik$lakepoptype, levels = c('extant species', 'ancestor DE', 'ancestor DEF'))

p <- ggplot()
p <- p + geom_histogram(data = dfoil.lik, aes(lakepoptype, fill = riverpop),
                        stat = 'count', colour = 'grey20')
p <- p + scale_fill_manual(name = 'Riverine source ',
                           values = sapply(c('Cgui', 'Cmam'), getpopcol),
                           labels = c(expression(italic('C. guineensis ')),
                                      expression(italic('C. sp. Mamfé '))))
p <- p + scale_x_discrete(labels = c('extant species', 'ancestor DE', 'ancestor DEF'))
p <- p + labs(x = "Lake Ejagham lineage", y = "Number of blocks")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 16),
               axis.text.y = element_text(size = 16),
               axis.title.x = element_text(size = 18),
               axis.title.y = element_text(size = 18),
               legend.position = 'top',
               legend.title = element_text(size = 14, face = 'bold'),
               legend.text = element_text(size = 14),
               legend.margin = margin(0, 0, 0, 0),
               legend.box.margin = margin(-5, -5, -5, -5),
               legend.key.height = unit(0.5, "cm"),
               legend.key.width = unit(0.5, "cm"),
               plot.margin = margin(0.75, 1.25, 0.5, 0.5, "cm"))
p.dfoil <- p


##### FIG 5D - example block #####
scaffold <- 'NC_022214.1'; xlims <- c(6, 7); blocklims <- c(6220001, 6670001) / 1000000

#filter(blocks.hic, ID == 'NC_022214.1_6220001_FEA') %>% select(contains('age'))

p.fd.a <- ggman(fd, pop = c('Eja.Dec.Mam', 'Fus.Dec.Mam', 'Fus.Eja.Mam'), scaffold = scaffold,
                drawpoints = FALSE, drawlines = TRUE, colvar.lines = 'pop', cols.lines = 'brew',
                smoothpar = 0.1, smoothmethod = 'loess',
                legpos = 'top', legtextsize = 14, plot.title = FALSE,
                my.xlim = xlims, blocklims = blocklims, plot.xann = FALSE, saveplot = FALSE)
p.fd.a <- p.fd.a + theme(plot.margin = margin(0.5, 0.5, 0.3, 0.75, "cm"))
p.fd.au <- ggman(fd, pop = c('Gui.Mam.Dec', 'Gui.Mam.Eja', 'Gui.Mam.Fus'), scaffold = scaffold,
                 drawpoints = FALSE, drawlines = TRUE, colvar.lines = 'pop',
                 cols.lines = as.character(sapply(c('Cdec', 'Ceja', 'Cfus'), getpopcol)),
                 smoothpar = 0.1, smoothmethod = 'loess',
                 legpos = 'top', legtextsize = 14,
                 my.xlim = xlims, blocklims = blocklims,
                 plot.xann = TRUE, xlab = paste0('position (Mb) along scaffold ', scaffold),
                 plot.title = FALSE, saveplot = FALSE)
p.fd.au <- p.fd.au + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.75, "cm"))
manplots <- ggarrange(p.fd.a, p.fd.au, ncol = 1, nrow = 2, heights = c(1, 1.2))


##### COMBINE PANELS #####
plots <- ggarrange(p.ages, p.sharing, p.dfoil, manplots, ncol = 2, nrow = 2)
plots <- plots + draw_plot_label(label = c('A', 'B', 'C', 'D'), size = 24,
                                 x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))
ggexport(plots, filename = fig5.file, width = 1000, height = 900)
