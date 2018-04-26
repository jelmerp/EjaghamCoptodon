##### SET-UP #####

## Other scripts to load:
source('scripts/admixblocks/admixblocks_plotfun.R')

## Files to read:
admixblocks.hiConfidence.file <- 'analyses_output/admixblocks.highconfidence.txt'
admixblocks.likely.file <- 'analyses_output/admixblocks.likely.txt'
outputSummary.file <- 'analyses_output/gphocs_summary_theta_tau.txt'

## Files to write:
admixblocks.suppfig.file <- 'figures/SuppFig_admixblocks.png'

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


##### PANEL A: AGE OF ALL BLOCKS #####
p <- ggplot()
p <- p + geom_violin(data = blocks.all, aes(riverpop, age.mean / 1000, fill = lakepop),
                     position = position_dodge(width=0.9))
p <- p + geom_hline(yintercept = root.split.point / 1000, colour = getpopcol('root'))
p <- p + geom_hline(yintercept = DEF.split.point / 1000, colour = getpopcol('DEF', fcol = 'pop.gphocsfile2'))
p <- p + scale_x_discrete(labels = c(expression(italic('C. guineensis')),
                                     expression(italic('C. sp. Mamfé'))))
p <- p + scale_fill_manual(values = sapply(c('Cdec', 'Ceja', 'Cfus'), getpopcol),
                            name = 'Lake species ',
                            labels = c(expression(italic('C. deckerti ')),
                                       expression(italic('C. ejagham ')),
                                       expression(italic('C. fusiforme'))))
p <- p + scale_y_continuous(expand = c(0, 0))
p <- p + labs(x = 'Riverine source', y = "Age of block (ky)")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 16),
               axis.text.y = element_text(size = 16),
               axis.title.x = element_text(size = 18),
               axis.title.y = element_text(size = 18),
               legend.position = 'top',
               legend.title = element_text(size = 14, face = 'bold'),
               legend.text = element_text(size = 14),
               legend.key.height = unit(0.5, "cm"),
               legend.key.width = unit(0.5, "cm"),
               legend.margin = margin(0, 0, 0, 0),
               legend.box.margin = margin(0, 0, -5, 0),
               plot.margin = margin(0.5, 1.5, 0.5, 0.5, "cm"))
p <- ggarrange(p, ncol = 1, nrow = 1)
p <- p + draw_plot_label(label = 'DEF', x = 0.915, y = 0.135, size = 16, hjust = 0,
                    colour = getpopcol('DEF', fcol = 'pop.gphocsfile2'))
p <- p + draw_plot_label(label = 'root', x = 0.915, y = 0.225, size = 16, hjust = 0,
                    colour = getpopcol('root', fcol = 'pop.gphocsfile2'))
panelA <- p


##### PANEL B: AGE OF UNIQUE VS. SHARED BLOCKS #####
p <- ggplot()
p <- p + geom_violin(data = blocks.hic, aes(unq, age.mean / 1000, color = riverpop),
                     fill = "white", position = position_dodge(width=0.9))
p <- p + geom_point(data = blocks.hic, aes(unq, age.mean / 1000, color = riverpop),
                    position = position_jitterdodge(dodge.width = 0.9))
p <- p + geom_rect(data = tau2, alpha = 0.35,
                   aes(xmin = -Inf, xmax = Inf, ymin = min / 1000, ymax = max / 1000, fill = pop))
p <- p + geom_hline(yintercept = DEF.split.point / 1000, colour = 'darkgreen')
p <- p + geom_hline(yintercept = DE.split.point / 1000, colour = '#0072EE')
p <- p + scale_fill_manual(values = sapply(as.character(tau2$pop), getpopcol, fcol = 'pop.gphocsfile2'))
p <- p + scale_x_discrete(labels = c('unique', 'shared:\n 2 species', 'shared:\n 3 species'))
p <- p + scale_color_manual(values = sapply(c('Cgui', 'Cmam'), getpopcol),
                            name = 'Riverine source ',
                            labels = c(expression(italic('C. guineensis ')),
                                       expression(italic('C. sp. Mamfé'))))
p <- p + scale_y_continuous(expand = c(0, 0))
p <- p + labs(x = 'Block sharing among lake species', y = "Age of block (ky)")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 16),
               axis.text.y = element_text(size = 16),
               axis.title.x =  element_text(size = 18),
               axis.title.y = element_text(size = 18),
               legend.position = 'top',
               legend.title = element_text(size = 14, face = 'bold'),
               legend.text = element_text(size = 14),
               legend.key.height = unit(0.5, "cm"),
               legend.key.width = unit(0.5, "cm"),
               legend.margin = margin(0, 0, 0, 0),
               legend.box.margin = margin(0, 0, -5, 0),
               plot.margin = margin(0.5, 1.5, 0.5, 0.5, "cm"))
p <- p + guides(fill = FALSE)
p <- p + coord_cartesian(ylim = c(0, 4.5))
p <- ggarrange(p, ncol = 1, nrow = 1)
p <- p + draw_plot_label(label = 'DE', x = 0.915, y = 0.27, size = 16, hjust = 0,
                         colour = getpopcol('DE', fcol = 'pop.gphocsfile2'))
p <- p + draw_plot_label(label = 'DEF', x = 0.915, y = 0.39, size = 16, hjust = 0,
                         colour = getpopcol('DEF', fcol = 'pop.gphocsfile2'))
panelB <- p


##### COMBINE PANELS #####
panelA <- panelA + theme(plot.margin = margin(0.5, 0.5, 1, 0.5, "cm"))
panelB <- panelB + theme(plot.margin = margin(0.5, 0.5, 1, 0.5, "cm"))

plots <- ggarrange(panelA, panelB, ncol = 2, nrow = 1)
plots <- plots + draw_plot_label(label = c('A', 'B'), size = 24, x = c(0, 0.5), y = c(1, 1))
ggexport(plots, filename = admixblocks.suppfig.file, width = 1000, height = 600)
