##### SETUP #####

## Files to read:
popcols.file <- 'metadata/popcolours.txt'
raxml_tree.file <- 'analyses_output/RAxML_IC_Score_BranchLabels.EjaC.Ckot.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.ICAscores.MLtree'
saguaro_output.summary.file <- 'analyses_output/saguaro_output_summary.txt'
saguaro_cactus1.file <- 'analyses_output/saguaro_cactus1.tre'
saguaro_cactus3.file <- 'analyses_output/saguaro_cactus3.tre'

## Files to write:
fig2.file <- 'figures/Fig2.png'

## Libraries:
library(ape)
library(tidyverse)
library(ggtree)
library(ggpubr)
library(cowplot)

## Read general files:
inds <- read.table(popcols.file, header = TRUE, stringsAsFactors = FALSE)

## Functions:
getpopcol <- function(pop, fcol = 'pop.gphocsfile1') inds$pop.col[match(pop, inds[, fcol])]


##### PANEL A: TREE WITH ICA SCORES #####

## EjaC.Ckot:
raxml <- read.raxml(raxml_tree.file)
raxml@data$bootstrap <- raxml@data$bootstrap / 100
p <- ggtree(raxml, size = 1)
p <- p + geom_label2(aes(label = bootstrap, fill = bootstrap, subset=(node!=19)), size = 3)
p <- p + scale_fill_continuous(limits = c(0, 1), low = 'red', high = 'darkgreen', name = 'ICA score')
p <- p + theme_tree(legend.position = c(0.15, 0.2), legend.title = element_text(size = 12))
p <- p + theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))
p <- p + ggplot2::xlim(-0.03, 0.5)
p <- p + geom_cladelabel(node = 12, align = TRUE, color = 'black', fontsize = 5, offset = 0.01,
                         label = 'paste(italic("C. kottae"))', parse = TRUE)
p <- p + geom_cladelabel(node = 1, align = TRUE, color = 'orangered1', fontsize = 5, offset = 0.01,
                         label = 'paste(italic("C. sp. Mamfé"))', parse = TRUE)
p <- p + geom_cladelabel(node = 16, align = TRUE, color = 'orange2', fontsize = 5, offset = 0.01,
                         label = 'paste(italic("C. guineensis"))', parse = TRUE)
p <- p + geom_cladelabel(node = 18, align = TRUE, color = 'green3', fontsize = 5, offset = 0.01,
                         label = 'paste(italic("C. fusiforme"))', parse = TRUE)
p <- p + geom_cladelabel(node = 22, align = TRUE, color = 'turquoise2', fontsize = 5, offset = 0.01,
                         label = 'paste(italic("C. deckerti"))', parse = TRUE)
p <- p + geom_cladelabel(node = 21, align = TRUE, color = 'blue2', fontsize = 5, offset = 0.01,
                         label = 'paste(italic("C. ejagham"))', parse = TRUE)
p <- p + annotate(geom = "text", x = -0.02, y = 13, label = "A", fontface = 'bold', color = 'black', size = 8)
trueboot <- subset(p$data, isTip == FALSE)
trueboot$truebootstrap1 <- c(NA, '*', NA, '*', '*', '*', NA, '*', '*', '*')
trueboot$truebootstrap2 <- c(NA, NA, 70, NA, NA, NA, NA, NA, NA, NA)
p <- p + geom_text(data = trueboot, aes(label = trueboot$truebootstrap1), size = 7, hjust = 1.2, vjust = -0.05)
p <- p + geom_text(data = trueboot, aes(label = trueboot$truebootstrap2), size = 3.5, hjust = 1.2, vjust = -1.2)
p.ICA.EjaC.Ckot <- p


##### PANEL C: SAGUARO #####
cac <- read.table(saguaro_output.summary.file, header = TRUE, stringsAsFactors = FALSE)

## Radiation monophyly:
rad.mon.p <- sum(filter(cac, EjaC.clade.mono == 1)$percentage) # 90.01968 --> EjaC monophyletic as clade
rad.mon <- data.frame(monophyly = c('yes', 'no'), percentage = c(rad.mon.p, 100 - rad.mon.p))
rad.mon$monophyly <- factor(rad.mon$monophyly, levels = c('yes', 'no'))
perc.label <- c(paste0(round(rad.mon$percentage[1], 2), "%"))

p <- ggplot(rad.mon, aes(x = '', y = percentage, fill = monophyly))
p <- p + geom_bar(width = 1, stat = 'identity')
p <- p + coord_polar("y", start = 0)
p <- p + scale_fill_manual(name = "monophyly\nof radiation", values = c('green', 'red'))
p <- p + theme_minimal()
p <- p + theme(axis.text.x = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = 'top',
               legend.title = element_text(size = 14, face = 'bold'),
               legend.text = element_text(size = 14),
               legend.margin = margin(0, 0, 0, 0),
               legend.box.margin = margin(0, 0, -5, 0),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               plot.margin = margin(1, 0, 1, 0, "cm"))
p <- ggarrange(p, ncol = 1, nrow = 1)
p <- p + draw_plot_label(label = perc.label, hjust = 0,
                         size = 16, x = c(0.45), y = c(0.4))
p.rad <- p


## Species monophyly:
sp.mon.p <- sum(filter(cac, EjaC.species.mono == 1)$percentage) # 87.61261 --> EjaC species monophyletic
sp.mon <- data.frame(monophyly = c('yes', 'no'), percentage = c(sp.mon.p, 100 - sp.mon.p))
sp.mon$monophyly <- factor(sp.mon$monophyly, levels = c('yes', 'no'))
perc.label <- c(paste0(round(sp.mon$percentage[1], 2), "%"))

p <- ggplot(sp.mon, aes(x = '', y = percentage, fill = monophyly))
p <- p + geom_bar(width = 1, stat = 'identity')
p <- p + coord_polar("y", start = 0)
p <- p + scale_fill_manual(name = "monophyly\nof each species", values = c('green', 'red'))
p <- p + theme_minimal()
p <- p + theme(axis.text.x = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = 'top',
               legend.title = element_text(size = 14, face = 'bold'),
               legend.text = element_text(size = 14),
               legend.margin = margin(0, 0, 0, 0),
               legend.box.margin = margin(0, 0, -5, 0),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               plot.margin = margin(1, 0, 1, 0, "cm"))
p <- ggarrange(p, ncol = 1, nrow = 1)
p <- p + draw_plot_label(label = perc.label, hjust = 0,
                         size = 16, x = c(0.45), y = c(0.4))
p.sp <- p


## Tree example, radiation monophyly:
cac <- read.nexus(saguaro_cactus1.file)
pops <- gsub('Sgal', 'Cmam', gsub('Tgui', 'Cgui', substr(cac$tip.label, 1, 4)))
d <- data.frame(node = 1:(length(cac$tip.label) + Nnode(cac)),
                color = c(sapply(pops, getpopcol), 'black', 'black', 'black', 'black',
                          getpopcol('Cgui'), getpopcol('Cfus'), getpopcol('Cfus'), getpopcol('Cdec')))
p <- ggtree(cac, layout = 'unrooted', size = 1)  %<+% d + aes(color = I(color))
p <- p + theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))
p <- ggarrange(p, ncol = 1, nrow = 1)
p <- p + draw_plot_label(label = "Example:\ncactus 1 (2.49%)",
                         hjust = 0, x = 0.05, y = 1, size = 16)
p.c1 <- p


## Tree example, species monophyly:
cac <- read.nexus(saguaro_cactus3.file)
pops <- gsub('Sgal', 'Cmam', gsub('Tgui', 'Cgui', substr(cac$tip.label, 1, 4)))
d <- data.frame(node = 1:(length(cac$tip.label) + Nnode(cac)),
                color = c(sapply(pops, getpopcol), 'black', getpopcol('Cdec'), 'black', getpopcol('Cfus'),
                          getpopcol('Cfus'), 'black', getpopcol('Cgui'), getpopcol('Ceja')))
p <- ggtree(cac, layout = 'unrooted', size = 1)  %<+% d + aes(color = I(color))
p <- p + theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))
p <- ggarrange(p, ncol = 1, nrow = 1)
p <- p + draw_plot_label(label = "Example:\ncactus 3 (35.63%)",
                         hjust = 0, x = 0.05, y = 1, size = 16)
p.c3 <- p


## Combine plots:
(plots <- ggarrange(p.rad, p.c1, p.sp, p.c3, ncol = 2, nrow = 2))
plots <- plots + draw_plot_label(label = "C",
                                  x = 0.06, y = 0.99, hjust = 0, size = 26)
ggexport(plots, filename = fig2.file, width = 600, height = 600)

