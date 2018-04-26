##### SETUP #####

## Files to read:
popcols.file <- 'metadata/popcolours.txt'
raxml.withSgal.file <- 'analyses_output/RAxML_IC_Score_BranchLabels.EjaC.Dstat.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.ICAscores.MLtree'
raxml.withSgalCkot.file <- 'analyses_output/RAxML_IC_Score_BranchLabels.EjaC.Ckot.Dstat.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.ICAscores.MLtree'

## Files to write:
raxmltrees.suppfig.file <- 'figures/SuppFig_raxmlTrees.png'

## Libraries:
library(ape)
library(tidyverse)
library(ggtree)
library(gdata)

## Read general files:
inds <- read.table(popcols.file, header = TRUE, stringsAsFactors = FALSE)

## Functions:
getpopcol <- function(pop, fcol = 'pop.gphocsfile1') inds$pop.col[match(pop, inds[, fcol])]


##### PANEL A: WITH SGAL #####
raxml <- read.raxml(raxml.withSgal.file)
raxml@phylo$edge.length[2] <- 0.04 # Shorten very long S. galilaeus branch
raxml@data$bootstrap[which(raxml@data$bootstrap == 200)] <- -2 # Had to manually change this in file, otherwise not accepted by ggtree
raxml@data$bootstrap <- raxml@data$bootstrap / 100
p <- ggtree(raxml, size = 1)
p <- p + geom_label2(aes(label = bootstrap, fill = bootstrap, subset=(node!=c(20))), size = 3)
p <- p + scale_fill_continuous(low = 'red', high = 'darkgreen', name = 'ICA score')
p <- p + ggplot2::xlim(-0.005, 0.1)
p <- p + geom_cladelabel(node = 13, align = TRUE, color = 'black', fontsize = 5, offset = 0.003,
                         label = 'paste(italic("S. galilaeus"))', parse = TRUE)
p <- p + geom_cladelabel(node = 10, align = TRUE, color = 'orangered1', fontsize = 5, offset = 0.003,
                         label = 'paste(italic("C. sp. Mamfé"))', parse = TRUE)
p <- p + geom_cladelabel(node = 24, align = TRUE, color = 'orange2', fontsize = 5, offset = 0.003,
                         label = 'paste(italic("C. guineensis"))', parse = TRUE)
p <- p + geom_cladelabel(node = 19, align = TRUE, color = 'green3', fontsize = 5, offset = 0.003,
                         label = 'paste(italic("C. fusiforme"))', parse = TRUE)
p <- p + geom_cladelabel(node = 23, align = TRUE, color = 'turquoise2', fontsize = 5, offset = 0.003,
                         label = 'paste(italic("C. deckerti"))', parse = TRUE)
p <- p + geom_cladelabel(node = 22, align = TRUE, color = 'blue2', fontsize = 5, offset = 0.003,
                         label = 'paste(italic("C. ejagham"))', parse = TRUE)
p <- p + annotate(geom = "text", x = -0.005, y = 13.5, label = "A", fontface = 'bold', color = 'black', size = 8)
trueboot <- subset(p$data, isTip == FALSE)
trueboot$truebootstrap1 <- c(NA, NA, '*', '*', '*', '*', NA, '*', '*', '*', '*')
p <- p + geom_text(data = trueboot, aes(label = trueboot$truebootstrap1), size = 7, hjust = 1.2, vjust = -0.05)
SgalOnly.plot <- p


##### PANEL B: WITH SGAL & CKOT #####
raxml <- read.raxml(raxml.withSgalCkot.file)
raxml@phylo$edge.length[2] <- 0.06
raxml@data$bootstrap[which(raxml@data$bootstrap == 200)] <- -2
raxml@data$bootstrap <- raxml@data$bootstrap / 100
p2 <- ggtree(raxml, size = 1)
p2 <- p2 + geom_label2(aes(label = bootstrap, fill = bootstrap, subset=(node!=c(23))), size = 3)
p2 <- p2 + scale_fill_continuous(low = 'red', high = 'darkgreen', name = 'ICA score')
p2 <- p2 + ggplot2::xlim(-0.01, 0.15)
p2 <- p2 + geom_cladelabel(node = 15, align = TRUE, color = 'black', fontsize = 5, offset = 0.005,
                         label = 'paste(italic("S. galilaeus"))', parse = TRUE)
p2 <- p2 + geom_cladelabel(node = 28, align = TRUE, color = 'black', fontsize = 5, offset = 0.005,
                         label = 'paste(italic("C. kottae"))', parse = TRUE)
p2 <- p2 + geom_cladelabel(node = 8, align = TRUE, color = 'orangered1', fontsize = 5, offset = 0.005,
                         label = 'paste(italic("C. sp. Mamfé"))', parse = TRUE)
p2 <- p2 + geom_cladelabel(node = 27, align = TRUE, color = 'orange2', fontsize = 5, offset = 0.005,
                         label = 'paste(italic("C. guineensis"))', parse = TRUE)
p2 <- p2 + geom_cladelabel(node = 22, align = TRUE, color = 'green3', fontsize = 5, offset = 0.005,
                         label = 'paste(italic("C. fusiforme"))', parse = TRUE)
p2 <- p2 + geom_cladelabel(node = 26, align = TRUE, color = 'turquoise2', fontsize = 5, offset = 0.005,
                         label = 'paste(italic("C. deckerti"))', parse = TRUE)
p2 <- p2 + geom_cladelabel(node = 25, align = TRUE, color = 'blue2', fontsize = 5, offset = 0.005,
                         label = 'paste(italic("C. ejagham"))', parse = TRUE)
p2 <- p2 + annotate(geom = "text", x = -0.01, y = 15.5, label = "B", fontface = 'bold', color = 'black', size = 8)
trueboot2 <- subset(p2$data, isTip == FALSE)
trueboot2$truebootstrap1 <- c(NA, NA, '*', '*', NA, '*', '*', NA, '*', '*', '*','*','*')
trueboot2$truebootstrap2 <- c(NA, NA, NA, NA, 88, NA, NA, NA, NA, NA, NA, NA, NA)
p2 <- p2 + geom_text(data = trueboot2, aes(label = trueboot2$truebootstrap1), size = 7, hjust = 1.2, vjust = -0.05)
p2 <- p2 + geom_text(data = trueboot2, aes(label = trueboot2$truebootstrap2), size = 3.5, hjust = 1.2, vjust = -1.2)
SgalCkot.plot <- p2


##### COMBINE PLOTS #####
plots <- ggarrange(SgalOnly.plot, SgalCkot.plot, ncol = 2, nrow = 1)
ggexport(plots, filename = raxmltrees.suppfig.file, width = 750, height = 400)