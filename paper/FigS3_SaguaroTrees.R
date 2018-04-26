##### SETUP #####

## Files to read:
popcols.file <- 'metadata/popcolours.txt'
cactus3.file <- 'analyses_output/saguaro_cactus3.tre'
cactus20.file <- 'analyses_output/saguaro_cactus20.tre'

## Files to write:
cactus3.figure <- 'figures/SuppFig_cactus3.png'
cactus20.figure <- 'figures/SuppFig_cactus20.png'

## Libraries:
library(ape)
library(tidyverse)
library(ggtree)
library(gdata)

## Read general files:
inds <- read.table(popcols.file, header = TRUE, stringsAsFactors = FALSE)

## Functions:
getpopcol <- function(pop, fcol = 'pop.gphocsfile1') inds$pop.col[match(pop, inds[, fcol])]


##### CACTUS 3 #####
c3 <- read.nexus(cactus3.file)
pops <- gsub('Sgal', 'Cmam', gsub('Tgui', 'Cgui', substr(c3$tip.label, 1, 4)))
d <- data.frame(node = 1:(length(c3$tip.label) + Nnode(c3)),
                color = c(sapply(pops, getpopcol), 'black', getpopcol('Cdec'), 'black', getpopcol('Cfus'),
                          getpopcol('Cfus'), 'black', getpopcol('Cgui'), getpopcol('Ceja')))

p <- ggtree(c3, layout = 'unrooted', size = 1)  %<+% d + aes(color=I(color))
p <- p + geom_tiplab(aes(subset=(node==17)), label = 'paste(italic("C."), italic(" guineensis"))', vjust = -2, offset = 0, parse = TRUE, size = 5, color = 'orange2')
p <- p + geom_tiplab(aes(subset=(node==8)), label = 'paste(italic("C."), italic(" sp. Mamfé"))', offset = 0, parse = TRUE, size = 5, color = 'orangered1')
p <- p + geom_tiplab(aes(subset=(node==6)), label = 'paste(italic("C."), italic(" fusiforme"))', hjust = -0.15, vjust = 0, parse = TRUE, size = 5, color = 'green3')
p <- p + geom_tiplab(aes(subset=(node==1)), label = 'paste(italic("C."), italic(" deckerti"))', vjust = -0.5, hjust = 0.8, parse = TRUE, size = 5, color = 'turquoise2')
p <- p + geom_tiplab(aes(subset=(node==4)), label = 'paste(italic("C."), italic(" ejagham"))', vjust = -2, hjust = 0.85, parse = TRUE, size = 5, color = 'blue2')
p <- p + ggplot2::xlim(-0.1, 0.4)
p <- p + theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))
p <- p + annotate(geom = "text", x = 0.02, y = 0.04, label = "Cactus 3 (35.63%)", color = 'black', size = 6.5)
ggsave(cactus3.figure, p, width = 7, height = 5)


##### CACTUS 20 #####
c20 <- read.nexus(cactus20.file)

pops <- gsub('Sgal', 'Cmam', gsub('Tgui', 'Cgui', substr(c20$tip.label, 1, 4)))
d <- data.frame(node = 1:(length(c20$tip.label) + Nnode(c20)),
                color = c(sapply(pops, getpopcol), 'black', getpopcol('Cfus'), getpopcol('Cfus'), 'black',
                          getpopcol('Cdec'), getpopcol('Ceja'), 'black', getpopcol('Cgui')))
p <- ggtree(c20, layout = 'unrooted', size = 1)  %<+% d + aes(color=I(color))

p <- p + geom_tiplab(aes(subset=(node==18)), label = 'paste(italic("C."), italic(" guineensis"))', hjust = -0.1, vjust = 1, parse = TRUE, size = 5, color = 'orange2')
p <- p + geom_tiplab(aes(subset=(node==8)), label = 'paste(italic("C."), italic(" sp. Mamfé"))', hjust = 0.5, vjust = -0.5, parse = TRUE, size = 5, color = 'orangered1')
p <- p + geom_tiplab(aes(subset=(node==6)), label = 'paste(italic("C."), italic(" fusiforme"))', hjust = 0, vjust = 2, parse = TRUE, size = 5, color = 'green3')
p <- p + geom_tiplab(aes(subset=(node==2)), label = 'paste(italic("C."), italic(" deckerti"))', hjust = 1.2, vjust = -0.7, parse = TRUE, size = 5, color = 'turquoise2')
p <- p + geom_tiplab(aes(subset=(node==3)), label = 'paste(italic("C."), italic(" ejagham"))', hjust = 0, vjust = -0.5, parse = TRUE, size = 5, color = 'blue2')
p <- p + ggplot2::xlim(-0.6, 0.3)
p <- p + theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))
p <- p + annotate(geom = "text", x = -0.47, y = 0.14, label = "Cactus 20 (49.38%)", color = 'black', size = 6.5)
ggsave(cactus20.figure, p, width = 7, height = 5)
