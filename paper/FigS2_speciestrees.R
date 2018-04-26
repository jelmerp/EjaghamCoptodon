##### SET-UP #####

## Files to read:
popcols.file <- 'metadata/popcolours.txt'
phylonet.speciestree.file <- 'analyses_output/EjaC.Ckot.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb.ST_MDC.phylonet.output.txt'
phylonet.mltree.file <- 'analyses_output/EjaC.Dstat.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb.NMIG0.ML.phylonet.output.txt'
mpest.tree.file <- 'analyses_output/mpest_EjaC.Ckot.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb.tre'
astral.treefile <- 'analyses_output/astral_EjaC.Cgal.Cgui.DP5.GQ20.MAXMISS0.9.MAF0.01.minsites250.Window100kb_allGeneTrees.tre'

## Files to write:
speciestrees.suppfig.file <- 'figures/SuppFig_speciesTrees.png'

## Libraries:
library(ape)
library(tidyverse)
library(ggtree)
library(gdata)

## Read general files:
inds <- read.table(popcols.file, header = TRUE, stringsAsFactors = FALSE)

## Functions:
getpopcol <- function(pop, fcol = 'pop.gphocsfile1') inds$pop.col[match(pop, inds[, fcol])]


##### FUNCTIONS #####
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


##### PANEL A: PHYLONET SPECIES TREE #####
tree <- read.raxml(phylonet.speciestree.file)
p <- ggtree(tree, size = 1)
p <- p + geom_label2(aes(label = bootstrap, fill = bootstrap, subset=(node!=19)), size = 3)
p <- p + scale_fill_continuous(limits = c(0, 100), low = 'red', high = 'darkgreen', name = 'ICA score')
p <- p + geom_tiplab(aes(subset=(node==1)), label = 'paste(italic("C."), italic(" kottae (/ S. galilaeus)"))', parse = TRUE, size = 5)
p <- p + geom_tiplab(aes(subset=(node==2)), label = 'paste(italic("C."), italic(" guineensis"))', parse = TRUE, size = 5, color = 'orange2')
p <- p + geom_tiplab(aes(subset=(node==3)), label = 'paste(italic("C."), italic(" sp. Mamfé"))', parse = TRUE, size = 5, color = 'orangered1')
p <- p + geom_tiplab(aes(subset=(node==4)), label = 'paste(italic("C."), italic(" fusiforme"))', parse = TRUE, size = 5, color = 'green3')
p <- p + geom_tiplab(aes(subset=(node==6)), label = 'paste(italic("C."), italic(" deckerti"))', parse = TRUE, size = 5, color = 'turquoise2')
p <- p + geom_tiplab(aes(subset=(node==5)), label = 'paste(italic("C."), italic(" ejagham"))', parse = TRUE, size = 5, color = 'blue2')
p <- p + ggplot2::xlim(-250, 10000)
p <- p + annotate(geom = "text", x = -250, y = 6.75, label = "A", color = 'black', fontface = 'bold', size = 8)
p <- p + annotate(geom = "text", x = 4500, y = 6.6, label = "Phylonet species tree (MDC)", color = 'black', size = 6.5)
p.phylonet.speciestree.EjaC.Ckot <- p


##### PANEL B: PHYLONET ML TREE #####
tree <- read.tree(phylonet.mltree.file)
p <- ggtree(tree, size = 1)
p <- p + geom_tiplab(aes(subset=(node==1)), label = 'paste(italic("S."), italic(" galilaeus (/ C. kottae)"))', parse = TRUE, size = 5)
p <- p + geom_tiplab(aes(subset=(node==2)), label = 'paste(italic("C."), italic(" guineensis"))', parse = TRUE, size = 5, color = 'orange2')
p <- p + geom_tiplab(aes(subset=(node==3)), label = 'paste(italic("C."), italic(" sp. Mamfé"))', parse = TRUE, size = 5, color = 'orangered1')
p <- p + geom_tiplab(aes(subset=(node==4)), label = 'paste(italic("C."), italic(" fusiforme"))', parse = TRUE, size = 5, color = 'green3')
p <- p + geom_tiplab(aes(subset=(node==6)), label = 'paste(italic("C."), italic(" deckerti"))', parse = TRUE, size = 5, color = 'turquoise2')
p <- p + geom_tiplab(aes(subset=(node==5)), label = 'paste(italic("C."), italic(" ejagham"))', parse = TRUE, size = 5, color = 'blue2')
p <- p + ggplot2::xlim(-0.3, 9)
p <- p + annotate(geom = "text", x = -0.3, y = 6.75, label = "B", color = 'black', fontface = 'bold', size = 8)
p <- p + annotate(geom = "text", x = 4.5, y = 6.6, label = "Phylonet ML network (no reticulations)", color = 'black', size = 6.5)
p.phylonetML.Cgal <- p


##### PANEL C: MPEST #####
trees.mpest <- read.tree(mpest.tree.file)
tree <- trees.mpest[[1]]
tree$edge.length[which(tree$edge.length == 9)] <- 1
p <- ggtree(tree, size = 1)
p <- p + geom_tiplab(aes(subset=(node==6)), label = 'paste(italic("C."), italic(" kottae"))', parse = TRUE, size = 5)
p <- p + geom_tiplab(aes(subset=(node==1)), label = 'paste(italic("C."), italic(" guineensis"))', parse = TRUE, size = 5, color = 'orange2')
p <- p + geom_tiplab(aes(subset=(node==2)), label = 'paste(italic("C."), italic(" sp. Mamfé"))', parse = TRUE, size = 5, color = 'orangered1')
p <- p + geom_tiplab(aes(subset=(node==3)), label = 'paste(italic("C."), italic(" fusiforme"))', parse = TRUE, size = 5, color = 'green3')
p <- p + geom_tiplab(aes(subset=(node==4)), label = 'paste(italic("C."), italic(" deckerti"))', parse = TRUE, size = 5, color = 'turquoise2')
p <- p + geom_tiplab(aes(subset=(node==5)), label = 'paste(italic("C."), italic(" ejagham"))', parse = TRUE, size = 5, color = 'blue2')
p <- p + ggplot2::xlim(-0.2, 8)
p <- p + annotate(geom = "text", x = -0.2, y = 6.75, label = "C", color = 'black', fontface = 'bold', size = 8)
p <- p + annotate(geom = "text", x = 3.5, y = 6.5, label = "MP-EST species tree", color = 'black', size = 6.5)
p.mpest <- p


##### PANEL D: ASTRAL #####
tree.astral <- read.tree(astral.treefile)
tree.astral$edge.length[-which(tree.astral$edge.length == 0)] <- 1
p <- ggtree(tree.astral, layout = 'unrooted', size = 1)
p <- p + geom_tiplab(aes(subset=(node==4)), label = 'paste(italic("C."), italic(" guineensis"))', offset = 0.07, parse = TRUE, size = 5, color = 'orange2')
p <- p + geom_tiplab(aes(subset=(node==5)), label = 'paste(italic("C."), italic(" sp. Mamfé"))', offset = 0.03, parse = TRUE, size = 5, color = 'orangered1')
p <- p + geom_tiplab(aes(subset=(node==3)), label = 'paste(italic("C."), italic(" fusiforme"))', vjust = -0.5, parse = TRUE, size = 5, color = 'green3')
p <- p + geom_tiplab(aes(subset=(node==1)), label = 'paste(italic("C."), italic(" deckerti"))', offset = 0.03, parse = TRUE, size = 5, color = 'turquoise2')
p <- p + geom_tiplab(aes(subset=(node==2)), label = 'paste(italic("C."), italic(" ejagham"))', offset = 0.03, parse = TRUE, size = 5, color = 'blue2')
p <- p + ggplot2::xlim(-1.8, 1.4)
p <- p + annotate(geom = "text", x = -1.8, y = 1.85, label = "D", color = 'black', fontface = 'bold', size = 8)
p <- p + annotate(geom = "text", x = -0.2, y = 1.7, label = "ASTRAL species tree", color = 'black', size = 6.5)
p.astral <- p

##### SUPPLEMENTARY SPECIES TREE FIGURE #####
png(speciestrees.suppfig.file, width = 900, height = 750)
multiplot(p.phylonet.speciestree.EjaC.Ckot, p.mpest, p.phylonetML.Cgal, p.astral, cols = 2)
dev.off()
