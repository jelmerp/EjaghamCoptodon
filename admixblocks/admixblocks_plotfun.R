##### SET-UP #####

## Files to read:
fd.file <- 'analyses_output/fd/fd_EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01.txt'
fd.triplets.file <- 'analyses/admixblocks/fd/input/trios_fd.txt'
popcols.file <- 'metadata/popcolours.txt'
scaffolds.file <- 'metadata/scaffolds_withLength.txt'

## Libraries:
library(gdata)
library(reshape2)
library(RColorBrewer)
library(data.table)
library(tidyverse)

## Pop colours:
getpopcol <- function(pop, fcol = 'pop.gphocsfile1') inds$pop.col[match(pop, inds[, fcol])]
inds <- read.table(popcols.file, header = TRUE, stringsAsFactors = FALSE)

## Triplets:
triplets <- read.table(fd.triplets.file, header= TRUE, as.is = TRUE)
triplets.long <- paste0(triplets$popA, '.', triplets$popB, '.', triplets$popC)
triplets.short <- gsub('\\.', '', triplets.long)
triplets.short <- gsub('Cfus', 'F', triplets.short)
triplets.short <- gsub('Cdec', 'D', triplets.short)
triplets.short <- gsub('Ceja', 'E', triplets.short)
triplets.short <- gsub('Cmam', 'A', triplets.short)
triplets.short <- gsub('Cgui', 'U', triplets.short)

## Scaffolds:
scaffolds <- read.table(scaffolds.file, header = TRUE, as.is = TRUE)
scaffolds.gt1Mb <- scaffolds[scaffolds$scaffold.length > 1000000, ]
scaffolds.NC <- scaffolds[grep('NC_', scaffolds$scaffold.name), ]
scaffolds.NC <- scaffolds.NC[-which(scaffolds.NC$scaffold.name == 'NC_013663.1'), ]$scaffold.name

## Load data:
fd.file <- as.data.frame(fread(fd.filename))

## Misc:
popnames.long <- c('Cdec', 'Ceja', 'Cfus', 'Cgui', 'Cmam')
popnames.short <- c('D', 'E', 'F', 'U', 'A')


##### MANHATTAN PLOT #####
## Function for a single Manhattan plot using ggplot:
ggman <- function(df, xvar = 'start', yvar = 'fd', colvar.lines = NULL, colvar.points = NULL,
                  scaffold = NULL, pop = NULL,
                  drawpoints = TRUE, cols.points = 'brew',
                  drawlines = FALSE, cols.lines = 'brew', smoothmethod = 'loess', smoothpar = 0.2,
                  my.ymin = NULL, my.ymax = NULL, my.xlims = NULL,
                  blocklims = NULL, hline = NULL,
                  xlab = NULL, ylab = NULL, legplot = TRUE, legpos = 'top', legtextsize = 15,
                  plot.title = NULL, plot.xann = TRUE, plot.xtitle = FALSE,
                  plotfile = NULL, plot.folder = 'figures/', filetype = 'png',
                  printplot = TRUE, saveplot = TRUE, file.open = FALSE) {

  df <- as.data.frame(df)

  if(xvar == 'start') df$start <- df$start / 1000000

  if(is.null(scaffold)) {
    scaffold <- df$scaffold[1]
    cat('No scaffold specified, plotting the first scaffold:', scaffold, '\n')
  }

  if(scaffold != 'all') df <- df[df$scaffold == scaffold, ]
  if(scaffold == 'all' & xvar == 'start') {
    df$window.nr <- seq(1, nrow(df))
    xvar <- 'window.nr'
  }

  if(!is.null(pop)) df <- df[df$pop %in% pop, ]
  if(is.null(pop)) pop <- df$pop[1]

  if(is.null(plot.title)) plot.title <- paste(scaffold, paste(pop, collapse = '_'))

  popnames <- unique(df$pop)
  if(cols.points[1] == 'popcols') cols.points <- inds$pop.col[match(popnames, inds$pop.gphocsfile1)]
  if(cols.lines[1] == 'popcols') cols.lines <- inds$pop.col[match(popnames, inds$pop.gphocsfile1)]
  if(cols.points[1] == 'brew') cols.points <- brewer.pal(7, "Set1")
  if(cols.lines[1] == 'brew') cols.lines <- brewer.pal(7, "Set1")

  if(!is.null(my.xlims)) {
    if(min(df$start) < my.xlims[1]) df <- df[-which(df$start < my.xlims[1]), ]
    if(max(df$start) > my.xlims[2]) df <- df[-which(df$start > my.xlims[2]), ]
  }

  if(is.null(my.ymax)) {
    my.ymax <- max(df[, yvar], na.rm = TRUE) + (0.2 * max(df[, yvar], na.rm = TRUE))
    if(my.ymax > 1) my.ymax <- 1
  }
  if(is.null(my.ymin)) my.ymin <- min(df[, yvar], na.rm = TRUE) - abs((0.2 * min(df[, yvar], na.rm = TRUE)))


  p <- ggplot()

  if(!is.null(colvar.points)) {
    if(drawpoints == TRUE) {
      p <- p + geom_point(aes_string(x = xvar, y = yvar, fill = colvar.points), data = df,
                          colour = 'white', shape = 21, size = 2)
      p <- p + scale_fill_manual(values = cols.points[1:length(unique(df[, colvar.points]))])
    }
    if(drawpoints == 'blocksOnly') {
      #cat('Drawing points for blocks only...\n')
      p <- p + geom_point(aes_string(x = xvar, y = yvar, fill = colvar.points), data = subset(df, is.block == TRUE),
                          colour = 'white', shape = 21, size = 2)
      p <- p + scale_fill_manual(values = cols.points[1:length(unique(df[, colvar.points]))])
      p <- p + guides(fill = FALSE)
    }
  }

  if(is.null(colvar.points)) {
    if(drawpoints == TRUE) p <- p + geom_point(aes_string(x = xvar, y = yvar), data = df, colour = cols.points[1])
    if(drawpoints == 'blocksOnly') p <- p + geom_point(aes_string(x = xvar, y = yvar), data = subset(df, is.block == TRUE),
                                                       colour = cols.points[1])
  }

  if(!is.null(colvar.lines)) {
    if(drawlines == TRUE) {
      p <- p + geom_smooth(aes_string(x = xvar, y = yvar, colour = colvar.lines), data = df,
                           se = FALSE, inherit.aes = TRUE, size = 2, method = smoothmethod, span = smoothpar)
      if(yvar != 'fd') p <- p + scale_color_manual(values = cols.lines[1:length(unique(df[, colvar.lines]))])
      if(yvar == 'fd') {
        rawvals <- unique(df[, colvar.lines])

        popA1 <- sapply(strsplit(rawvals, split = '\\.'), '[', 1)[1]
        popA2 <- sapply(strsplit(rawvals, split = '\\.'), '[', 1)[2]
        popA3 <- sapply(strsplit(rawvals, split = '\\.'), '[', 1)[3]
        popB1 <- sapply(strsplit(rawvals, split = '\\.'), '[', 2)[1]
        popB2 <- sapply(strsplit(rawvals, split = '\\.'), '[', 2)[2]
        popB3 <- sapply(strsplit(rawvals, split = '\\.'), '[', 2)[3]
        popC1 <- sapply(strsplit(rawvals, split = '\\.'), '[', 3)[1]
        popC2 <- sapply(strsplit(rawvals, split = '\\.'), '[', 3)[2]
        popC3 <- sapply(strsplit(rawvals, split = '\\.'), '[', 3)[3]

        p <- p + scale_color_manual(values = cols.lines[1:length(unique(df[, colvar.lines]))],
                                    breaks = unique(df[, colvar.lines]),
                                    labels = c(bquote(paste(.(popA1), "-", bold(.(popB1)), "-", bold(.(popC1)))),
                                               bquote(paste(.(popA2), "-", bold(.(popB2)), "-", bold(.(popC2)))),
                                               bquote(paste(.(popA3), "-", bold(.(popB3)), "-", bold(.(popC3))))))
      }
    }
  }

  if(is.null(colvar.lines)) if(drawlines == TRUE)
    p <- p + geom_smooth(aes_string(x = xvar, y = yvar, colour = colvar.lines), data = df,
                         se = FALSE, fill = "NA", inherit.aes = FALSE, size = 2, method = smoothmethod, span = smoothpar)

  if(!is.null(blocklims)) {
    p <- p + geom_vline(xintercept = blocklims[1], colour = 'grey10', linetype = 2)
    p <- p + geom_vline(xintercept = blocklims[2], colour = 'grey10', linetype = 2)
  }
  if(!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, colour = 'grey10', linetype = 1)
  }


  if(yvar == 'fst') ylab <- expression(F[ST])
  if(yvar == 'fd') ylab <- expression(f[d])

  p <- p + labs(title = plot.title)
  if(!is.null(xlab)) p <- p + labs(x = xlab)
  if(!is.null(ylab)) p <- p + labs(y = ylab)

  p <- p + theme_bw()

  p <- p + theme(axis.text.y = element_text(size = 16))
  if(plot.xann == TRUE) p <- p + theme(axis.text.x = element_text(size = 16))
  if(plot.xann == FALSE) p <- p + theme(axis.text.x = element_blank())

  if(is.null(xlab)) p <- p + theme(axis.title.x = element_blank())
  if(!is.null(xlab)) p <- p + theme(axis.title.x = element_text(size = 18))
  p <- p + theme(axis.title.y = element_text(size = 18))

  if(plot.title != FALSE) p <- p + theme(plot.title = element_text(size = 20, hjust = 0.5))
  if(plot.title == FALSE) p <- p + theme(plot.title = element_blank())

  p <- p + theme(plot.margin = unit(c(1, 1.2, 1, 0.5), "lines")) # top, right, ..

  if(my.ymax < 0.5) my.ticks <- c(0.1, 0.2, 0.3, 0.4)
  if(my.ymax > 0.5) my.ticks <- c(0.2, 0.4, 0.6, 0.8)

  p <- p + scale_y_continuous(expand = c(0, 0), limits = c(my.ymin, my.ymax), breaks = my.ticks)
  if(!is.null(my.xlims)) p <- p + scale_x_continuous(expand = c(0, 0), limits = my.xlims)
  if(is.null(my.xlims)) p <- p + scale_x_continuous(expand = c(0, 0))

  if(legplot == TRUE) {
    p <- p + theme(legend.title = element_text(size = legtextsize, face = 'bold'))
    p <- p + theme(legend.title = element_blank())
    p <- p + theme(legend.text = element_text(size = legtextsize))
    p <- p + theme(legend.key.height = unit(0.5, "cm"))
    p <- p + theme(legend.key.width = unit(0.5, "cm"))
    p <- p + theme(legend.position = legpos)
    p <- p + theme(legend.margin = margin(0, 0, 0, 0))
    p <- p + theme(legend.box.margin = margin(-5, -5, -5, -5))
  }
  if(legplot == FALSE) p <- p + theme(legend.position = "none")

  pop.filename <- paste0(pop, collapse = '_')
  if(is.null(plotfile)) plotfile <- paste0(scaffold, '.', yvar, '.', pop.filename)
  plotfile.full <- paste0(plot.folder, '/', plotfile, '.', filetype)

  if(saveplot == TRUE) {
    cat('saving:', plotfile.full, '\n')
    ggsave(filename = plotfile.full, plot = p, width = 10, height = 8)
  }
  if(file.open == TRUE) system(paste("xdg-open", plotfile.full))

  if(printplot == TRUE) print(p)
  return(p)
}
