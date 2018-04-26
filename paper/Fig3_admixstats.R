##### SET-UP #####

## Files to read:
inds <- read.table(popcols.file, header = TRUE, stringsAsFactors = FALSE)

## Admixtools and dfoil output:
dstats.file <- 'analyses_output/admixtools_EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01.EjaC.Dstat.dmode.out'
f4ratio.file <- 'analyses_output/admixtools_EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01.f4ratio.out'
dfoil.noCdec.file <- 'analyses_output/dfoil_EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01.noCdec'
dfoil.noCeja.file <- 'analyses_output/dfoil_EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01.noCeja'
dfoil.noCfus.file <- 'analyses_output/dfoil_EjaC.Dstat.DP5.GQ20.MAXMISS0.5.MAF0.01.noCfus'

## Trees illustrating the method:
d1.tree.file <- 'metadata/d1.tree.png'
d2.tree.file <- 'metadata/d2.tree.png'
f4r.tree.file <- 'metadata/f4r.tree.eps'
dfoil.tree.file <- 'metadata/dfolil.tree.png'

## Files to write:
fig3.file <- 'figures/Fig3_admixstats.png'

## Libraries:
library(gdata)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(png)
library(grid)


##### FUNCTIONS #####
getpopcol <- function(pop, fcol = 'pop.gphocsfile1') inds$pop.col[match(pop, inds[, fcol])]

return.dfmode <- function(running.mode = 'dmode', filename) {
  output <- read.delim(filename, sep = "", header = FALSE, as.is = TRUE)
  output <- output[, -1] # Get rid of "result:" column
  if(running.mode == 'dmode') colnames(output) <- c('popA', 'popB', 'popC', 'popD', 'D', 'se', 'Z', 'BABA', 'ABBA', 'nrSNPS')
  if(running.mode == 'fmode') colnames(output) <- c('popA', 'popB', 'popC', 'popD', 'f4', 'se', 'Z', 'BABA', 'ABBA', 'nrSNPS')
  return(output)
}

return.f4ratio <- function(filename) {
  output <- read.delim(outputfile, sep = "", header = FALSE)
  output <- output[, -c(1, 6:8, 10)]
  colnames(output) <- c('popA', 'popD', 'popX', 'popC', 'popB', 'alpha', 'se', 'Z')
  output <- output[, c('popA', 'popB', 'popX', 'popC', 'popD', 'alpha', 'se', 'Z')]
  output$Z2 <- abs(1 - output$alpha) / output$se
  output <- output %>% mutate(alpha = round(alpha, 3), se = round(se, 5), Z = round(Z, 2), Z2 = round(Z2, 2))
  return(output)
}

read.dfoil.out <- function(pop, filename) {
  dfoil <- read.delim(filename, header = TRUE, as.is = TRUE)
  dfoil$coord <- NULL
  dfoil <- dplyr::rename(dfoil, pop = X.chrom)
  dfoil$pop <- pop
  dfoil <- mutate(dfoil, T12 = round(T12, 4), T34 = round(T34, 4), T1234 = round(T1234, 4),
                  DFO = round(DFO_stat, 3), DFO_chisq = round(DFO_chisq),
                  DIL = round(DIL_stat, 3), DIL_chisq = round(DIL_chisq),
                  DFI = round(DFI_stat, 3), DFI_chisq = round(DFI_chisq),
                  DOL = round(DOL_stat, 3), DOL_chisq = round(DOL_chisq))
  return(dfoil)
}


##### PANEL A: D-STATISTICS #####
## Prep data:
d <- return.dfmode(dstats.file)
popcombs <- c('Cdec.Ceja.Cmam', 'Cfus.Cdec.Cmam', 'Cfus.Ceja.Cmam', 'Cgui.Cmam.Cdec', 'Cgui.Cmam.Ceja', 'Cgui.Cmam.Cfus')
d$popcomb <- paste0(d$popA, '.', d$popB, '.', d$popC)
d$riverpop <- ifelse(d$popA %in% c('Cmam', 'Cgui'), d$popB, d$popC)
d$lakepop <- ifelse(d$popA %in% c('Cdec', 'Ceja', 'Cfus'), d$popB, d$popC)
d <- filter(d, popcomb %in% popcombs) %>% arrange(popcomb)

## Plot:
p <- ggplot(d, aes(x = factor(popcomb), y = -D, colour = popcomb))
p <- p + geom_pointrange(aes(ymax = -d$D + d$se, ymin = -d$D - d$se))
p <- p + geom_hline(yintercept = 0, colour = 'grey40')
p <- p + coord_flip()
p <- p + scale_colour_manual(values = as.character(sapply(d$lakepop, getpopcol)))
p <- p + scale_x_discrete(breaks = popcombs,
                          labels = c(expression(paste("Dec-", "Eja", "-", "Mam")),
                                     expression(paste("Fus-", bold("Dec"), "-", bold("Mam"))),
                                     expression(paste("Fus-", bold("Eja"), "-", bold("Mam"))),
                                     expression(paste("Gui-", bold("Mam"), "-", bold("Dec"))),
                                     expression(paste("Gui-", bold("Mam"), "-", bold("Eja"))),
                                     expression(paste("Gui-", bold("Mam"), "-", bold("Fus")))))
p <- p + labs(y = "D")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 16),
               axis.text.y = element_text(size = 16),
               axis.title.x = element_text(size = 18),
               axis.title.y = element_blank(),
               plot.title = element_text(size = 18, hjust = 2),
               plot.margin = margin(1, 1, 1.2, 0.8, "cm"))
p <- p + guides(colour = FALSE)
d.plot <- p


##### PANEL B: F4-RATIO #####
f4r <- return.f4ratio(f4ratio.file)
f4r <- f4r %>% filter(popC == 'Cmam', popA == 'Cfus', popB == 'Cdec')
f4r$pop <- f4r$popB
f4r2 <- f4r
f4r2$alpha <- 1 - f4r2$alpha
f4r2$pop <- 'Cmam'
f4r <- rbind(f4r, f4r2)
f4r$popcomb <- paste0(f4r$popB, '.', f4r$popX)
f4r <- f4r %>% arrange(popcomb, pop)
f4r$pop <- factor(f4r$pop, levels = c('Cmam', 'Cdec', 'Ceja'))
f4r$se[2] <- NA

p <- ggplot(f4r, aes(x = factor(popcomb), y = alpha, fill = pop))
p <- p + geom_bar(stat = "identity", position = 'stack', colour = 'black')
p <- p + coord_flip()
p <- p + geom_errorbar(aes(ymax = f4r$alpha + f4r$se,
                           ymin = f4r$alpha - f4r$se), width = 0.25, colour = 'grey20')
p <- p + scale_fill_manual(name = "", values = c('orangered1', 'turquoise2'),
                           labels = c(expression(italic("C. sp. Mamfé  ")), expression(italic('C. deckerti'))))
p <- p + scale_x_discrete(breaks = c("Cdec.Ceja"), labels = c('C. ejagham'), expand = c(0, 0))
p <- p + labs(y = bquote(italic(f[4]) ~ "ratio ancestry proportions"))
p <- p + ggtitle(expression(italic("C. ejagham")))
p <- p + theme_bw()
p <- p + theme(plot.title = element_text(size = 20, hjust = 0.5, vjust = -1),
               axis.title.x = element_text(size = 18),
               axis.title.y = element_blank(),
               axis.text.x = element_text(size = 16),
               axis.text.y = element_blank(),
               panel.border = element_rect(colour = 'grey80'),
               legend.position = 'bottom',
               legend.title = element_text(size = 16, face = 'bold'),
               legend.text = element_text(size = 16),
               legend.key.height = unit(0.5, "cm"),
               legend.key.width = unit(0.5, "cm"),
               plot.margin = margin(0, 1, 0.1, 1.5, "cm"))
f4r.plot <- p


##### PANEL C: DFOIL #####
## Prep data:
pops <- c('noCdec', 'noCeja', 'noCfus')
dfoil.Sgui <- do.call(rbind, lapply(pops, read.dfoil.out, filename = dfoil.file))
df <- melt(dfoil.Sgui, id.vars = c('pop'), measure.vars = c('DFO', 'DIL', 'DFI', 'DOL'))

## Create plot:
p <- ggplot(df, aes(x = factor(pop), y = value, fill = variable))
p <- p + geom_bar(stat = "identity", position = 'dodge', colour = 'black')
p <- p + scale_fill_discrete(name = expression(D[FOIL] ~ statistic),
                             labels = c(expression(D[FO]), expression(D[IL]), expression(D[FI]), expression(D[OL])))
p <- p + scale_x_discrete(labels = c('Eja-Fus', 'Dec-Fus', 'Dec-Eja'))
p <- p + scale_y_continuous(limits = c(-0.25, 0.11))
p <- p + labs(x = "p1-p2", y = expression(value ~ of ~ D[FOIL] ~ statistic))
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16))
p <- p + theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))
p <- p + theme(legend.position = 'top')
p <- p + theme(legend.title = element_text(size = 15, face = 'bold'), legend.text = element_text(size = 15))
p <- p + theme(legend.key.height = unit(0.5, "cm"), legend.key.width = unit(0.5, "cm"))
dd <- ggplot_build(p)$data[[1]]
p <- p + annotate(geom = "text", x = dd$x[1], y = 0.09, label = "0", color = dd$fill[1], size = 6)
p <- p + annotate(geom = "text", x = dd$x[2], y = 0.09, label = "0", color = dd$fill[2], size = 6)
p <- p + annotate(geom = "text", x = dd$x[3], y = 0.09, label = "-", color = dd$fill[3], size = 8)
p <- p + annotate(geom = "text", x = dd$x[4], y = 0.09, label = "-", color = dd$fill[4], size = 8)
p <- p + annotate(geom = "text", x = dd$x[5], y = 0.09, label = "0", color = dd$fill[5], size = 6)
p <- p + annotate(geom = "text", x = dd$x[6], y = 0.09, label = "0", color = dd$fill[6], size = 6)
p <- p + annotate(geom = "text", x = dd$x[7], y = 0.09, label = "-", color = dd$fill[7], size = 8)
p <- p + annotate(geom = "text", x = dd$x[8], y = 0.09, label = "-", color = dd$fill[8], size = 8)
p <- p + annotate(geom = "text", x = dd$x[9], y = 0.09, label = "+", color = dd$fill[9], size = 8)
p <- p + annotate(geom = "text", x = dd$x[10], y = 0.09, label = "+", color = dd$fill[10], size = 8)
p <- p + annotate(geom = "text", x = dd$x[11], y = 0.09, label = "-", color = dd$fill[11], size = 8)
p <- p + annotate(geom = "text", x = dd$x[12], y = 0.09, label = "-", color = dd$fill[12], size = 8)
p <- p + annotate(geom = "text", x = 1, y = 0.11, label = "p12<->p4", color = 'grey20', size = 6)
p <- p + annotate(geom = "text", x = 2, y = 0.11, label = "p12<->p4", color = 'grey20', size = 6)
p <- p + theme(plot.margin = margin(1, 0.1, 0.1, 1.2, "cm"))
dfoil.plot <- p


##### IMPORT TREE FIGURES  #####
d1.tree <- rasterGrob(readPNG(d1.tree.file), interpolate = TRUE)
d2.tree <- rasterGrob(readPNG(d2.tree.file), interpolate = TRUE)
f4r.tree <- rasterGrob(readPNG(f4r.tree.file), interpolate = TRUE)
dfoil.tree <- rasterGrob(readPNG(dfoil.tree.file), interpolate = TRUE)

d1.tree <- ggplot(data.frame()) + geom_blank() + annotation_custom(d1.tree, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
d2.tree <- ggplot(data.frame()) + geom_blank() + annotation_custom(d2.tree, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
f4r.tree <- ggplot(data.frame()) + geom_blank() + annotation_custom(f4r.tree, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
dfoil.tree <- ggplot(data.frame()) + geom_blank() + annotation_custom(dfoil.tree, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)


##### COMBINE PLOTS  #####
plots1 <- ggdraw() +
  draw_plot(d.plot, x = 0, y = 0.35, width = 0.7, height = 0.65) +
  draw_plot(f4r.plot, x = 0, y = 0, width = 0.7, height = 0.35) +
  draw_plot(d1.tree, x = 0.65, y = 0.71, width = 0.35, height = 0.23) +
  draw_plot(d2.tree, x = 0.65, y = 0.47, width = 0.35, height = 0.23) +
  draw_plot(f4r.tree, x = 0.65, y = -0.05, width = 0.42, height = 0.47)
plots2 <- ggdraw() +
  draw_plot(dfoil.plot, x = 0, y = 0.25, width = 1, height = 0.75) +
  draw_plot(dfoil.tree, x = 0.08, y = 0.0, width = 1, height = 0.25)

plots <- ggarrange(plots1, plots2, ncol = 2, nrow = 1, widths = c(1.3, 1))
plots <- plots + draw_plot_label(label = c("A", "C", "B"), size = 24,
                                 x = c(0, 0, 0.55), y = c(1, 0.42, 1))

ggexport(plots, filename = fig3.file, width = 900, height = 650)