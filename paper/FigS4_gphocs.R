###### SET-UP #####

## Other scripts:
source('scripts/gphocs/gphocs_analyze_fun.R')

## Files to write:
gphocs.suppfig.file <- 'figures/SuppFig_Gphocs.png'

## Libraries:
library(ggpubr)
library(cowplot)


##### PANEL A: 2Nm FROM GUINEENSIS ######
Logsel <- filter(Log, var == 'm.prop', migfrom == 'Gui', migto != 'Gui' & migto != 'Mam',
                 migtoType3 %in% c('single', 'all'))
p.2Nm.U <- vplot(Logsel,
                 xvar = 'migto', fillvar = 'cn', colvar = 'migtoType3', shapevar = 'cn', yvar = 'val',
                 linecols = NULL, rm.violins = TRUE, hpdline.width = 1, y.max = 'max.hpd',
                 legpos = 'top', rm.leg.fill = TRUE,
                 xlab = 'migration target', ylab = 'population migration rate (2Nm)',
                 saveplot = FALSE)
p.2Nm.U <- p.2Nm.U + scale_colour_manual(values = brewer.pal(4, name = 'Set1')[c(2,4)],
                                         name = "migration in run\n(from Gui to:)")


##### PANEL A: 2Nm WITHIN-RADIATION ######
Logsel <- filter(Log, var == '2Nm' & migto != 'Gui' & migto != 'Mam',
                 migfromRun %in% c('DE', 'Dec', 'Eja', 'Fus'))
p.2Nm.W <- vplot(Logsel,
                 xvar = 'migto', fillvar = 'cn', colvar = 'migfrom', yvar = 'val',
                 linecols = 'pop.cols', y.max = 'max.hpd', rm.violins = TRUE, hpdline.width = 1,
                 legpos = 'top', legcolname = "migration source:",
                 xlab = 'migration target', ylab = 'population migration rate (2Nm)',
                 saveplot = FALSE)


##### PANEL C: MIG FROM GUI VS MAM ######
Logsel <- filter(Log, var == '2Nm' & migRun == 'A2anc.U2anc')
p.2Nm.AU <- vplot(Logsel,
                  xvar = 'migto', fillvar = 'cn', colvar = 'migfrom', yvar = 'val',
                  linecols = 'pop.cols', y.max = 'max.hpd', rm.violins = TRUE, hpdline.width = 1,
                  legpos = 'top', legcolname = "migration source:",
                  xlab = 'migration target', ylab = 'population migration rate (2Nm)',
                  saveplot = FALSE)


##### COMBINE PANELS #####
p.2Nm.U <- p.2Nm.U + theme(plot.margin = margin(0.5, 0.75, 0.5, 0.5, "cm"))
p.2Nm.AU <- p.2Nm.AU + theme(plot.margin = margin(0.5, 0.75, 0.5, 0.5, "cm"))
p.2Nm.W <- p.2Nm.W + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

plots <- ggarrange(p.2Nm.U, p.2Nm.AU, p.2Nm.W, ncol = 3, nrow = 1)
plots <- plots + draw_plot_label(label = c("A", "B", "C"), size = 24, x = c(0, 0.33, 0.67), y = c(1, 1, 1))
ggexport(plots, filename = gphocs.suppfig.file, width = 1100, height = 450)
