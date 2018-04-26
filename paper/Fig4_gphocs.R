##### SET-UP #####

## Other scripts:
source('scripts/gphocs/gphocs_analyze_fun.R')

## Files to write:
fig4.file <- 'figures/Fig4_Gphocs.png'

## Libraries:
library(ggpubr)
library(cowplot)

## Settings: migration bands to show for tau and theta:
migsel <- c('none', 'single', 'ancestral/extant', 'all')


#### PANEL A: DEMOGRAPHY OVERVIEW ####
tt <- tt.prep(filter(Log, migRun == 'Mam2all'), river.popsize = 2)
adj.horz <- c(0, 0, 0, 0, 0, -1.5, 0, 0, 0)
p.dem <- dplot(tt, ann.pops = TRUE, popnames.size = 5, popnames.adj.vert = 0.3,
               popnames.adj.horz = adj.horz, drawconnection = FALSE, saveplot = FALSE)


#### PANEL B: TAU ####
Logsel <- filter(Log, var == 'tau' & migfromRun %in% c('none', 'Mam'))
p.tau <- vplot(Logsel,
               xvar = 'pop', fillvar = 'cn', colvar = 'migtoType3', shapevar = 'cn', yvar = 'cval',
               linecols = NULL, rm.violins = TRUE, hpdline.width = 1, y.max = 'max.hpd',
               legpos = 'top', rm.leg.fill = TRUE, legend.nrow = 2,
               xlab = paste('population'), ylab = cvar('tau'),
               plot.title = NULL, saveplot = FALSE)
p.tau <- p.tau + scale_colour_manual(values = brewer.pal(4, name = 'Set1')[1:4],
                                     name = "migration in run\n(from Mam to:)")
p.tau <- p.tau + theme(plot.margin = margin(0.75, 0.75, 0.5, 0.5, "cm"))


#### PANEL C: THETA ####
Logsel <- filter(Log, var == 'theta' & pop %in% c('Dec', 'Eja', 'Fus', 'DE', 'DEF') & migfromRun %in% c('none', 'Mam'))
p.theta <- vplot(Logsel,
                 xvar = 'pop', fillvar = 'cn', colvar = 'migtoType3', shapevar = 'cn', yvar = 'cval',
                 linecols = NULL, rm.violins = TRUE, hpdline.width = 1, y.max = 'max.hpd',
                 legpos = 'top', rm.leg.fill = TRUE, legend.nrow = 2,
                 xlab = paste('population'), ylab = cvar('theta'),
                 plot.title = NULL, saveplot = FALSE)
p.theta <- p.theta + scale_colour_manual(values = brewer.pal(4, name = 'Set1')[1:4],
                                         name = "migration in run\n(from Mam to:)")
p.theta <- p.theta + theme(plot.margin = margin(0.75, 0.75, 0.5, 0.75, "cm"))


#### PANEL D: 2NM ####
Logsel <- filter(Log, var == '2Nm' & migfromRun == 'Mam' & migfrom == 'Mam' & migto != 'Gui')
p.2Nm <- vplot(Logsel,
               xvar = 'migto', fillvar = 'cn', colvar = 'migtoType3', shapevar = 'cn', yvar = 'val',
               linecols = NULL, rm.violins = TRUE, hpdline.width = 1, y.max = 'max.hpd',
               legpos = 'top', rm.leg.fill = TRUE, legend.nrow = 2,
               xlab = paste('migration target'), ylab = 'population migration rate (2Nm)',
               saveplot = FALSE)
p.2Nm <- p.2Nm + scale_colour_manual(values = brewer.pal(4, name = 'Set1')[2:4],
                                     name = "migration in run\n(from Mam to:)")
p.2Nm <- p.2Nm + theme(plot.margin = margin(0.5, 0.75, 0.75, 0.75, "cm"))


#### PANEL E: TOTAL MIGRATION RATE ####
Logsel <- filter(Log, var == 'm' & migfromRun == 'Mam' & migfrom == 'Mam' & migto != 'Gui')
p.mtot <- vplot(Logsel,
                xvar = 'migto', fillvar = 'cn', colvar = 'migtoType3', shapevar = 'cn', yvar = 'cval',
                linecols = NULL, rm.violins = TRUE, hpdline.width = 1, y.max = 'max.hpd',
                legpos = 'top', rm.leg.fill = TRUE, legend.nrow = 2,
                xlab = paste('migration target'), ylab = 'M (total migration rate)',
                saveplot = FALSE)
p.mtot <- p.mtot + scale_colour_manual(values = brewer.pal(4, name = 'Set1')[2:4],
                                       name = "migration in run\n(from Mam to:)")
p.mtot <- p.mtot + theme(plot.margin = margin(0.5, 0.75, 0.75, 0.75, "cm"))


#### PANEL F: PROPORTION OF MIGRANTS ####
Logsel <- filter(Log, var == 'm.prop' & migfromRun == 'Mam' & migfrom == 'Mam' & migto != 'Gui')
p.mProp <- vplot(Logsel,
                 xvar = 'migto', fillvar = 'cn', colvar = 'migtoType3', shapevar = 'cn', yvar = 'val',
                 linecols = NULL, rm.violins = TRUE, hpdline.width = 1, y.max = 'max.hpd',
                 legpos = 'top', rm.leg.fill = TRUE, legend.nrow = 2,
                 xlab = paste('migration target'), ylab = expression(paste('migrant percentage')),
                 saveplot = FALSE)
p.mProp <- p.mProp + scale_colour_manual(values = brewer.pal(4, name = 'Set1')[2:4],
                                         name = "migration in run\n(from Mam to:)")
p.mProp <- p.mProp + theme(plot.margin = margin(0.5, 0.75, 0.75, 0.5, "cm"))


##### COMBINE PANELS #####
plots <- ggarrange(p.dem, p.tau, p.theta, p.2Nm, p.mtot, p.mProp, ncol = 3, nrow = 2)
plots <- plots + draw_plot_label(label = c('A', 'B', 'C', 'D', 'E', 'F'), size = 24,
                                 x = c(0, 0.33, 0.66, 0, 0.33, 0.66), y = c(1, 1, 1, 0.5, 0.5, 0.5))
ggexport(plots, filename = fig4.file, width = 1200, height = 900)
