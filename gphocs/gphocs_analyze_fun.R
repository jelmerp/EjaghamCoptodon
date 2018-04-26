#### SET-UP ####

## Files to read:
popcols.file <- 'metadata/popcolours.txt'
gphocs_output.file <- 'analyses_output/gphocs_mergedOutput.txt'

## Libraries:
library(data.table)
library(TeachingDemos)
library(plyr)
library(reshape2)
library(tidyverse)
library(RColorBrewer)

## Further set-up:
Log <- as.data.frame(fread(gphocs_output.file, stringsAsFactors = TRUE))

inds <- read.table(popcols.file, header = TRUE, stringsAsFactors = FALSE)
inds.arr <- arrange(inds[!is.na(inds$pop.rank), ], pop.rank)

kidpops <- c('Dec', 'Eja', 'Fus', 'DE', 'DEF', 'Gui', 'Mam', 'UA', 'AU')
parentpops <- c('DE', 'DE', 'DEF', 'DEF', 'root', 'AU', 'AU', 'root', 'root')
currentpops <- c('Dec', 'Eja', 'Fus', 'Mam', 'Gui')
ancpops <- c('DE', 'DEF', 'AU', 'root')
pops <- c('Dec', 'Eja', 'DE', 'Fus', 'DEF', 'Mam', 'Gui',  'UA', 'root')


#### DATA PROCESSING FUNCTIONS ####

## Master: merge log files into a single df and edit:
mergelogs <- function(logfolder, burnIn, lastSample, subSample) {

  logfiles <- list.files(paste0(logfolder, '/raw'), pattern = '\\.log')

  Log <- lapply(logfiles, cutPrepLog, logfolder = logfolder,
                burnIn = burnIn, lastSample = lastSample, subSample = subSample)

  if(length(Log) > 1) Log <- do.call(rbind, Log) else Log <- Log[[1]]
  cat('Done processing individuals logs.\n')

  Log$migtoRun <- factor(Log$migtoRun, levels = inds.arr$pop.short)
  Log$migfromRun <- factor(Log$migfromRun, levels = inds.arr$pop.short)
  Log$migto <- factor(Log$migto, levels = inds.arr$pop.short)
  Log$pop <- factor(Log$pop, levels = inds.arr$pop.short)

  Log$phyl <- factor('model')
  Log$cn <- factor('aa')

  Log$migfromType <- factor(inds$migfrom.type[match(Log$migRun, inds$pop.short)],
                            levels = c('non', 'single', 'mult'))
  Log$migtoType <- factor(inds$migto.type[match(Log$migRun, inds$pop.short)])
  Log$migtoType2 <- factor(inds$migto.type2[match(Log$migRun, inds$pop.short)])
  Log$migtoType3 <- factor(gsub('anc|cur', 'anc\\|cur', Log$migtoType2), levels = c('none', 'single', 'anc|cur', 'all'))

  Log <- Log[, c('Sample', 'var', 'val', 'cval', 'pop', 'migfrom', 'migto', 'rep', 'migRun',
                 'migfromRun', 'migtoRun', 'phyl', 'cn',
                 'migfromType', 'migtoType', 'migtoType2', 'migtoType3')]
  return(Log)
}


## Step 3: Master to first apply cutLog, then prepOneLog:
cutPrepLog <- function(logfile, logfolder, burnIn, lastSample, subSample,
                       return.log = TRUE, write.cutLog = TRUE) {

  Log <- cutLog(logfile = logfile, logfolder = logfolder, burnIn = burnIn,
                lastSample = lastSample, subSample = subSample,
                return.log = return.log, write.cutLog = write.cutLog)

  Log <- prepOneLog(oneLog = Log, logfile = logfile)

  return(Log)
}

## Step 2: Melt dataframe and prep vars:
prepOneLog <- function(oneLog, logfile) {

  colnames(oneLog) <- gsub('\\.\\.', 2, colnames(oneLog))
  colnames(oneLog) <- gsub('Cgal', 'Cmam', colnames(oneLog))
  colnames(oneLog)[(ncol(oneLog)-2):ncol(oneLog)] <- c('mut_NA', 'dataLd_NA', 'FullLd_NA')

  migcols <- grep('m_', colnames(oneLog))
  for(migcol in migcols) {
    m <- oneLog[, migcol]
    migpattern <- unlist(strsplit(colnames(oneLog)[migcol], split = '_'))[2]
    migto <- unlist(strsplit(migpattern, split = '2'))[2]
    migto.column <- grep(paste0('theta_', migto, '$'), colnames(oneLog))
    th <- oneLog[, migto.column]

    ## Population migration rate:
    colname.mig2 <- paste0('2Nm_', migpattern)
    oneLog$newcolumn1 <- (m * m.scale) * (th * t.scale) / 4
    colnames(oneLog)[grep('newcolumn1', colnames(oneLog))] <- colname.mig2

    ## Proportion of migrants:
    colname.mig3 <- paste0('m.prop_', migpattern)
    oneLog$newcolumn2 <- (m * m.scale) * mutrate.gen
    colnames(oneLog)[grep('newcolumn2', colnames(oneLog))] <- colname.mig3

    ## Total migration rate:
    # colname.mig3 <- paste0('m.prop_', migpattern)
    # oneLog$newcolumn2 <- (m * m.scale) * mutrate.gen
    # colnames(oneLog)[grep('newcolumn2', colnames(oneLog))] <- colname.mig3
  }

  mlog <- melt(oneLog, id = 'Sample') %>% separate(variable, sep = '_', into = c('var', 'population'))
  mlog$population <- factor(mapvalues(mlog$population, from = inds$pop.gphocsfile1,
                                      to = inds$pop.short, warn_missing = FALSE))
  mlog$pop <- mlog$population

  mlog <- mlog %>% separate(pop, sep = '2', into = c('migfrom', 'migto'))
  mlog$migfrom[which(is.na(mlog$migto))] <- NA
  mlog$pop[which(!is.na(mlog$migto))] <- NA
  mlog <- dplyr::rename(mlog, val = value)
  mlog <- mlog[, c('Sample', 'var', 'val', 'pop', 'migfrom', 'migto')]

  mlog$migfrom <- factor(mlog$migfrom)
  mlog$migto <- factor(mlog$migto)

  migRun <- mapvalues(unlist(strsplit(logfile, split = '_'))[1],
                      from = inds$pop.gphocsfile2, to = inds$pop.short, warn_missing = FALSE)
  mlog$migRun <- factor(migRun)
  mlog$rep <- factor(as.integer(gsub('.*rep([0-9])\\.log', '\\1', logfile)))
  migRunSplit <- unlist(strsplit(migRun, split = '2'))
  mlog$migfromRun <- factor(migRunSplit[1])
  mlog$migtoRun <- factor(ifelse(length(migRunSplit) == 2, migRunSplit[2], migRunSplit[1]))

  mlog <- add.cvalue(mlog)

  return(mlog)
}

## Step 1: Open one logfile and cut off burn-in and/or last samples:
cutLog <- function(logfile, logfolder, burnIn, lastSample, subSample,
                   return.log = TRUE, write.cutLog = TRUE) {

  Log <- read.table(paste0(logfolder, '/raw/', logfile), header = TRUE)

  ## Remove burn-in:
  if(burnIn > 0) Log <- Log[-which(Log$Sample < burnIn), ]

  ## Remove final samples:
  if(!is.null(lastSample)) if(any(Log$Sample > lastSample)) Log <- Log[-which(Log$Sample > lastSample), ]

  ## Subsample:
  Log <- Log[seq(from = 1, to = nrow(Log), by = subSample), ]

  cat(logfile, "\tnrows:", nrow(Log), "\tLast sample:", max(Log$Sample),
      "\tBurnIn:", burnIn, "\tSubsample:", subSample, '\n')

  if(write.cutLog == TRUE) {
    if(!dir.exists(paste0(logfolder, '/cut/'))) dir.create(paste0(logfolder, '/cut/'))
    write.table(Log, paste0(logfolder, '/cut/', logfile), sep = '\t', quote = FALSE, row.names = FALSE)
  }
  if(return.log == TRUE) return(Log)
}


#### HELPER FUNCTIONS ####

## Add converted demographic values:
add.cvalue <- function(Log) {
  # Log <- mlog

  Log$cval <- NA

  ## Migration rate:
  smr <<- Log %>% filter(var == 'tau') %>% dplyr::group_by(pop) %>% dplyr::summarise(tau = mean(val))
  focalpops <- as.character(unique(Log$migto))
  focalpops <- focalpops[!is.na(focalpops)]
  if(length(focalpops) > 1) newmig <- do.call(rbind, lapply(focalpops, getlifespan, Log = Log, smr = smr))
  if(length(focalpops) == 1) newmig <- getlifespan(focalpops, Log = Log, smr = smr)

  if(length(focalpops) >= 1) Log$cval[newmig$frows] <- (Log$val[newmig$frows] * m.scale) * (newmig$lifespan * t.scale)

  ## Tau & theta:
  Log$cval[Log$var == 'theta'] <- (Log$val[Log$var == 'theta'] * t.scale) / (4 * mutrate.gen) # Ne = theta / 4 Ne u
  Log$cval[Log$var == 'tau'] <- (Log$val[Log$var == 'tau'] * t.scale) / (mutrate.gen) # time = tau / mutrate

  return(Log)
}

## Get lifespan for a pop in a specific run:
getlifespan <- function(focalpop, Log, smr) {
  #focalpop <- focalpops

  if(is.null(focalpop)) focalpop <- unlist(strsplit(migRun, split = '2'))[2]

  tau.parent <- smr$tau[smr$pop == getparent(focalpop)]

  if(focalpop %in% currentpops) lifespan <- tau.parent else lifespan <- tau.parent - smr$tau[smr$pop == focalpop]

  frows <- which(Log$var == 'm' & Log$migto == focalpop)

  return(data.frame(frows, lifespan))
}

## Get parent pop:
getparent <- function(kidpop) parentpops[match(kidpop, kidpops)]

## Get intermediate colour:
get.midpoint <- function(col1, col2) {
  col <- rgb(red = (col2rgb(col1)[1] + col2rgb(col2)[1]) / 2,
             green = (col2rgb(col1)[2] + col2rgb(col2)[2]) / 2,
             blue = (col2rgb(col1)[3] + col2rgb(col2)[3]) / 2,
             maxColorValue = 255)
  return(col)
}

## Functions for Highest Posterior Density (HPD):
library(TeachingDemos) # has emp.hpd function
hpd.min <- function(x) emp.hpd(x)[1]
hpd.max <- function(x) emp.hpd(x)[2]


## Get full ("gf") population name:
gf <- function(short) inds$pop.full[match(short, inds$pop.short)]

## Get converted variable name (e.g. theta --> Ne):
cvar <- function(variable) {
  if(variable == 'theta') return(expression(N[e] ~ "(in 1000s)"))
  if(variable == 'tau') return("divergence time (ka ago)")
}

## Summarize theta and tau parameters:
tt.prep <- function(Logg, x.even = FALSE, river.even = FALSE, river.popsize = 1,
                    summary.provided = FALSE) {
  # subset(Log, migRun == 'Mam2all'); x.even = FALSE; summary.provided = FALSE

  if(summary.provided == FALSE) {
    tt <- subset(Logg, var %in% c('theta', 'tau')) %>% group_by(pop, var) %>%
      dplyr::summarise(cval.mean = mean(cval / 1000)) %>%
      dcast(pop ~ var)
  }

  tt$NeToScale <- 1
  tt$x.min <- NA

  if(x.even != TRUE) {
    tt$NeToScale[match(c('Mam', 'Gui', 'AU', 'root'), tt$pop)] <- 0

    if(river.even != TRUE) {
    tt$theta[match('Mam', tt$pop)] <- river.popsize
    tt$theta[match('Gui', tt$pop)] <- river.popsize
    tt$theta[match('AU', tt$pop)] <- river.popsize
    tt$theta[match('root', tt$pop)] <- river.popsize
    }
  }

  if(x.even == TRUE) {
    tt$theta <- c(1, 1, 1.5, 1, 1, 1, 1, 1, 1)
    tt$tau <- c(9.5, 8, 6, 3, NA, NA, NA, NA, NA)
  }

  tt$x.min[match('Gui', tt$pop)] <- 0.5
  tt$x.min[match('AU', tt$pop)] <- round(tt$x.min[match('Gui', tt$pop)] + tt$theta[match('Gui', tt$pop)], 2)
  tt$x.min[match('Mam', tt$pop)] <- round(tt$x.min[match('AU', tt$pop)] + tt$theta[match('AU', tt$pop)], 2)
  tt$x.min[match('root', tt$pop)] <- round(tt$x.min[match('AU', tt$pop)] + tt$theta[match('AU', tt$pop)], 2)

  tt$x.min[match('Fus', tt$pop)] <- round(tt$x.min[match('Mam', tt$pop)] + tt$theta[match('Mam', tt$pop)], 2) + 0.25
  tt$x.min[match('DEF', tt$pop)] <- round(tt$x.min[match('Fus', tt$pop)] + tt$theta[match('Fus', tt$pop)], 2)
  tt$x.min[match('DE', tt$pop)] <- round(tt$x.min[match('DEF', tt$pop)] + tt$theta[match('DEF', tt$pop)], 2)
  tt$x.min[match('Dec', tt$pop)] <- round(tt$x.min[match('DE', tt$pop)] + tt$theta[match('DE', tt$pop)], 2)
  tt$x.min[match('Eja', tt$pop)] <- round(tt$x.min[match('DE', tt$pop)] - tt$theta[match('Eja', tt$pop)], 2)

  tt$x.max <- round(tt$x.min + tt$theta, 2)
  tt$x.max[match('root', tt$pop)] <- tt$x.min[match('DEF', tt$pop)]

  tt$y.min <- round(ifelse(tt$pop %in% currentpops, 0, tt$tau), 2)
  tt$y.max <- round(ifelse(tt$pop == 'root', tt$y.min + 0.7, tt$tau[match(getparent(tt$pop), tt$pop)]), 2)

  tt$popcol <- inds$pop.col[match(tt$pop, inds$pop.short)]

  return(tt)
}

## Summarize migration parameters:
m.prep <- function(Log) {
  m <- dplyr::filter(Log, var == 'm') %>% group_by(migfrom, migto, var) %>% dplyr::summarise(value = mean(cval, na.rm = TRUE))
  return(m)
}


#### PLOTTING FUNCTIONS ####

## Demography plot:
dplot <- function(tt, m = NULL, x.min = NULL, y.max = NULL, x.even = FALSE,
                  col.scale = TRUE, drawconnection = TRUE,
                  ann.pops = TRUE, cur.below.xaxis = FALSE, popnames.size = 7,
                  popnames.adj.vert = 0.05, popnames.adj.horz = 0,
                  xlab = expression(N[e] ~ "(1 tick mark = 1000)"), ylab = 'time (ka ago)', plot.title = NULL,
                  plotfolder = NULL, saveplot = TRUE, filetype = 'pdf',
                  plotfilewidth = 6, plotfileheight = 6, file.open = TRUE, filename = NULL) {

  if(saveplot == FALSE) file.open <- FALSE

  p <- ggplot()
  if(x.even == TRUE | col.scale == FALSE) {
    p <- p + geom_rect(data = tt, colour = 'grey20',
                       aes(xmin = x.min, xmax = x.max, ymin = y.min, ymax = y.max, fill = popcol))
  } else {
    p <- p + geom_rect(data = tt,
                       aes(xmin = x.min, xmax = x.max, ymin = y.min, ymax = y.max, fill = popcol, color = factor(NeToScale)))
    p <- p + scale_colour_manual(breaks = c(0, 1), values = c('grey40', 'grey10'))
    p <- p + guides(colour = FALSE)
  }

  p <- p + scale_fill_identity()

  p <- p + theme_bw()
  p <- p + theme(axis.text.y = element_text(size = 18))
  if(is.null(xlab)) p <- p + theme(axis.title.x = element_blank())
  if(!is.null(xlab)) p <- p + theme(axis.title.x = element_text(size = 18))
  p <- p + theme(axis.title.y = element_text(size = 20))
  p <- p + theme(plot.title = element_text(size = 20, hjust = 0.5))
  p <- p + theme(plot.margin = unit(c(1, 1.2, 1, 0.5), "lines")) # top, right, ..

  p <- p + scale_y_continuous(expand = c(0, 0), breaks = seq(from = 0, to = 25, by = 1))
  p <- p + scale_x_continuous(expand = c(0.01, 0.01), breaks = seq(from = 0, to = 25, by = 1))

  if(ann.pops == TRUE) {
    if(cur.below.xaxis == TRUE) {
      x.breaks <- (tt[tt$pop %in% currentpops, ]$x.min + tt[tt$pop %in% currentpops, ]$x.max) / 2
      x.labs <- tt[tt$pop %in% currentpops, 'pop']
      p <- p + scale_x_continuous(expand = c(0, 0), breaks = x.breaks, labels = x.labs)
      p <- p + theme(axis.text.x = element_blank())

      x.locs <- (tt[tt$pop %in% ancpops, ]$x.min + tt[tt$pop %in% ancpops, ]$x.max) / 2
      y.locs <- tt[tt$pop %in% ancpops, ]$y.min + popnames.adj.vert
      anclabs <- tt[tt$pop %in% ancpops, 'pop']
      p <- p + annotate(geom = "text", x = x.locs, y = y.locs,
                        label = anclabs, color = "black", size = popnames.size)
    }
    if(cur.below.xaxis == FALSE) {
      x.locs <- ((tt$x.min + tt$x.max) / 2) + popnames.adj.horz
      y.locs <- tt$y.min + popnames.adj.vert
      p <- p + annotate(geom = "text", x = x.locs, y = y.locs,
                        label = tt$pop, color = "black", size = popnames.size)
      p <- p + theme(axis.text.x = element_blank())
    }
  }
  if(ann.pops == FALSE) p <- p + theme(axis.text.x = element_blank())

  if(x.even == FALSE & drawconnection == TRUE) {
    line.y <- tt$tau[match('root', tt$pop)]
    line.left <- tt$x.max[match('AU', tt$pop)]
    line.right <- tt$x.min[match('DEF', tt$pop)]
    p <- p + geom_segment(aes(x = line.left, y = line.y, xend = line.right, yend = line.y), colour = "grey20")
  }

  if(x.even == TRUE & drawconnection == TRUE) {
    line.y <- tt$tau[match('root', tt$pop)]
    line.left <- tt$x.max[match('DEF', tt$pop)]
    line.right <- tt$x.min[match('AU', tt$pop)]
    p <- p + geom_segment(aes(x = line.left, y = line.y, xend = line.right, yend = line.y), colour = "grey20")
  }

  if(!is.null(plot.title)) p <- p + labs(title = plot.title)
  if(!is.null(xlab)) p <- p + labs(x = xlab)
  if(!is.null(ylab)) p <- p + labs(y = ylab)

  if(!is.null(x.min)) p <- p + coord_cartesian(xlim  = c(x.min, Eja_leftlim + Eja_Ne + 0.5))
  if(!is.null(y.max)) p <- p + coord_cartesian(ylim = c(0, y.max))

  if(is.null(plotfolder)) plotfolder <- paste0('analyses/gphocs/plots')
  plotfile <- paste0(plotfolder, '/', filename, '.', filetype)
  if(saveplot == TRUE) ggsave(filename = plotfile, plot = p, width = plotfilewidth, height = plotfileheight)
  if(file.open == TRUE) system(paste("xdg-open", plotfile))

  plot(p)
  return(p)
}


## Violin plot:
vplot <- function(data, xvar, yvar = 'val', fillvar = 'cn', colvar = 'cn', shapevar = 'cn',
                  linecols = 'grey30', fillcols = NULL,
                  shade = TRUE, shadecol = 'grey80',
                  y.min = 0, y.max = 'max.hpd', ylims.dft = FALSE, ymax.expand = 0.05,
                  facetwrap = NULL, rm.violins = FALSE,
                  statsum = TRUE, meandot.size = 2, hpdline.width = 0.8,
                  legpos = 'none', legfillname = NULL, legcolname = NULL, legend.nrow = 1,
                  rm.leg.col = TRUE, rm.leg.fill = FALSE, rm.leg.shape = FALSE,
                  plot.title = NULL, xlab = '', ylab = '',
                  saveplot = TRUE, filename = NULL, filetype = 'png',
                  plotfilewidth = 8, plotfileheight = 8, plotfolder = NULL, file.open = TRUE) {

  if(saveplot == FALSE) file.open <- FALSE

  if(nrow(data) == 0) stop("Error: no rows in data")

  if(unique(data$var) %in% c('theta', 'tau') & yvar == 'cval') data$cval <- data$cval / 1000
  if(unique(data$var) == 'm.prop') data$val <- data$val * 100

  ## Create base plot:
  p <- ggplot()

  if(rm.violins == FALSE) {
    p <- p + geom_violin(data = data, aes_string(x = xvar, y = yvar, fill = fillvar, colour = colvar))
  } else {
    p <- p + geom_blank(data = data, aes_string(x = xvar, y = yvar, fill = fillvar, colour = colvar))
  }

  if(shade == TRUE) {
    nrvars <- length(unique(data[, xvar]))
    rect_left = c(seq(from = 0.5, to = nrvars, by = 2))
    rectangles <- data.frame(x.min = rect_left, x.max = rect_left + 1)
    if(length(shadecol) == 1) shadecol <- rep(shadecol, nrow(rectangles))
    p <- p + geom_rect(data = rectangles, aes(xmin = x.min, xmax = x.max, ymin = -Inf, ymax = Inf),
                       fill = shadecol, colour = shadecol, alpha = 0.5)
  }

  if(rm.violins == FALSE) {
    p <- p + geom_violin(data = data, aes_string(x = xvar, y = yvar, fill = fillvar, colour = colvar))
  } else {
    p <- p + geom_blank(data = data, aes_string(x = xvar, y = yvar, fill = fillvar, colour = colvar))
  }

  ## Colours and fills:
  if(length(fillcols) == 1) if(fillcols == 'pop.cols') {
      levels.sorted <- levels(data[, fillvar])[sort(match(unique(data[, fillvar]), levels(data[, fillvar])))]
      fillcols <- inds.arr$pop.col[match(levels.sorted, inds.arr$pop.short)]
  }
  if(length(fillcols) == 1) if(fillcols != 'pop.cols') fillcols <- rep(fillcols, length(unique(data[, fillvar])))

  if(length(linecols) == 1) if(linecols == 'pop.cols') {
    levels.sorted <- levels(data[, colvar])[sort(match(unique(data[, colvar]), levels(data[, colvar])))]
    linecols <- inds.arr$pop.col[match(levels.sorted, inds.arr$pop.short)]
  }
  if(length(linecols) == 1) if(linecols != 'pop.cols') linecols <- rep(linecols, length(unique(data[, colvar])))

  if(is.null(linecols)) linecols <- brewer.pal(n = length(unique(data[, colvar])), name = 'Set1')
  if(is.null(fillcols)) fillcols <- brewer.pal(n = length(unique(data[, colvar])), name = 'Set1')

  p <- p + scale_colour_manual(values = linecols, name = legcolname)
  p <- p + scale_fill_manual(values = fillcols, name = legfillname)

  ## Facet wrap:
  if(!is.null(facetwrap)) p <- p + facet_wrap(~ rep, nrow = 2)

  ## Compute and show data summaries:
  if(statsum == TRUE) {
    p <- p + stat_summary(data = data, aes_string(x = xvar, y = yvar, fill = fillvar, color = colvar, width = 0.4),
                          fun.ymin = hpd.min, fun.ymax = hpd.max, geom = "errorbar",
                          size = hpdline.width, position = position_dodge(width = 0.9))
    p <- p + stat_summary(data = data, aes_string(x = xvar, y = yvar, fill = fillvar, color = colvar, shape = shapevar),
                          fun.y = mean, geom = "point",
                          size = meandot.size, position = position_dodge(width = 0.9))
  }

  ## Axis limits:
  if(ylims.dft == FALSE) {
    if(y.max %in% c('max.hpd', 'max.value')) {
      if(y.max == 'max.hpd' & yvar == 'val') df <- group_by_(data, fillvar, colvar, xvar) %>% dplyr::summarise(Max = hpd.max(val))
      if(y.max == 'max.value' & yvar == 'val') df <- group_by_(data, fillvar, colvar, xvar) %>% dplyr::summarise(Max = max(val))
      if(y.max == 'max.hpd' & yvar == 'cval') df <- group_by_(data, fillvar, colvar, xvar) %>% dplyr::summarise(Max = hpd.max(cval))
      if(y.max == 'max.value' & yvar == 'cval') df <- group_by_(data, fillvar, colvar, xvar) %>% dplyr::summarise(Max = max(cval))
      max <- as.numeric(max(df$Max))
      y.max <- max + (ymax.expand * max)
    }
    p <- p + coord_cartesian(ylim = c(y.min, y.max))
  }
  p <- p + scale_y_continuous(expand = c(0, 0))

  ## Axis and plot titles/labels:
  if(!is.null(plot.title)) p <- p + labs(title = plot.title)
  p <- p + labs(x = xlab, y = ylab)

  ## Formatting:
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(size = 18),
                 axis.text.y = element_text(size = 18),
                 axis.title.x = element_text(size = 20),
                 axis.title.y = element_text(size = 20),
                 plot.title = element_text(size = 20, hjust = 0.5))

  ## Legend formatting:
  p <- p + theme(legend.title = element_text(size = 15, face = 'bold'),
                 legend.text = element_text(size = 15),
                 legend.key.height = unit(0.5, "cm"),
                 legend.key.width = unit(0.5, "cm"),
                 legend.background = element_rect(fill = "grey90", colour = "grey30"),
                 legend.key = element_rect(fill = "grey90"),
                 legend.position = legpos)

  ## Remove legend for constants:
  if(colvar == 'cn' | rm.leg.col == TRUE) p <- p + guides(colour = FALSE)
  if(fillvar == 'cn' | rm.leg.fill == TRUE) p <- p + guides(fill = FALSE)
  if(shapevar == 'cn' | rm.leg.shape == TRUE) p <- p + guides(shape = FALSE)

  ## Edit legend if not plotting violins:
  if(rm.violins == TRUE) p <- p + guides(colour = guide_legend(override.aes = list(linetype = 1, shape = 16)))

  ## Legend across multiple rows:
  p <- p + guides(colour = guide_legend(nrow = legend.nrow, byrow = TRUE))

  ## Save plot:
  if(is.null(plotfolder)) plotfolder <- paste0('analyses/gphocs/plots')
  if(!dir.exists(plotfolder)) dir.create(plotfolder)
  plotfile <- paste0(plotfolder, '/', filename, '.', filetype)

  if(saveplot == TRUE) {
    cat('Saving plot:', plotfile, '\n')
    ggsave(filename = plotfile, plot = p, width = plotfilewidth, height = plotfileheight)
  }
  if(file.open == TRUE) system(paste("xdg-open", plotfile))

  print(p)
  return(p)
}

