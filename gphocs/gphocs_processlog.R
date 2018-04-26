task <- commandArgs()[1] # task <- 'process' # task <- 'cut' # task <- 'load'
logfolder <- commandArgs()[2] # logfolder <- 'output/EjaC.Cmam_0307' #logfolder <- 'output/EjaC.Cmam.Cgui_0307' #logfolder <- 'output/EjaC.Cgui_0308'
if(length(commandArgs()) > 2) burnIn <- as.integer(commandArgs()[3]) else burnIn <- 250000
if(length(commandArgs()) > 3) subSample <- as.integer(commandArgs()[4]) else subSample <- 50
if(length(commandArgs()) > 4) lastSample <- as.integer(commandArgs()[5]) else lastSample <- NULL

if(length(commandArgs()) > 5) m.scale <- commandArgs()[6] else m.scale <- 1000
if(length(commandArgs()) > 6) t.scale <- commandArgs()[7] else t.scale <- 0.0001
if(length(commandArgs()) > 7) mutrate.gen <- commandArgs()[8] else mutrate.gen <- 7.5e-9
if(length(commandArgs()) > 8) gentime <- commandArgs()[9] else gentime <- 1

# logfolder <- 'output/seqfile170617/focal'; burnIn <- 200000; subSample <- 5; lastSample <- NULL
# m.scale <- 1000; t.scale <- 0.0001; mutrate.gen <- 6.6e-8; gentime <- 1

cat('task:', task, '\n')
cat('logfolder:', logfolder, '\n')
cat('burnIn:', burnIn, '\n')
cat('lastSample:', lastSample, '\n')
cat('subSample:', subSample, '\n')
cat('m.scale:', m.scale, '\n')
cat('t.scale:', t.scale, '\n')
cat('mutrate:', mutrate.gen, '\n')
cat('gentime:', gentime, '\n')

## Populations etc:
mutrate.yr <- mutrate.gen / gentime

model <- unlist(strsplit(unlist(strsplit(logfolder, split = '/'))[2], split = '_'))[1]
if(model == 'EjaC.Cmam') migfrom.pop <- 'Cmam'
if(model == 'EjaC.Cgui') migfrom.pop <- 'Cgui'
if(model == 'EjaC.Cmam.Cgui') migfrom.pop <- 'Cmam/Cgui'

kidpops <- c('Dec', 'Eja', 'Fus', 'DE', 'DEF', 'Gui', 'Mam', 'UA')
parentpops <- c('DE', 'DE', 'DEF', 'DEF', 'root', 'AU', 'AU', 'root')
if(model != 'EjaC.Cmam.Cgui' & model != 'EjaC.AU') parentpops[grep('DEFG', parentpops)] <- 'root'
currentpops <- c('Dec', 'Eja', 'Fus', 'Mam', 'Gui')
ancpops <- c('DE', 'DEF', 'AU', 'root')
pops <- c('Dec', 'Eja', 'DE', 'Fus', 'DEF', 'Mam', 'Gui',  'UA', 'root')


#### PROCESS OR LOAD LOGDATA ####
## Process:
if(task == 'process') {
  Log <- mergelogs(logfolder = logfolder, burnIn = burnIn, lastSample = lastSample, subSample = subSample)
  saveRDS(Log, file = paste0(logfolder, '/mergedlog.RDS'))
}

## Load previously processed:
if(task == 'load') Log <- readRDS(paste0(logfolder, '/mergedlog.RDS'))

## Cut only:
if(task == 'cut') {
  for(logfolder in logfolders) {
    aap <- sapply(list.files(logfolder, pattern = 'log'), cutLog,
                  logfolder = logfolder, burnIn = burnIn, lastSample = lastSample, return.log = FALSE)
  }
}


