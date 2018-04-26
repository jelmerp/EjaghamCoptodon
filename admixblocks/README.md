#### ADMIXTURE BLOCK ANALYSES ####

## Step 1
Follow the steps in fd_master.sh to get per-window Fd-stats. This will call fd_run.sh and fd_submaster.sh.

## Step 2
Process the raw Fd output using the R scripts fd_process1_sumstats.R through fd_process5_filter.R to obtain admixture blocks.
Using hybridcheck_master.sh, fd_process3_runHC.R can be run non-interactively.

admixblocks_plotfun.R contains functions for Manhattan-style plots, and is used by the scripts to produce Fig 5 and Fig S5.


