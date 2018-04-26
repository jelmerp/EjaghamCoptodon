#### DFOIL ON PREDEFINED ADMIXTURE BLOCKS ####

## Step 1
Use dfoil_blocks_master.sh to run DFOIL on the admixture blocks inferred previously (see admixblocks folder). This will call dfoil_inputprep.sh (to prepare input files) and dfoil_blocks_run.sh (to actually run dfoil).

## Step 2
Run dfoil_blocks_1_process.R and dfoil_blocks_2_filter.R to process and filter the output from Step 1.

#### GENOME-WIDE DFOIL ####
Use dfoil_entireGenome.sh to run DFOIL on the entire genome. This will call dfoil_inputprep.sh (to prepare input files) and dfoil_blocks_run.sh (to actually run dfoil).
