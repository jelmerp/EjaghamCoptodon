#!/bin/bash -l

CONTROL_FILE=$1
NPROC=$2
GPHOCS=software/gphocs_1.3/G-PhoCS/bin/G-PhoCS
  
date
echo "Script: gphocs_5_run.sh"
echo "Control file: $CONTROL_FILE"
echo "Number of processors: $NPROC"

## Run:
export OMP_NUM_THREADS=$NPROC
$GPHOCS $CONTROL_FILE -n $NPROC
