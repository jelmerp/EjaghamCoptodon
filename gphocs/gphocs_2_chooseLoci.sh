#!/bin/bash -l

SCAFFOLD=$1

date
echo "Script: gphocs_avoidHardMaskByScaffold.sh"
echo "Scaffold: $SCAFFOLD"

Rscript scripts/gphocs/gphocs_2_chooseLoci.R $SCAFFOLD $PWD

echo "Done with script."
date