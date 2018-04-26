#!/bin/bash -l

## Command-line args:
NEXUS_IN=$1
NEXUS_OUT=$2

## Software:
SPLITSTREE=/nas02/home/j/e/jelmerp/splitstree4/SplitsTree

## Report:
date
echo "Script: splistree_run.sh"
echo "Nexus input: $NEXUS_IN"
echo "Nexus output: $NEXUS_OUT"

## Run Splitstree:
echo "Running Splitstree..."
$SPLITSTREE -g -i $NEXUS_IN -x "UPDATE; SAVE REPLACE=yes FILE=$NEXUS_OUT.tmp; QUIT"

## Remove sequence from Nexus output:
echo "Removing actual sequence from Nexus output..."
START=$(grep -n "BEGIN Characters;" $NEXUS_OUT.tmp | cut -f1 -d:)
END=$(grep -n "END;.*Characters" $NEXUS_OUT.tmp | cut -f1 -d:)
sed "$START,${END}d" $NEXUS_OUT.tmp > $NEXUS_OUT

rm $NEXUS_OUT.tmp

echo "Done with script."
date