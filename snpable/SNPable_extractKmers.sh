#!/bin/bash

date
echo "Script: SNPable_extractKmers.sh..."

REFERENCE=$1
SNPABLE=$2

cd kmers
$SNPABLE/splitfa $REFERENCE 50 | split -l 20000000

echo "Done with script."
date