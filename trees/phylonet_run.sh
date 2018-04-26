#!/bin/bash -l

TREEFILE="$1"
NEXUS_IN="$2"
NEXUS_OUT="$3"
PHYLONET_CMD="$4"
TAXONMAP="$5"
FOCALNETWORK_NAME="$6"
FOCALNETWORK="$7"
MAXMEM=$8

if [ -z $MAXMEM ]; then MAXMEM=4; fi

PHYLONET=software/PhyloNet_3.6.1.jar

date
echo "Script: phylonet_run.sh"
echo "Tree file: $TREEFILE"
echo "Nexus in: $NEXUS_IN"
echo "Nexus out: $NEXUS_OUT"
echo "Phylonet command: "$PHYLONET_CMD""
echo "Taxon map: "$TAXONMAP""
echo "Focal network (for CalGTProb): "$FOCALNETWORK""
echo "Focal network name (for CalGTProb): "$FOCALNETWORK_NAME""
echo "Max memory usage: "$MAXMEM""

echo "Preparing nexus file..."
Rscript scripts/phylonet/phylonet_prepnexus.R "$TREEFILE" "$NEXUS_IN" "$NEXUS_OUT" "$PHYLONET_CMD" "$TAXONMAP" "$FOCALNETWORK_NAME" "$FOCALNETWORK"

echo "Running phylonet..."
java -Xmx${MAXMEM}g -XX:+UseConcMarkSweepGC -jar $PHYLONET $NEXUS_IN
echo "Done."

grep "Dendroscope" $NEXUS_OUT*txt | cut -d ":" -f 2 | sed 's/^[ \t]*//;s/[ \t]*$//' > $NEXUS_OUT.phylonet.output.dendro
# sed command gets rid of leading and lagging whitespace

echo "Done with script."
date
