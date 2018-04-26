#### PRODUCE PHYLOGENETIC TREES USING SEVERAL METHODS ####

## RaxML - entire genome
Use raxml_master.sh which will call raxml_run.sh

## RaxML - gene trees by window
First, generate window coordinates using raxml_gnrWindowCoords.R. Then, use raxml.window.master.sh which will call raxml.window.sub.sh, raxml.window.sub2.sh, raxml_run.sh.

## Phylonet
See phylonet_master.sh to first prepare input files using phylonet_prepnexus.R and then run Phylonet using phylonet_run.sh.

## ASTRAL
See astral.sh

## MP-EST
See mpest.sh

## Splitstree
See splitstree_master.sh, which will call splitstree_run.sh
