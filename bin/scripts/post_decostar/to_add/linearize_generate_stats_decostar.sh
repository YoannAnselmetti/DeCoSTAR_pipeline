# Cedric Chauve, December 20, 2016
# Linearizing and analyzing evolution for the genomes of the decostar datasets

#! /bin/bash

# $1: DeCoSTAR OUTPUT directory
# $2: DeCoSTAR OUTPUT files prefix
# $3: Gene file with gene family ID
# $4: Linearization threshold value

longPREF=$1/$2
PREF=$2_Lin$4
scriptLINEAR=bin/scripts/post_decostar/code/linearize_assign_scj.sh

$scriptLINEAR ${longPREF}.reconciliations.newick ${longPREF}.speciesTree.newick ${longPREF}.adjacencies.txt $3 $4 $1 $PREF M1 python2