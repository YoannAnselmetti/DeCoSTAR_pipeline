#!/bin/bash
###
###   Goal:
###      Get number of files in a directory (browse directory tree)
###
###   INPUT:
###      1- gene_TAG-species_name association file
###         (data/INPUT_DATA/name_geneID_18anopheles)
###      2- INPUT file to modify
###         (data/GENE_TREES/unrooted_trees-ASTEI+filt.nwk)
###      3- OUTPUT file modified
###         (data/GENE_TREES/unrooted_trees-ASTEI+filt+spe.nwk)
###      4- Character separator between species name and gene TAG ID
###         (@)
###
###   OUTPUT:
###      - Number of files contained in INPUT directory
###
###   Name: 16bisA-add_speciesName_to_geneID_for_RAW_trees.sh	Author: Yoann Anselmetti
###   Creation date: 2016/09/06									Last modification: 2016/11/15
###

name_geneID="data/INPUT_DATA/name_geneID_18anopheles"
inputF="data/GENE_TREES/unrooted_trees-ASTEI+filt.nwk"
outputF="data/GENE_TREES/unrooted_trees-ASTEI+filt+spe.nwk"
sep="@"

cp $inputF $outputF

while read line; do
	name=$(echo $line|cut -d" " -f1)
	ID=$(echo $line|cut -d" " -f2)
	echo "sed -i \"s/$ID/$name$sep$ID/g\" $outputF"
	sed -i "s/$ID/$name$sep$ID/g" $outputF
done < ${name_geneID}
