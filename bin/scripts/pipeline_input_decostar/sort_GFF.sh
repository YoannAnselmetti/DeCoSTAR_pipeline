#!/bin/bash
###
###   Goal:
###      Script to sort GFF files that will be used to create GENE files
###
###   INPUT:
###      1- INPUT directory containing GFF files unsorted
###         (data/INPUT_DATA/GFF)
###      2- OUTPUT directory where will be stored GFF files sorted
###         (data/GFF_to_GENE_files/sorted_GFF)
###   OUTPUT:
###      - OUTPUT directory containing GFF files sorting by gene order in genome species (GGF sorted file name = species name)
###
###   Name: sort_GFF.sh             Author: Yoann Anselmetti
###   Creation date: 2015/07/21     Last modification: 2017/10/18
###


if [ -d $2 ]; then
   rm -rf $2/*
else
   mkdir -p $2
fi

# Course all GFF files unsorted in INPUT directory
for GFFfile in $(ls $1)
do
   fullGFFfile="$1/"$GFFfile
   echo $fullGFFfile
   name_species=$(echo $GFFfile | cut -d. -f1)
   if [[ ! -z $name_species ]];
	then
		echo -e "\t=> "$name_species
		# Sort GFF files and store it: OUTPUT directory with file name = species name
		sort -k9d,9d -k4n,4n $fullGFFfile > $2/$name_species"_sorted.gff3"
	fi
done
