#!/bin/bash
###
###   Goal: 
###      Script to sort GENE files
###
###   INPUT:
###      1- INPUT directory containing GENE files unsorted
###         (data/GFF_to_GENE_files/GENE)
###      2- OUTPUT directory where will be stored sorted GENE files
###         (data/GFF_to_GENE_files/sorted_GENE)
###   OUTPUT:             
###      - OUTPUT directory containing GENE files sorted by gene order on the genome of species
###         
###   Name: sort_GENE.sh            Author: Yoann Anselmetti
###   Creation date: 2015/07/22     Last modification: 2017/10/18
###         


if [ -d $2 ]; then
   rm -rf $2/*
else
   mkdir -p $2
fi

# Course all GFF files unsorted in INPUT directory
for GENEfile in $(ls $1)
do
   fullGENEfile="$1/"$GENEfile
   echo $fullGENEfile
   name_species=$(echo $GENEfile | cut -d. -f1)
   # If $name_species is not empty 
   if [[ ! -z $name_species ]]; then
      echo -e "\t=> "$name_species
      # Sort GFF files and store it OUTPUT directiory with file name = species name
      sort -k2d,2d -k5n,5n -k6rn,6rn $fullGENEfile > $2/$name_species"_sorted.txt"
   fi
done
