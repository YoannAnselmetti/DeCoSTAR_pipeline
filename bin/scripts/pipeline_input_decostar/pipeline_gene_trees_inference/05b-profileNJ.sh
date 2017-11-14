#!/bin/bash
###                                                                       
###   Goal:                                                               
###      Contract weak edges of gene trees and resolve multifurcation     
###      by minimizing number of GDup and GLos                            
###                                                                       
###   INPUT:                                                              
###      1- PATH of profileNJ executive file                              
###         (Software/TREE_INFERENCE/profileNJ/bin/profileNJ)             
###      2- Species tree file                                             
###         (DECOSTAR/DATA_DeCoSTAR/INPUT_DATA/Anopheles_species_tree_X_topology.nwk)
###         (DECOSTAR/DATA_DeCoSTAR/INPUT_DATA/Anopheles_species_tree_WG_topology.nwk)
###      3- INPUT unrooted gene tree file for current Gene Family ID      
###         (DECOSTAR/DATA_DeCoSTAR/GENE_TREES/CDS/bootstrap_support/profileNJ/UNROOTED_GENE_TREES/gene_tree_GF0000001.nwk)
###      4- OUTPUT DIRECTORY to store resulting file of profileNJ         
###         (DECOSTAR/DATA_DeCoSTAR/GENE_TREES/CDS/bootstrap_support/profileNJ/X_topo/PROC_GENE_TREES_100/)
###         (DECOSTAR/DATA_DeCoSTAR/GENE_TREES/CDS/bootstrap_support/profileNJ/WG_topo/PROC_GENE_TREES_100/)
###      5- Support threshold to choose to contract or not a branch tree  
###         (Ex: 100)                                                     
###      6- Species position compare to the separator with gene ID prefix/postfix)
###          (prefix)                                                     
###                                                                       
###   OUTPUT:                                                             
###      - Contract weaks branches and resolve multifurcation             
###        by minimizing GDup and GLos                                    
###                                                                       
###   Name: 05b-profileNJ.sh                 Author: Yoann Anselmetti     
###   Creation date: 2016/01/06              Last modification: 2016/11/15
###

startTime=`date +%s`



separator="@"

if [ ! -d "$4" ]; then
	mkdir -p $4
fi

if ! [[ $6 == "prefix" || $6 == "postfix" ]]; then
    echo "ERROR: prefix/postfix boolean should be equal to \"prefix\" or \"postfix\" and not: $6" >&2;
    exit 1
fi

re='^[0-9]+$'
if ! [[ ${5} =~ $re ]]; then
   echo "ERROR: parameter 5 (minimum support for profileNJ) is not a number: ${5}" >&2;
   exit 1
else
   if (( "${5}" > "100" )); then 
      echo "ERROR: Threshold value: ${5} > 100" >&2;
      exit 1
   fi
fi

GF_ID=$(echo $3 | awk -F/ '{print $NF}' | sed 's/gene_tree_//g' | sed 's/.nwk//')

echo -e "Processing Gene family $GF_ID:"  
echo -e "\t$1 -s $2 -g $3 -o $4/${GF_ID}.nwk --spos $6 -r best --sep $separator --firstbest --cost 2 1 --seuil $5"
$1 -s $2 -g $3 -o $4/${GF_ID}.nwk --spos $6 -r best --sep $separator --firstbest --cost 2 1 --seuil $5



endTime=`date +%s`
runtime=$(($endTime-$startTime))

secondes=""
minutes=""
hours=$(($runtime/3600))
modH=$(($runtime%3600))
if [ $modH == 0 ]; then
    minutes="00"
    secondes="00"
else
	minutes=$(($modH/60))
	modM=$(($modH%60))
    if [ $modM == 0 ]; then
	    secondes="00"
	else
		secondes=$modM
	fi
fi

echo "Execution time: "$hours":"$minutes":"$secondes