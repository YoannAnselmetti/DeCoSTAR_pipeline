#!/bin/bash
###                                                                                      
###   Goal:                                                                              
###      Filter gene trees containing included genes                                     
###                                                                                      
###   INPUT:                                                                             
###      1- INPUT gene trees file to filter with included genes                          
###         (data/INPUT_DATA/unrooted_raw_trees.nwk)        
###      2- Inclusion gene file                                                          
###         (data/GFF_to_GENE_files/with_filter/ALL_species_Inclusion_file)
###      3- OUTPUT gene trees file                                                         
###         (data/GENE_TREES/unrooted_trees_filtered.nwk)   
###                                                                                      
###   OUTPUT:                                                                            
###      - Gene trees file without included gene with species name corresponding to the gene ID                                                                
###                                                                                      
###   Name: filter_geneTrees_with_includedGenes.sh     Author: Yoann Anselmetti       
###   Creation date: 2015/09/22                        Last modification: 2017/10/18  
###                                                                                      

##########
#### REMOVED GENE TREES CONTAINING INCLUDED GENES
##########

# Get list of gene included in another gene from inclusion gene file ($2)
echo -e -n "\tGet list of genes inluded in others genes ..."
cut -f2 $2 > $2"_buffer"
echo "DONE"

# Removed gene trees containing included genes and store it in OUTPUT gene trees file $3
echo -e -n "\tDelete gene trees containing included genes ..."
grep -Fv -f $2"_buffer" $1 > $3
echo "DONE"

rm $2"_buffer"

# # Command-line to see difference between gene trees with filter and without filter and have size distribution of gene trees present in "unrooted_trees-ASTEI+filt.nwk" but not in "unrooted_trees-ASTEI-filt.nwk":
# diff data/GENE_TREES/UNROOTED/unrooted_trees-ASTEI-filt.nwk data/GENE_TREES/UNROOTED/unrooted_trees-ASTEI+filt.nwk |grep \> |grep -o -n ":"| cut -d : -f 1 | uniq -c|awk -F" " '{print $(NF-1)}'|sort -n|uniq -c