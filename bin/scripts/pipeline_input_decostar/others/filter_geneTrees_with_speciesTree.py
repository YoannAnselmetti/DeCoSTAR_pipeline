#! /usr/bin/env python2
# -*- coding: utf-8 -*-
###                                                                            
###   Goal:                                                                    
###      parse and filter GENE file with gene trees file (and rewrite GENE family ID) to obtain a file with gene info for ADseq instance creation                                                     
###                                                                            
###   INPUT:                                                                   
###      1- INPUT GENE file => OUTPUT file of GFF->GENE pipeline               
###         (data/GFF_to_GENE_files/CDS/with_filter/ALL_species_GENE_file_with_GF)
###      2- Gene trees file                                                    
###         (data/GENE_TREES/trees_DeCoSTAR_Xtopo.nwk)
###         (data/GENE_TREES/trees_DeCoSTAR_WGtopo.nwk)
###      3- OUTPUT annotation gene file                                        
###         (data/data_DeCoSTAR/GENE_file)                           
###      4- separator used in gene trees between species name and gene ID       
###         (@)                                                                
###      5- Boolean to know if we have to write gene that are not in gene trees 
###         (Y/y/Yes/yes: if want to write gene not present in gene trees)     
###         (N/n/No/no: if DON'T want to write gene not present in gene trees) 
###                                                                            
###   OUTPUT:                                                                  
###      - Create annotation gene file for ARt-DeCo_seq instance creation      
###                                                                            
###   Name: filter_GENE_with_geneTrees.py      Author: Yoann Anselmetti     
###   Creation date: 2015/11/11                   Last modification: 2017/11/07
###                                                                            

from sys import argv
from re import search, match
from os import close, path, makedirs
from collections import namedtuple   #New in version 2.6
import subprocess
from datetime import datetime
import errno

from ete3 import Tree

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else:
         raise

################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   species_file=argv[1]
   GT_infile=argv[2]
   GT_outfile=argv[3]
   sep=argv[4]


   species_tree=Tree(species_file)
   list_species=species_tree.get_leaf_names()

   print list_species

   input_file=open(GT_infile,"r")
   output_file=open(GT_outfile,"w")
   for tree_line in input_file:
      tree_str=tree_line.replace("\n","")
      # print tree_str
      tree=Tree(tree_str)
      list_genes=list()
      for spe_gene in tree.get_leaf_names():
         species=spe_gene.split(sep)[0]
         if species in list_species:
            list_genes.append(spe_gene)
      if len(list_genes)>1:
         tree.prune(list_genes)
         output_file.write(tree.write(format=9)+"\n")
   input_file.close()
   output_file.close()

   # Get and print execution time of the script
   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
