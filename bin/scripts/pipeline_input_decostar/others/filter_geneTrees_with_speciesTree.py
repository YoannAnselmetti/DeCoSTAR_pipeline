#! /usr/bin/env python3
# -*- coding: utf-8 -*-
###
###   Goal:
###      Remove genes from gene trees taht don't corresponds to species present in species tree
###
###   INPUT:
###      1- Species tree file
###      2- INPUT Gene trees file
###      3- OUTPUT Gene trees file
###      4- separator used in gene ID between species name and gene ID
###
###   OUTPUT:                                                                  
###      - New gene trees file limited to species present in species tree
###
###   Name: filter_geneTrees_with_speciesTrees.py      Author: Yoann Anselmetti     
###   Creation date: 2015/11/11                        Last modification: 2020/11/05
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

   print(list_species)

   input_file=open(GT_infile,"r")
   output_file=open(GT_outfile,"w")
   for tree_line in input_file:
      tree_str=tree_line.replace("\n","")
      # print(tree_str)
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
