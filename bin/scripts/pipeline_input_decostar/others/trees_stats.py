#! /usr/bin/env python2
# -*- coding: utf-8 -*-
###                                                                            
###   Goal:                                                                    
###      Read gene trees files and compute global statistics
###                                                                            
###   INPUT:                                                                   
###      1- Gene trees file                                                 
###         (data/GENE_TREES/trees_DeCoSTAR_Xtopo.nwk)
###                                                                            
###   OUTPUT:                                                                  
###      - Compute general statitics on a gene trees file     
###                                                                            
###   Name: trees_stats.py          Author: Yoann Anselmetti     
###   Creation date: 2018/05/28     Last modification: 2019/05/28
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
   gene_trees_file=argv[1]


   uniq_list_leaves_nb=set()
   input_file=open(gene_trees_file,"r")
   for tree_line in input_file:
      tree_str=tree_line.replace("\n","")
      # print tree_str
      tree=Tree(tree_str)

      leaves_nb=0
      for leaf in tree.get_leaf_names():
         leaves_nb+=1

      ### Store number of leaves in current tree
      uniq_list_leaves_nb.add(leaves_nb)

   input_file.close()

   print sorted(uniq_list_leaves_nb)


   # Get and print execution time of the script
   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
