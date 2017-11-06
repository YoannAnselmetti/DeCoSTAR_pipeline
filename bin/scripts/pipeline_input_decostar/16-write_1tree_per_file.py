#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:                                                                         
###      Create FASTA file of all genes of given species                            
###   INPUT:                                                                        
###      1- INPUT gene trees file                                                   
###         (data/GENE_TREES/trees_DeCoSTAR_Xtopo.nwk)
###         (data/GENE_TREES/trees_DeCoSTAR_WGtopo.nwk)
###      2- Directory containing RAW gene trees of Rob Waterhouse (to name gene trees)
###         (data/INPUT_DATA/OG_CDS_newtrees)                     
###      3- OUTPUT directory where gene tree files will be stored                   
###         (data/data_DeCoSTAR/decostar/Xtopo_pNJ/DeCoSTAR_Anopheles_Xtopo_gene_trees)            
###         (data/data_DeCoSTAR/decostar/WGtopo/DeCoSTAR_Anopheles_WGtopo_gene_trees)           
###      4- Path written in gene trees file for DeCo*                               
###         (data/data_DeCoSTAR/decostar/Xtopo_pNJ/DeCoSTAR_Anopheles_Xtopo_gene_trees)                              
###         (data/data_DeCoSTAR/decostar/WGtopo/DeCoSTAR_Anopheles_WGtopo_gene_trees)                             
###      5- OUTPUT trees file for DeCo*                                             
###         (data/data_DeCoSTAR/decostar/Xtopo_pNJ/distrib_DeCoSTAR_Anopheles_Xtopo_gene_trees.txt)
###         (data/data_DeCoSTAR/decostar/WGtopo/distrib_DeCoSTAR_Anopheles_WGtopo_gene_trees.txt)
###      6- Character separator between species name and gene ID                    
###         (@)                                                                     
###      7- prefix/postfix boolean                                                  
###         (prefix or postfix)                                                     
###                                                                                 
###   OUTPUT:	(RUN in ~)                                                       
###      - OUTPUT FASTA file containing gene sequences of selected species          
###                                                                                 
###   Name: 16-write_1tree_per_file.py            Author: Yoann Anselmetti    
###   Creation date: 2016/09/02                   Last modification: 2016/11/16
###


from sys import argv
from os import close, path, makedirs, listdir
from datetime import datetime
from ete3 import Tree

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise


if __name__ == '__main__':

   start_time = datetime.now()

   input_file=open(argv[1],"r")
   GT_dir=argv[2]
   output_dir=argv[3]
   path_written=argv[4]
   output_trees_file=argv[5]
   separator="@"
   order_bool="prefix"

   OUTPUT_DIR=path.dirname(output_trees_file)
   # To be sure than directory have no "/" to the end of the path
   OUTPUT_DIR=path.normpath(OUTPUT_DIR)
   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   # To be sure than directory have no "/" to the end of the path
   output_dir=path.normpath(output_dir)
   # Create output_dir if not existing
   if not path.exists(output_dir):
      mkdir_p(output_dir)


   # Get list of file of RAW gene trees   
   raw_gene_trees = [f for f in listdir(GT_dir) if path.isfile(path.join(GT_dir,f))]

   # Store association gene ID with gene family ID in "dict_geneID_gfID"
   dict_geneID_gfID={}
   for gt in raw_gene_trees:
      tree_file=GT_dir+"/"+gt
      str_tree=""

      raw_tree_file=open(tree_file,'r')
      str_tree=raw_tree_file.read()

      if str_tree!="();":
         gfID=gt.split(".")[1]
         tree=Tree(tree_file)
         # Get list of extant genes in current gene tree
         for gene in tree.get_leaf_names():
            dict_geneID_gfID[gene]=gfID

   i=0
   # Read gene trees file to split in 1 tree / file with ID from raw gene trees
   output_trees=open(output_trees_file,"w")
   for tree_line in input_file:
      tree=Tree(tree_line)
      # Get list of extant genes in current gene tree
      for leaf in tree.get_leaf_names():
         # print leaf
         gene=""
         if separator in leaf:
            if order_bool=="prefix":
               gene=leaf.split(separator)[1]
            elif order_bool=="postfix":
               gene=leaf.split(separator)[0]
            else:
               exit("ERROR, parameter 7 should be equal to \"postfix\" or \"prefix\" !!!")
         else:
            gene=leaf
            # print gene
         if gene in dict_geneID_gfID:
            gfID=dict_geneID_gfID[gene]
            # print gfID
            i+=1
            file_name=gfID+".nwk"
            output_treefile=open(output_dir+"/"+file_name,"w")
            output_treefile.write(tree_line)
            output_treefile.close()
            output_trees.write(path_written+"/"+gfID+".nwk\n")
            break
         else:
            exit("ERROR, gene "+gene+" is NOT present in raw gene trees of directory "+GT_dir)
   input_file.close()
   output_trees.close()
   print "There are "+str(i)+" gene families in OUTPUT directory "+output_dir