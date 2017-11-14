#!/usr/bin/env python
# -*- coding: utf-8 -*-
###                                                                       
###   Goal:                                                               
###      Create gene tree file and distance matrix file                   
###      to root gene tree with profileNJ from RAxML gene trees           
###                                                                       
###   INPUT:                                                              
###      1- INPUT DIRECTORY containing files Best ML tree compute by RAxML
###         (DATA/DATA_WaterHouse/GENE_TREES/PROT/bootstrap/RAxML/GF0000001)
###      2- OUTPUT directory path where gene trees files will be stored   
###         (DATA/DATA_WaterHouse/GENE_TREES/PROT/bootstrap/profileNJ/UNROOTED_GENE_TREES)                                          
###      3- Bootstrap or SH-like support compute with RAxML boolean       
###         (bootstrap | SH-like)                                         
###      4- gene_TAG-species_name association file                         
###         (WORKING_PROJECT/fixed_params_160407/INPUT_DATA/name_geneID_18anopheles)
###      6- prefix/postfix boolean                                        
###         (prefix or postfix)                                           
###                                                                       
###   OUTPUT:                                                             
###      - Gene trees files and distance matrix files for each gene tree  
###        in gene trees file use as INPUT                                
###                                                                       
###   Name: 05a-modifGT.py                   Author: Yoann Anselmetti
###   Creation date: 2015/09/21              Last modification: 2017/03/29
###



from sys import argv, path
from re import search, sub, match
from os import close, path, makedirs, rename, listdir
from datetime import datetime

from Bio import Phylo   # Biopython library
from ete3 import Tree


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise

def file_len(fname):
   with open(fname) as f:
      for i, l in enumerate(f):
         pass
   return i+1


def replaceStringInFile(filePath,old_pattern,new_pattern):
   "replaces all string by a regex substitution"
   tempName = filePath+'~~~'
   inputFile = open(filePath)
   outputFile = open(tempName,'w')
   fContent = unicode(inputFile.read(), "utf-8")

   outText = sub(old_pattern, new_pattern, fContent)
   outputFile.write((outText.encode("utf-8")))

   outputFile.close()
   inputFile.close()

   rename(tempName, filePath)


################
###   MAIN
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   DIR_RAxML=argv[1]
   DIR_gene_trees=argv[2]
   bool_supp_raxml=argv[3]
   speciesName_geneID_file=open(argv[4],"r")
   order_bool=argv[5]

   sep="@"

   # Store association Gene-Species in "dict_gene_species"
   dict_geneID_speciesName={}
   for gene in speciesName_geneID_file:
      r=search("^([^\t ]*)[\t ]*([^\t ]*)\n$",gene)
      if r:
         speciesName=r.group(1)
         geneID=r.group(2)
         dict_geneID_speciesName[geneID]=speciesName
   speciesName_geneID_file.close()


   # To be sure than directory have no "/" to the end of the path
   DIR_RAxML=path.normpath(DIR_RAxML)
   # To be sure than directory have no "/" to the end of the path
   DIR_gene_trees=path.normpath(DIR_gene_trees)


   # Create DIR_gene_trees if not existing
   if not path.exists(DIR_gene_trees):
      mkdir_p(DIR_gene_trees)


   # Get Gene Family ID 
   GF_ID=path.basename(DIR_RAxML)

   # Get path of Best ML tree with support (SH-like or bootstrap)
   BestML_tree=""
   if bool_supp_raxml=="bootstrap":
      BestML_tree="RAxML_bipartitionsBranchLabels."+GF_ID+"_treeBS"
   elif bool_supp_raxml=="SH-like":
      BestML_tree="RAxML_fastTreeSH_Support."+GF_ID+"_treeSH"
   else:
      exit("\n!!! ERROR, RAxML support boolean, should be \"bootstrap\" or \"SH-like\" and not "+bool_supp_raxml+" !!!")

   # TEST if a Best ML tree have been computed for the current Gene family ID
   if not BestML_tree in listdir(DIR_RAxML):
      print "NO BestML Tree file for Gene Family: "+GF_ID
   else:
      print "Processing Gene Family "+GF_ID+":"
      # Path for OUTPUT gene tree & matrix distance files
      # OUTPUT_DIST=DIR_dist_matrix+"/dist_matrix_"+GF_ID+".dist"
      OUTPUT_tree=DIR_gene_trees+"/gene_tree_"+GF_ID+".nwk"

      # Get tree with RAxML support
      INPUT_tree=DIR_RAxML+"/"+BestML_tree

      # Write INPUT_tree (Best ML RAxML tree) without branch length in OUTPUT_tree
      in_tree=open(INPUT_tree,'r')
      out_tree=open(OUTPUT_tree,'w')
      # Remove branch length from INPUT RAxML gene tree
      tree_with_ete3_format = sub('\)(:[0-9.]+)\[([0-9]+)\]', r')\2\1', in_tree.read())
      # Add species name to gene 
      out_tree.write(tree_with_ete3_format)
      out_tree.close()
      in_tree.close()


   #########################
   ### with ETE3 module  ### => Function to prune gene trees!!!
   #########################

      list_genes=list()
      tree=Tree(OUTPUT_tree)
      # Get list of extant genes in current gene tree
      list_genes=tree.get_leaf_names()



      # Browse list_genes to create FASTA file of the current gene tree
      for gene in list_genes:
         print "\t"+gene
         gene_present=False
         species=""
         for geneID in dict_geneID_speciesName:
            if geneID in gene:
               gene_present=True
               species=dict_geneID_speciesName[geneID]
               break
         if not gene_present:
            exit("ERROR, species of gene "+gene+" is not present in "+str(argv[4])+" !!!!")
         else:
            old_pattern=gene
            new_pattern=""
            if order_bool=="prefix":
               new_pattern=species+sep+gene
            elif order_bool=="postfix":
               new_pattern=gene+sep+species
            else:
               exit("ERROR, parameter 5 should be equal to \"prefix\" or \"postfix\" !!!")
            replaceStringInFile(OUTPUT_tree,old_pattern,new_pattern)


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
