#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###                                                                       
###   Goal:                                                               
###      To add gene family ID to the GENE file containing all species of the dataset
###                                                                           
###   INPUT:                                                              
###      1- INPUT GENE file for ALL species                               
###         (data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file)
###      2- Gene trees file || Directory containing 1 file/gene tree || Gene family file
###         (data/INPUT_DATA/OG_CDS_newtrees)           
###      3- OUTPUT gene file with gene families ID                                        
###         (data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF)
###      4- Separator between species name and gene ID
###         (@)
###                                                                       
###   OUTPUT:                                                             
###      - GENE file with gene family ID               
###                                                                       
###   Name: add_geneFamilyID.py           Author: Yoann Anselmetti    
###   Creation date: 2016/09/09           Last modification: 2020/11/05
###                                                                       

from sys import argv
from re import search
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

def get_GF_ID(i):
   GF_ID=""
   i+=1
   if i>=1000000 and i<=9999999:
      GF_ID="GF"+str(i)
   elif i>=100000:
      GF_ID="GF0"+str(i)
   elif i>=10000:
      GF_ID="GF00"+str(i)
   elif i>=1000:
      GF_ID="GF000"+str(i)
   elif i>=100:
      GF_ID="GF0000"+str(i)
   elif i>=10:
      GF_ID="GF00000"+str(i)
   elif i>=1:
      GF_ID="GF000000"+str(i)
   else:
      exit("ERROR too much gene families!!!")
   return i,GF_ID


################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   INPUT_gene_file=open(argv[1],"r")
   GT=argv[2]
   OUTPUT_gene_file=argv[3]
   sep=argv[4]

   # To be sure than directory have no "/" to the end of the path
   OUTPUT_dir=path.dirname(OUTPUT_gene_file)

   # Remove OUTPUT_dir if existing
   if not path.exists(OUTPUT_dir):
      mkdir_p(OUTPUT_dir)


   dict_geneID_gfID=dict()
   if path.isdir(GT):
      # Get list of file of RAW gene trees 
      gene_trees = [f for f in listdir(GT) if path.isfile(path.join(GT,f))]

      # Store association gene ID with gene family ID in "dict_geneID_gfID" from raw gene trees file
      for gt in gene_trees:
         GT_file=GT+"/"+gt
         gene_tree_file=open(GT_file, 'r')

         str_tree=gene_tree_file.read()
         if str_tree!="();":
            gfID=gt.split(".")[1]
            tree=Tree(GT_file)
            # Get list of extant genes in current gene tree
            for gene in tree.get_leaf_names():
               dict_geneID_gfID[gene]=gfID

   elif path.isfile(GT):
      input_file=open(GT,"r")
      bool_file=""
      i=0
      for line in input_file:
         if not bool_file:
            r=search("^\(.*\);\n$",line)
            if r:
               bool_file="GT"
               print("File containing gene cluster is a gene TREES file")
            else:
               bool_file="GF"
               print("File containing gene cluster is a gene FAMILIES file")

         if bool_file=="GF":
            r=search("^([^\t]*)\t([^\t\n]*)\n$",line)
            if r:
               gene_family=r.group(1)
               gene=r.group(2)
               dict_geneID_gfID[gene]=gene_family

            else:
               exit("ERROR: wrong format for gene cluster!!!")

         elif bool_file=="GT":
            tree_str=line.replace("\n","")
            # print(tree_str)
            tree=Tree(tree_str)
            i,gfID=get_GF_ID(i)
            for gene in tree.get_leaf_names():
               dict_geneID_gfID[gene]=gfID

   else:
      exit("ERROR, parameter 2:\n\t"+GT+"\nshould correspond to a gene trees or families file or a directory containing gene tree files!!!")



   # Read GENE file and modify it
   output_file=open(OUTPUT_gene_file,"w")
   output_file.write("#species\tctg/scaff/chr\tgene_family\tgene\torientation_gene\tstart_gene\tend_gene\t#exons\texons_position\n")
   for line in INPUT_gene_file:
      r=search('^([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t\n]*)\n$', line)
      if r:
         # print(line)
         species=r.group(1)
         contig=r.group(2)
         gene=r.group(3)
         orientation=r.group(4)
         start_old=r.group(5)
         stop_old=r.group(6)
         nb_exon_old=r.group(7)
         exon_pos_old=r.group(8)
         start_new=r.group(9)
         stop_new=r.group(10)
         nb_exon_new=r.group(11)
         exon_pos_new=r.group(12)

         geneID=species+sep+gene

         if species!="#species":
            gfID=""
            if geneID in dict_geneID_gfID:
               gfID=dict_geneID_gfID[geneID]
            else:
               print("WARNING: gene "+gene+" is not present in gene families/trees!!!")
               gfID="NA"

            # if start_old!=start_new or stop_old!=stop_new:
            #    print(line)

            output_file.write(species+"\t"+contig+"\t"+gfID+"\t"+gene+"\t"+orientation+"\t"+start_new+"\t"+stop_new+"\t"+nb_exon_new+"\t"+exon_pos_new+"\n")

      else:
         exit("ERROR, line:\n\t"+line+"\t of file "+argv[1]+" is incorrectly written!!!")


   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))