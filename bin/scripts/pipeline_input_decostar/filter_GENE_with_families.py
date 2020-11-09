#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###                                                                       
###   Goal:                                                               
###      Discard genes that are not present in gene trees/families.
###      Step done before to detect included genes => Reduce the number of included genes that will be discarded
###                                                                       
###   INPUT:                                                              
###      1- Gene TREES/FAMILIES file                                   
###         (data/INPUT_DATA/unrooted_raw_trees.nwk)    
###      2- INPUT directory containing sorted GENE files                  
###         (data/GFF_to_GENE_files/sorted_GENE) 
###      3- OUTPUT directory path where results will be stored            
###         (data/GFF_to_GENE_files/filtered_GENE)
###      4- Separator between species name and gene ID
###         (@)
###
###   OUTPUT:                                                             
###      - OUTPUT directory containing sorted and filtered genes file for each species                                               
###                                                                       
###   Name: filter_GENE_with_families.py     Author: Yoann Anselmetti
###   Creation date: 2016/09/09              Last modification: 2020/11/05
###                                                                       

from sys import argv, exit
from re import search
from os import close, path, makedirs, listdir
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


def get_list_genes_in_GF(dict_species_geneList,input_file,sep):
   list_genes=list()
   bool_file=""
   for line in input_file:
      if not bool_file:
         bool_file="OK"
         r=search("^\(",line)
         if r:
            bool_GT=True
            print("\t=> File containing gene cluster is a gene TREES file")
         else:
            bool_GT=False
            print("\t=> File containing gene cluster is a gene FAMILIES file")

      # If file is a newick, NHX or tree file 
      if bool_GT:
         # Store gene trees in ETE3 "Tree object"
         tree_str=line.replace("\n","")
         tree=Tree(tree_str)

         # For each leaf in current tree: extract species name and Gene_ID
         for gene in tree.get_leaf_names(): 
            # print(gene)           
            species=gene.split(sep)[0]
            gene_ID=gene.split(sep)[1]

            # Store Gene_ID in dict_species_geneList
            if not species in dict_species_geneList:
               dict_species_geneList[species]=list()
            dict_species_geneList[species].append(gene_ID)

      else:
         # Test if the file is Tab-separated with 2 columns (=> Gene family file format) 
         r=search("^([^\t]*)\t([^\t\n]*)\n$",line)
         if r:
            GF_ID=r.group(1)
            gene=r.group(2)

            species=gene.split(sep)[0]
            gene_ID=gene.split(sep)[1]

            # Store Gene_ID in dict_species_geneList
            if not species in dict_species_geneList:
               dict_species_geneList[species]=list()
            dict_species_geneList[species].append(gene_ID)
         else:
            exit("ERROR: Wrong format for gene cluster => Should de a Tab-separated files of 2 columns: c1: gene family ID, c2: gene ($(species_name)$sep$(gene_ID))!!!")

         


if __name__ == '__main__':

   start_time = datetime.now()

   GF_file=argv[1]
   INPUT_dir=argv[2]
   OUTPUT_dir=argv[3]
   sep=argv[4]

   # Create OUTPUT_dir if not existing
   mkdir_p(OUTPUT_dir)

   #######
   ### STORE GENES PRESENT IN GENE TREES/FAMILIES / SPECIES NAME
   #######
   print("1/ STORE gene present in gene families/trees:")
   dict_species_geneList=dict()
   input_file=open(GF_file,"r")
   get_list_genes_in_GF(dict_species_geneList,input_file,sep)
   input_file.close()


   #######
   ### FILTER GENES PRESENT IN GENE FILE THAT ARE NOT PRESENT IN GENE TREES/FAMILIES
   #######
   # Get list of GENE files present in INPUT directory
   list_files=listdir(INPUT_dir)

   print("2/ Filter genes that are present not present in gene families/trees considered:") 
   dict_species_geneFilt=dict()
   # Browse GENE files contained in INPUT directory
   for gene_file in sorted(list_files):
      # print(gene_file)
      r_spe=search('^(.*)_sorted.txt$', gene_file)
      if r_spe:
         name_spe=r_spe.group(1)
         print("\t"+name_spe)

         output_file=open(OUTPUT_dir+"/"+name_spe+"_filtered.txt","w")
         input_file=open(INPUT_dir+"/"+gene_file)
         # Browse sorted GENE file of "name_spe"
         output_file.write("#species\tctg/scaff\tgene\torientation_gene\tstart_gene\tend_gene\t#exons\texons_position\n")
         for line in input_file:
            r=search('^([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([+-])[\t ]*([0-9]*)[\t ]*([0-9]*)[\t ]*([0-9]*)[\t ]*([:0-9-]*)\n$', line)
            if r:
               # print(line)
               species=r.group(1)
               contig=r.group(2)
               gene=r.group(3)
               orientation=r.group(4)
               start=r.group(5)
               stop=r.group(6)
               nb_exon=r.group(7)
               exon_pos=r.group(8)

               if species!="#species":
                  # print(species+" "+gene)
                  if gene in dict_species_geneList[species]:
                     output_file.write(line)
                     dict_species_geneList[species].remove(gene)
                  else:
                     if not species in dict_species_geneFilt:
                        dict_species_geneFilt[species]=list()
                     dict_species_geneFilt[species].append(gene)

   print("\nNumber of genes filtered for each species:")
   for species in sorted(dict_species_geneFilt):
      print("\t- "+species+":\t"+str(len(dict_species_geneFilt[species]))+" genes filtered")

   print("\nNumber of genes present in gene trees/families but not in GFF for each species (should be equal to 0):")
   for species in sorted(dict_species_geneList):
      print("\t- "+species+":\t"+str(len(dict_species_geneList[species]))+" genes in gene trees/families but not in GFF")


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
