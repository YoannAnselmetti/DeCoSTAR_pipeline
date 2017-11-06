#!/usr/bin/env python
# -*- coding: utf-8 -*-
###                                                                       
###   Goal:                                                               
###      Discard genes that are niot present in gene trees/families (before to detect included genes to limit the number of included genes that will be discarded)
###                                                                       
###   INPUT:                                                              
###      1- Gene TREES/FAMILIES file                                   
###         (data/INPUT_DATA/unrooted_raw_trees.nwk)
###      2- gene_TAG-species_name association file                        
###         (data/INPUT_DATA/name_geneID_18anopheles)     
###      3- INPUT directory containing sorted GENE files                  
###         (data/GFF_to_GENE_files/sorted_GENE) 
###      4- OUTPUT directory path where results will be stored            
###         (data/GFF_to_GENE_files/filtered_GENE)
###   OUTPUT:                                                             
###      - OUTPUT directory containing sorted and filtered genes file for each species                                               
###                                                                       
###   Name: filter_GENE_with_families.py     Author: Yoann Anselmetti
###   Creation date: 2016/09/09              Last modification: 2017/10/18
###                                                                       


from sys import argv, exit
from re import search, match
from os import close, path, makedirs, listdir, mkdir
from datetime import datetime
import errno
from ete3 import Tree

def uniq(seq):
   # Not order preserving
   keys = {}
   for e in seq:
       keys[e] = 1
   return keys.keys()

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise

def mean(table):
    return sum(table, 0.0) / len(table)

def variance(table):
    m=mean(table)
    return mean([(x-m)**2 for x in table])

def SD(table):
   return variance(table)**0.5


def get_list_genes_in_GF(dict_species_geneList,input_file):
   list_genes=list()
   bool_file=""
   for line in input_file:
      if not bool_file:
         r=search("^\(",line)
         if r:
            bool_file="GT"
            print "File containing gene cluster is a gene TREES file"
         else:
            bool_file="GF"
            print "File containing gene cluster is a gene FAMILIES file"

      if bool_file=="GF":
         r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         if r:
            ODBMOZ2_Level=r.group(1)
            ODBMOZ2_OG_ID=r.group(2)
            Protein_id=r.group(3)
            Gene_ID=r.group(4)
            Organism=r.group(5)
            UniProt_Species=r.group(6)
            UniProt_ACC=r.group(7)
            UniProt_Description=r.group(8)
            InterPro_domains=r.group(9)

            if ODBMOZ2_Level!="ODBMOZ2_Level":
               species=""
               gene_present=False
               # Check if species of the gene is present in "speciesName_geneID_file". Else we have to prune current gene tree to filter this gene
               for geneID in dict_geneID_speciesName:
                  if geneID in Gene_ID:
                     gene_present=True
                     species=dict_geneID_speciesName[geneID]
                     break

               if gene_present:
                  if not species in dict_species_geneList:
                     dict_species_geneList[species]=list()
                  dict_species_geneList[species].append(Gene_ID)
         else:
            exit("ERROR: wrong format for gene cluster!!!")

      elif bool_file=="GT":
         tree_str=line.replace("\n","")
         # print tree_str
         tree=Tree(tree_str)

         for gene in tree.get_leaf_names():
            gene_present=False
            for geneID in dict_geneID_speciesName:
               if geneID in gene:
                  gene_present=True
                  species=dict_geneID_speciesName[geneID]
                  break

            if gene_present:
               if not species in dict_species_geneList:
                  dict_species_geneList[species]=list()
               dict_species_geneList[species].append(gene)


if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   GF_file=argv[1]
   speciesName_geneID_file=open(argv[2],"r")
   INPUT_dir=argv[3]
   OUTPUT_dir=argv[4]


   # To be sure than directory have no "/" to the end of the path
   OUTPUT_dir=path.normpath(OUTPUT_dir)

   # To be sure than directory have no "/" to the end of the path
   INPUT_dir=path.normpath(INPUT_dir)

   # Create OUTPUT_dir if not existing
   if not path.exists(OUTPUT_dir):
      mkdir(OUTPUT_dir)


   # Store association Gene-Species in "dict_gene_species"
   dict_geneID_speciesName={}
   for gene in speciesName_geneID_file:
      r=search("^([^\t ]*)[\t ]*([^\t ]*)\n$",gene)
      if r:
         speciesName=r.group(1)
         geneID=r.group(2)
         dict_geneID_speciesName[geneID]=speciesName
   speciesName_geneID_file.close()


#################################################################
### STORE GENES PRESENT IN GENE TREES/FAMILIES / SPECIES NAME ###
#################################################################
   dict_species_geneList=dict()
   input_file=open(GF_file,"r")
   get_list_genes_in_GF(dict_species_geneList,input_file)
   input_file.close()


#########################################################
### GET LIST OF GENE FILES PRESENT IN INPUT DIRECTORY ###
#########################################################
   list_files=listdir(INPUT_dir)


#####################################################################################
### FILTER GENES PRESENT IN GENE FILE TAHT ARE NOT PRESENT IN GENE TREES/FAMILIES ###
#####################################################################################
   dict_species_geneFiltered=dict()
   # Browse GENE files contained in INPUT directory
   for gene_file in sorted(list_files):
      # print gene_file
      r_spe=search('^(.*)_sorted.txt$', gene_file)
      if r_spe:
         name_spe=r_spe.group(1)
         print name_spe

         output_file=open(OUTPUT_dir+"/"+name_spe+"_filtered.txt","w")
         input_file=open(INPUT_dir+"/"+gene_file)
         # Browse sorted GENE file of "name_spe"
         output_file.write("#species\tctg/scaff\tgene\torientation_gene\tstart_gene\tend_gene\t#exons\texons_position\n")
         for line in input_file:
            r=search('^([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([+-])[\t ]*([0-9]*)[\t ]*([0-9]*)[\t ]*([0-9]*)[\t ]*([:0-9-]*)\n$', line)
            if r:
               # print line
               species=r.group(1)
               contig=r.group(2)
               gene=r.group(3)
               orientation=r.group(4)
               start=r.group(5)
               stop=r.group(6)
               nb_exon=r.group(7)
               exon_pos=r.group(8)

               if species!="#species":
                  if gene in dict_species_geneList[species]:
                     output_file.write(line)
                     dict_species_geneList[species].remove(gene)
                  else:
                     if not species in dict_species_geneFiltered:
                        dict_species_geneFiltered[species]=list()
                     dict_species_geneFiltered[species].append(gene)

   print "\nNumber of genes filtered for each species:"
   for species in sorted(dict_species_geneFiltered):
      print "\t- "+species+":\t"+str(len(dict_species_geneFiltered[species]))+" genes filtered"


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
