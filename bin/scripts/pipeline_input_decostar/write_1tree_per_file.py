#! /usr/bin/env python2
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create FASTA file of all genes of given species
###   INPUT:
###      1- INPUT RAW gene trees file filtered with species name (-ASTEI+filt+spe)
###         (data/GENE_TREES/XXX.nwk)
###      2- GENE file contained in gene trees used for Anopheles dataset
###         (data/data_DeCoSTAR/GENE_file)
###      3- OUTPUT directory where gene tree files will be stored
###         (data/RAW_Anopheles_gene_trees)
###      4- Path written in gene trees file for DeCo*
###         (data/RAW_Anopheles_gene_trees)
###      5- OUTPUT trees file for DeCo*
###         (data/distrib_RAW_Anopheles_gene_trees.txt)
###      6- Character separator between species name and gene ID
###         (@)
###
###   OUTPUT:	(RUN in less than 1 minute)
###      - Gene trees files used as INPUT by DeCoST
###
###   Name: write_1tree_per_file.py    Author: Yoann Anselmetti
###   Creation date: 2016/09/02        Last modification: 2017/11/07
###


from sys import argv
from os import close, path, makedirs, listdir
from re import search
from datetime import datetime
from ete3 import Tree
import errno

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else:
         raise


if __name__ == '__main__':

   start_time = datetime.now()

   input_file=open(argv[1],"r")
   GENE_file=argv[2]
   output_dir=argv[3]
   path_written=argv[4]
   output_trees_file=argv[5]
   separator=argv[6]

   OUTPUT_DIR=path.dirname(output_trees_file)
   # To be sure than directory have no "/" to the end of the path
   OUTPUT_DIR=path.normpath(OUTPUT_DIR)
   # Create OUTPUT_DIR if not existing
   mkdir_p(OUTPUT_DIR)

   # To be sure than directory have no "/" to the end of the path
   output_dir=path.normpath(output_dir)
   # Create output_dir if not existing
   mkdir_p(output_dir)


   # Store association gene ID with gene family ID in "dict_geneID_gfID"
   dict_geneID_gfID={}
   input_gene=open(GENE_file,"r")
   for gene in input_gene:
      r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$', gene)
      if r:
         name_species=r.group(1)
         contig_ID=r.group(2)
         gf_ID=r.group(3)
         gene_ID=r.group(4)
         orientation=r.group(5)
         start=r.group(6)
         stop=r.group(7)
         exon_nb=r.group(8)
         exon_pos=r.group(9)

         if name_species!="#species":
            dict_geneID_gfID[gene_ID]=gf_ID
   input_gene.close()


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
            gene=leaf.split(separator)[1]
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
         # else:
         #    print "Gene family of gene "+gene+" is no more present after gene trees inference pipeline"
         break

   input_file.close()
   output_trees.close()
   print "There are "+str(i)+" gene families in OUTPUT directory "+output_dir
