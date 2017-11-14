#! /usr/bin/env python
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

from ete3 import Tree

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise


################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   GENE_file="data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF"
   GT_file="data/GENE_TREES/trees_DeCoSTAR_Xtopo.nwk"
   OUTPUT_annot_file="data/data_DeCoSTAR/GENE_file"
   sep="@"
   bool_gene="no"


   bool_GENE=False
   if (bool_gene=="Y" or bool_gene=="y" or bool_gene=="Yes" or bool_gene=="yes"):
      bool_GENE=True
   elif (bool_gene=="N" or bool_gene=="n" or bool_gene=="No" or bool_gene=="no"):
      bool_GENE=False
   else:
      exit("!!! ERROR, bool_gene (4th parameter) is incorrectly written: "+bool_gene+" should be \"yes\" or \"no\" !!!")

   OUTPUT_DIR=path.dirname(path.realpath(OUTPUT_annot_file))
   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   # STRUCTURE for gene edge and adjacency (EDGE==ADJ)
   GENE=namedtuple("GENE",["spe","ctg","gf","id","ori","start","end","exon_nb","exon_pos"])


############################################
### INDEXATION OF GENES INFOS BY GENE ID ###
############################################
   print "INDEX of gene info by gene ID ...",
   dict_ID_gene={}
   with open(GENE_file,'r') as ortho_file:
      for line in ortho_file:
         # If INPUT file: INPUT_DATA/ALL_species_GENE_file
         r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         if r:
            species=r.group(1)
            contig=r.group(2)
            gf_ID=r.group(3)
            gene_ID=r.group(4)
            orientation=r.group(5)
            start_pos=r.group(6)
            end_pos=r.group(7)
            exon_nb=r.group(8)
            exon_list=r.group(9)

            if species!="#species":
               gene_info=GENE(species,contig,gf_ID,gene_ID,orientation,start_pos,end_pos,exon_nb,exon_list)
               dict_ID_gene[gene_ID]=gene_info

         else:
            exit("\n!!! ERROR, in line:\n\t"+line+" of file "+GENE_file+" is not to the expected format\n!!!")
   ortho_file.close()
   print "DONE"



   ########################
   ### with ETE3 module ###
   ########################
   print "Browse gene trees file "+GT_file+" to filter gene annotation file:"
   i=0
   output_annot=open(OUTPUT_annot_file,"w")
   input_file=open(GT_file,"r")
   list_trees=list()
   # Browse gene trees file line by line
   for tree_line in input_file:
      list_genes=list()
      tree_str=tree_line.replace("\n","")
      # print tree_str
      tree=Tree(tree_str)
      # Get list of extant genes in current gene tree
      for spe_gene in tree.get_leaf_names():
         gene=spe_gene.split("@")[1]
         if gene in list_genes:
            exit("\n!!! ERROR: Gene "+gene+" is present in several trees in gene trees file: "+GT_file+" !!!")
         else:
            if gene in dict_ID_gene:
               list_genes.append(gene)
               output_annot.write(dict_ID_gene[gene].spe+"\t"+dict_ID_gene[gene].ctg+"\t"+dict_ID_gene[gene].gf+"\t"+dict_ID_gene[gene].id+"\t"+dict_ID_gene[gene].ori+"\t"+dict_ID_gene[gene].start+"\t"+dict_ID_gene[gene].end+"\t"+dict_ID_gene[gene].exon_nb+"\t"+dict_ID_gene[gene].exon_pos+"\n")
               dict_ID_gene.pop(gene,None)
            else:
               print "\tGene "+gene+" is not present in GENE file "+GENE_file+" !!!"

   print "\n\t=> There are "+str(len(dict_ID_gene))+" genes present in GENE file "+GENE_file+" but not present in gene trees file "+GT_file+" !!!\n"
   if bool_GENE:
      for gene in dict_ID_gene:
         output_annot.write(dict_ID_gene[gene].spe+"\t"+dict_ID_gene[gene].ctg+"\tNA\t"+dict_ID_gene[gene].id+"\t"+dict_ID_gene[gene].ori+"\t"+dict_ID_gene[gene].start+"\t"+dict_ID_gene[gene].end+"\t"+dict_ID_gene[gene].exon_nb+"\t"+dict_ID_gene[gene].exon_pos+"\n")
   dict_ID_gene.clear()
   output_annot.close()

   # Sort OUTPUT annotation gene file by species, then contigs and genes.
   output_file=open("buffer_file",'w')
   output_file.write("#species\tctg\tgene_family\tgene\tgene_orientation\tstart_gene\tend_gene\t#exons\texons_position\n")
   output_file.close()
   command_line="sort -k1d,1d -k2d,2d -k6n,6n "+OUTPUT_annot_file+" >> buffer_file; mv buffer_file "+OUTPUT_annot_file
   subprocess.call(command_line,shell=True)

   # Get and print execution time of the script
   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
