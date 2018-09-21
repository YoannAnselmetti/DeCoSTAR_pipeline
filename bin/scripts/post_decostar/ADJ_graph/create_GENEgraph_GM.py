#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal
###      Create adjacencies graph file where node corresponds to gene ID
###      from a file containing assignment of part of contigs/scaffolds of species on its chromosomes
###      (RHmap, linkage map, genetic map, chromosome map, pairwise alignment, ...)
###
###   INPUT:
###      1- File containing scaffolds assignment to chromosome (AGP file)
###         (9passeriformes_dataset/data/INPUT_DATA/ADJ_graph/Zbor_scaff_anchored_in_Tgut_chr_050618_filt.agp)
###      2- CTG file
###         (9passeriformes_dataset/data/data_DeCoSTAR/CTG_file)
###      3- Species corresponding to the INPUT DOT file contig graph
###         (Zosterops_borbonicus)
###      4- OUTPUT file where Contig graph will be stored
###         (9passeriformes_dataset/results/ADJ_graph/DOT/GENE/Zosterops_borbonicus/Zbor_CTG_filt.dot)
###      5- OUTPUT SVG file where Contig graph will be stored
###         (9passeriformes_dataset/results/ADJ_graph/SVG/GENE/Zosterops_borbonicus/Zbor_CTG_filt.svg)
###      6- Boolean to know if we want CTG links
###         (Y/y: Yes  |  N/n: No)
###
###   OUTPUT:
###      Adjacencies graph file where node corresponds to gene ID
###      (annotated black edge corresponds to contig/scaffold with their ID)
###
###   Name: create_GENEgraph_GM.py                Author: Yoann Anselmetti
###   Creation date: 2016/03/07                   Last modification: 2016/06/28
###

from sys import argv, stdout
from re import search
from os import close, path, makedirs
from datetime import datetime
from collections import namedtuple   #New in version 2.6
import errno
import pygraphviz as pgv
import subprocess

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
   RHmap_file=argv[1]
   CTG_file=argv[2]
   INPUT_species=argv[3]
   OUTPUT_dot=argv[4]
   OUTPUT_svg=argv[5]
   bool_links=argv[6]

   LINKS=False
   if bool_links=="n" or bool_links=="N":
      LINKS=False
   elif bool_links=="y" or bool_links=="Y":
      LINKS=True
   else:
      exit("!!! ERROR, 7th parameter should be y/y or N/n and not "+bool_links+"!!!")


   OUTPUT_DIR=path.dirname(path.realpath(OUTPUT_dot))
   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   # To be sure than directory have no "/" to the end of the path
   OUTPUT_DIR=path.normpath(OUTPUT_DIR)

   OUTPUT_DIR2=path.dirname(path.realpath(OUTPUT_svg))
   # Create OUTPUT_DIR2 if not existing
   if not path.exists(OUTPUT_DIR2):
      mkdir_p(OUTPUT_DIR2)

   # To be sure than directory have no "/" to the end of the path
   OUTPUT_DIR2=path.normpath(OUTPUT_DIR2)

   CTG=namedtuple("CTG",["g1","g2"])

   color_list=("blue","red","green","purple","gold","orangered","chocolate","blueviolet","indigo","fuchsia","crimson","darkgreen","chartreuse4","darkblue","salmon","darkturquoise","cyan","tomato","teal","orchid","maroon","lime","deeppink","darkorange","dodgerblue","coral","sienna","olive","goldenrod","darkcyan","mediumvioletred","darkred","darkgrey","brown","darkgoldenrod","darkolivegreen","darkseagreen","darkslategrey","darkviolet","darkslateblue","orange","yellow","olivedrab")
   color_size=len(color_list)

#####################################################
### GET GENE ANNOTATION AND STORE IT IN dict_gene ###
#####################################################
   dict_CTG=dict()
   # Browse contigs extremities file and store ctg in "dict_spe_ctg" by species name key
   ctg_file=open(CTG_file,'r')
   for line in ctg_file:
      r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r:
         species=r.group(1)
         ctg=r.group(2)
         size=r.group(3)
         gene_nb=r.group(4)
         gf1=r.group(5)
         g1=r.group(6)
         ori1=r.group(7)
         start_g1=r.group(8)
         gf2=r.group(9)
         g2=r.group(10)
         ori2=r.group(11)
         stop_g2=r.group(12)

         if species==INPUT_species:
            dict_CTG[ctg]=CTG(g1,g2)




####################################################
### BROWSE RH MAP TO WRITE CTG ADJACENCIES GRAPH ###
####################################################
   graph=pgv.AGraph(strict=False)
   graph.node_attr['penwidth']=5
   graph.edge_attr['penwidth']=5
   graph.node_attr['fontsize']=20
   graph.edge_attr['fontsize']=20
   graph.node_attr['fontcolor']="black"
   graph.edge_attr['fontcolor']="black"
   graph.node_attr['fontname']="times-bold"
   graph.edge_attr['fontname']="times-bold"
   # Browse new extant adjacencies file and store adj of selected species in "list_ADJ"
   map_file=open(RHmap_file,'r')
   chromosome=""
   gene=""
   i=0
   bool_first=True
   ctg_out=False
   for line in map_file:
      r_AGP=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r_AGP:
         chrom=r_AGP.group(1)
         start=r_AGP.group(2)
         end=r_AGP.group(3)
         nb=r_AGP.group(4)
         status=r_AGP.group(5)
         ctg=r_AGP.group(6)
         start_CTG=r_AGP.group(7)
         end_CTG=r_AGP.group(8)
         orientation=r_AGP.group(9)

         # Allow to consider CTG with multiple locations (where end is annotated with "_chrX" or "_X" where X is an integer)
         # short_ctg=ctg.split("_")[0]+" ("+orientation+")"
         short_ctg=ctg+" ("+orientation+")"

         if status!="N":
            # Contig that are considered by DeCoSTAR
            if ctg in dict_CTG:
               if chromosome!=chrom:
                  if not bool_first:
                     i+=1
                     if i>color_size-1:
                        i=0
                  bool_first=False
                  graph.add_node(chrom,shape="box",color=color_list[i])
                  if(orientation=="+"):
                     graph.add_node(dict_CTG[ctg].g1,color=color_list[i])
                     graph.add_node(dict_CTG[ctg].g2,color=color_list[i])
                     if ctg_out:
                        graph.add_edge(chrom,dict_CTG[ctg].g1,color="grey")
                        ctg_out=False
                     else:
                        graph.add_edge(chrom,dict_CTG[ctg].g1,color="black")
                     if LINKS:
                        if ctg_out:
                           graph.add_edge(dict_CTG[ctg].g1,dict_CTG[ctg].g2,color="grey",label=short_ctg)
                           ctg_out=False
                        else:
                           graph.add_edge(dict_CTG[ctg].g1,dict_CTG[ctg].g2,color="black",label=short_ctg)
                     gene=dict_CTG[ctg].g2
                  elif(orientation=="-"):
                     graph.add_node(dict_CTG[ctg].g1,color=color_list[i])
                     graph.add_node(dict_CTG[ctg].g2,color=color_list[i])
                     if ctg_out:
                        graph.add_edge(chrom,dict_CTG[ctg].g2,color="grey")
                        ctg_out=False
                     else:
                        graph.add_edge(chrom,dict_CTG[ctg].g2,color="black")
                     if LINKS:
                        if ctg_out:
                           graph.add_edge(dict_CTG[ctg].g2,dict_CTG[ctg].g1,color="grey",label=short_ctg)
                           ctg_out=False
                        else:
                           graph.add_edge(dict_CTG[ctg].g2,dict_CTG[ctg].g1,color="black",label=short_ctg)
                     gene=dict_CTG[ctg].g1
                  elif(orientation=="?" or orientation=="0" or orientation=="na"):
                     graph.add_node(short_ctg,color=color_list[i])
                     if ctg_out:
                        graph.add_edge(chrom,short_ctg,color="grey")
                        ctg_out=False
                     else:
                        graph.add_edge(chrom,short_ctg ,color="black")
                     gene=short_ctg
                  chromosome=chrom
               else:
                  
                  if(orientation=="+"):
                     graph.add_node(dict_CTG[ctg].g1,color=color_list[i])
                     graph.add_node(dict_CTG[ctg].g2,color=color_list[i])
                     if ctg_out:
                        graph.add_edge(gene,dict_CTG[ctg].g1,color="grey")
                        ctg_out=False
                     else:
                        graph.add_edge(gene,dict_CTG[ctg].g1,color="black")
                     if LINKS:
                        graph.add_edge(dict_CTG[ctg].g1,dict_CTG[ctg].g2,color="black",label=short_ctg)
                     gene=dict_CTG[ctg].g2
                  elif(orientation=="-"):
                     graph.add_node(dict_CTG[ctg].g1,color=color_list[i])
                     graph.add_node(dict_CTG[ctg].g2,color=color_list[i])
                     if ctg_out:
                        graph.add_edge(gene,dict_CTG[ctg].g2,color="grey")
                        ctg_out=False
                     else:
                        graph.add_edge(gene,dict_CTG[ctg].g2,color="black")
                     if LINKS:
                        graph.add_edge(dict_CTG[ctg].g2,dict_CTG[ctg].g1,color="black",label=short_ctg)
                     gene=dict_CTG[ctg].g1
                  elif(orientation=="?" or orientation=="0" or orientation=="na"):
                     graph.add_node(short_ctg,color=color_list[i])
                     if ctg_out:
                        graph.add_edge(gene,short_ctg,color="grey")
                        ctg_out=False
                     else:
                        graph.add_edge(gene,short_ctg ,color="black")
                     gene=short_ctg

            else: #If CTG is not considered by DeCoSTAR print it in "grey"
               ctg_out=True
               if chromosome!=chrom:
                  if not bool_first:
                     i+=1
                     if i>color_size-1:
                        i=0
                  bool_first=False
                  graph.add_node(chrom,shape="box",color=color_list[i])
                  graph.add_node(short_ctg,color="grey")
                  graph.add_edge(chrom,short_ctg,color="grey")
                  gene=short_ctg
                  chromosome=chrom
               else:
                  graph.add_node(short_ctg,color="grey")
                  graph.add_edge(gene,short_ctg,color="grey")
                  gene=short_ctg

      else:
         exit("!!! ERROR, line:\n\t"+line+" of file "+NewADJ_file+" is incorrectly written !!!")
   map_file.close()

   graph.write(OUTPUT_dot)

   command_line="dot -Tsvg "+OUTPUT_dot+" -o "+OUTPUT_svg
   subprocess.call(command_line,shell=True)

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
