#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create adjacencies graph file where node corresponds to contig/scaffold ID 
###      form a file containing assignment of part of contigs/scaffolds on chromosome
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
###         (9passeriformes_dataset/results/ADJ_graph/DOT/CTG/Zosterops_borbonicus/Zbor_CTG_filt.dot)
###      5- OUTPUT SVG file where Contig graph will be stored
###         (9passeriformes_dataset/results/ADJ_graph/SVG/CTG/Zosterops_borbonicus/Zbor_CTG_filt.svg)
###
###   OUTPUT:
###      - Directory for SAM files containing 1 Directory/Species
###
###   Name: create_CTGgraph.py                    Author: Yoann Anselmetti
###   Creation date: 2016/03/07                   Last modification: 2018/06/28
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
      else:
         raise


################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   AGP_file=argv[1]
   CTG_file=argv[2]
   INPUT_species=argv[3]
   OUTPUT_dot=argv[4]
   OUTPUT_svg=argv[5]


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

   ADJ=namedtuple("ADJ",["ctg1","ctg2","sup","links"])

   color_list=("blue","red","green","purple","gold","orangered","chocolate","blueviolet","indigo","fuchsia","crimson","darkgreen","chartreuse4","darkblue","salmon","darkturquoise","cyan","tomato","teal","orchid","maroon","lime","deeppink","darkorange","dodgerblue","coral","sienna","olive","goldenrod","darkcyan","mediumvioletred","darkred","darkgrey","brown","darkgoldenrod","darkolivegreen","darkseagreen","darkslategrey","darkviolet","darkslateblue","orange","yellow","olivedrab")
   color_size=len(color_list)

#####################################################
### GET GENE ANNOTATION AND STORE IT IN dict_gene ###
#####################################################
   list_CTG=list()
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
            # print line
            list_CTG.append(ctg)


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
   map_file=open(AGP_file,'r')
   chromosome=""
   contig=""
   i=0
   bool_first=True
   list_CTG_rm=list()
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
         ctg_short=ctg.split("_")[0]

         if status!="N":
            if ctg in list_CTG:
               if chromosome!=chrom:
                  if not bool_first:
                     i+=1
                     if i>color_size-1:
                        i=0
                  bool_first=False
                  graph.add_node(chrom,shape="box",color=color_list[i])
                  graph.add_node(ctg_short,color=color_list[i])
                  graph.add_edge(chrom,ctg_short,color="black")
                  contig=ctg_short
                  chromosome=chrom
               else:
                  graph.add_node(ctg_short,color=color_list[i])
                  graph.add_edge(contig,ctg_short,color="black")
                  contig=ctg_short
            else:
               list_CTG_rm.append(ctg)
      else:
         exit("!!! ERROR, line:\n\t"+line+" of file "+AGP_file+" is incorrectly written !!!")
   map_file.close()

   graph.write(OUTPUT_dot)

   command_line="dot -Tsvg "+OUTPUT_dot+" -o "+OUTPUT_svg
   subprocess.call(command_line,shell=True)

   print str(len(list_CTG_rm))+" CTG in AGP file don't contain genes considered by DeCoSTAR:"
   for ctg in list_CTG_rm:
      print "\t"+ctg



   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
