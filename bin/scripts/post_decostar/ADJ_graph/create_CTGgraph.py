#! /usr/bin/env python3
# -*- coding: utf-8 -*-
###
###   Goal:
###      Add adjacencies predicted by DeCoSTAR on chromosome graph where node corresponds to contig/scaffold ID
###
###   INPUT:
###      1- INPUT dot file for contig graph (with Genetic Map)
###         (27avian_dataset/results/ADJ_graph/DOT/CTG/Zosterops_borbonicus/Zosterops_borbonicus_CTGmap.dot)
###      2- New EXTANT adjacencies proposed by DeCo*
###         (27avian_dataset/results/decostar/ADseq+scaff_Boltz_kT0.1/DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Lin0.1_M2_new_extant_adjacencies_with_scaff)
###      3- Species corresponding to the INPUT DOT file contig graph
###         (Zosterops_borbonicus)
###      4- OUTPUT file where contig graph will be stored
###         (27avian_dataset/results/ADJ_graph/DOT/CTG/Zosterops_borbonicus/Zosterops_borbonicus_DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Lin0.1_M2_CTGmap.dot)
###      5- OUTPUT SVG file where Contig graph will be stored
###         (27avian_dataset/results/ADJ_graph/SVG/CTG/Zosterops_borbonicus/Zosterops_borbonicus_DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Lin0.1_M2_CTGmap.svg)
###      6- Minimum support
###         (Ex: 0.0)
###      7- Keep scaffolding adjacencies (same if DeClone support < Min support)
###         (y/Y: YES | n/N: NO)
###
###   OUTPUT:
###      - Adjacencies graph files in SVG and DOT format with node corresponding to contig/scaffold ID
###      - Combination of chromosome map with DeCoSTAR prediction (could be improved to take several methods: GOS-ASM for example)
###
###   Name: create_CTGgraph.py                    Author: Yoann Anselmetti
###   Creation date: 2016/03/07                   Last modification: 2020/11/09
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
   INPUT_DOT_file=argv[1]
   NewADJ_file=argv[2]
   INPUT_species=argv[3]
   OUTPUT_file=argv[4]
   OUTPUT_svg=argv[5]
   min_support=argv[6]
   bool_links=argv[7]

   LINKS=False
   if bool_links=="n" or bool_links=="N":
      LINKS=False
   elif bool_links=="y" or bool_links=="Y":
      LINKS=True
   else:
      exit("!!! ERROR, 7th parameter should be y/y or N/n and not "+bool_links+"!!!")

   OUTPUT_DIR=path.dirname(path.realpath(OUTPUT_file))
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

   min_support=float(min_support)

##################################################################
### BROWSE NEW EXTANT ADJACENCIES FILE PROPOSED BY ARt-DeClone ###
##################################################################
   list_ADJ=list()
   # Browse new extant adjacencies file and store adj of selected species in "list_ADJ"
   newadj_file=open(NewADJ_file,'r')
   for line in newadj_file:
      ctg1=""
      ctg2=""
      support="?"
      links="?"
      r_ADseqDC=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      r_ADseq=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      r_DC=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r_ADseqDC:
         species=r_ADseqDC.group(1)
         ctg1=r_ADseqDC.group(2)
         ctg2=r_ADseqDC.group(3)
         oriC1=r_ADseqDC.group(4)
         oriC2=r_ADseqDC.group(5)
         distC=r_ADseqDC.group(6)
         gf1=r_ADseqDC.group(7)
         gf2=r_ADseqDC.group(8)
         g1=r_ADseqDC.group(9)
         g2=r_ADseqDC.group(10)
         oriG1=r_ADseqDC.group(11)
         oriG2=r_ADseqDC.group(12)
         distG=r_ADseqDC.group(13)
         vscore=r_ADseqDC.group(14)
         dscore=r_ADseqDC.group(15)
         links=r_ADseqDC.group(16)
         support=r_ADseqDC.group(17)
      elif r_ADseq:
         species=r_ADseq.group(1)
         ctg1=r_ADseq.group(2)
         ctg2=r_ADseq.group(3)
         oriC1=r_ADseq.group(4)
         oriC2=r_ADseq.group(5)
         distC=r_ADseq.group(6)
         gf1=r_ADseq.group(7)
         gf2=r_ADseq.group(8)
         g1=r_ADseq.group(9)
         g2=r_ADseq.group(10)
         oriG1=r_ADseq.group(11)
         oriG2=r_ADseq.group(12)
         distG=r_ADseq.group(13)
         vscore=r_ADseq.group(14)
         dscore=r_ADseq.group(15)
         links=r_ADseq.group(16)
      elif r_DC:
         species=r_DC.group(1)
         ctg1=r_DC.group(2)
         ctg2=r_DC.group(3)
         oriC1=r_DC.group(4)
         oriC2=r_DC.group(5)
         distC=r_DC.group(6)
         gf1=r_DC.group(7)
         gf2=r_DC.group(8)
         g1=r_DC.group(9)
         g2=r_DC.group(10)
         oriG1=r_DC.group(11)
         oriG2=r_DC.group(12)
         distG=r_DC.group(13)
         support=r_DC.group(14)
      else:
         exit("!!! ERROR, line:\n\t"+line+" of file "+NewADJ_file+" is incorrectly written !!!")

      if species==INPUT_species:
         adj_info=ADJ(ctg1,ctg2,support,links)
         list_ADJ.append(adj_info)

   newadj_file.close()




#####################################################################
### BROWSE "list_ADJ" LIST TO CREATE CTG GRAPH FILE TO DOT FORMAT ###
#####################################################################
   graph=pgv.AGraph()
   if path.exists(INPUT_DOT_file):
      graph=pgv.AGraph(INPUT_DOT_file)
   else:
      graph.node_attr['penwidth']=5
      graph.edge_attr['penwidth']=5
      graph.node_attr['fontsize']=20
      graph.edge_attr['fontsize']=20
      graph.node_attr['fontcolor']="black"
      graph.edge_attr['fontcolor']="black"
      graph.node_attr['fontname']="times-bold"
      graph.edge_attr['fontname']="times-bold"
   for adj in list_ADJ:
      # ctg1=adj.ctg1.split("_")[0]
      # ctg2=adj.ctg2.split("_")[0]
      ctg1=adj.ctg1
      ctg2=adj.ctg2
      # print(adj)
      # If the ADJ cooresponds to a CTG
      if adj.sup=="?":
         # If ARt-DeCo results (adj without DeClone support)
         if adj.links=="N/A":
            if not ctg1 in graph:
               graph.add_node(ctg1,color="black")
            if not ctg2 in graph:
               graph.add_node(ctg2,color="black")
            graph.add_edge(ctg1,ctg2,color="blue")
         # If ADseq results (scaff adj without DeClone support)
         else:
            if not ctg1 in graph:
               graph.add_node(ctg1,color="black")
            if not ctg2 in graph:
               graph.add_node(ctg2,color="black")
            graph.add_edge(ctg1,ctg2,color="blue",label=adj.links)
      elif (float(adj.sup)>=min_support or (LINKS and adj.links!="?")):
         TAG=""
         if adj.links!="N/A":
            TAG=adj.sup+"\n("+adj.links+")"
         else:
            TAG=adj.sup
         if not ctg1 in graph:
            graph.add_node(ctg1,color="black")
         if not ctg2 in graph:
            graph.add_node(ctg2,color="black")
         if float(adj.sup)<=0.1:
            graph.add_edge(ctg1,ctg2,color="#FF0000",label=TAG)
         elif float(adj.sup)<=0.2:
            graph.add_edge(ctg1,ctg2,color="#E31100",label=TAG)
         elif float(adj.sup)<=0.3:
            graph.add_edge(ctg1,ctg2,color="#C62200",label=TAG)
         elif float(adj.sup)<=0.4:
            graph.add_edge(ctg1,ctg2,color="#AA3300",label=TAG)
         elif float(adj.sup)<=0.5:
            graph.add_edge(ctg1,ctg2,color="#8E4400",label=TAG)
         elif float(adj.sup)<=0.6:
            graph.add_edge(ctg1,ctg2,color="#715500",label=TAG)
         elif float(adj.sup)<=0.7:
            graph.add_edge(ctg1,ctg2,color="#556600",label=TAG)
         elif float(adj.sup)<=0.8:
            graph.add_edge(ctg1,ctg2,color="#397700",label=TAG)
         elif float(adj.sup)<=0.9:
            graph.add_edge(ctg1,ctg2,color="#1C8800",label=TAG)
         elif float(adj.sup)<=1.0:
            graph.add_edge(ctg1,ctg2,color="#009900",label=TAG)

   graph.write(OUTPUT_file)

   command_line="dot -Tsvg "+OUTPUT_file+" -o "+OUTPUT_svg
   subprocess.call(command_line,shell=True)

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time-start_time))
