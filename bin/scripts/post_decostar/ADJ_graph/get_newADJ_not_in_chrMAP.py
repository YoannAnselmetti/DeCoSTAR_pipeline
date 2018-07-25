#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create contigs file with gene located to the extremities of contigs
###
###   INPUT:
###      1- File containing scaffolds assignment to chromosome (AGP file)
###         (9passeriformes_dataset/data/INPUT_DATA/ADJ_graph/Zbor_scaff_anchored_in_Tgut_chr_050618_filt.agp)
###      2- CTG file
###         (9passeriformes_dataset/data/data_DeCoSTAR/CTG_file)
###      3- Species corresponding to the INPUT DOT file contig graph
###         (Zosterops_borbonicus)
###      4- New EXTANT adjacencies proposed by DeCo*
###         (9passeriformes_dataset/results/decostar/ADseq+scaff_Boltz_kT0.1/DeCoSTAR_9passeriformes_ADseq+scaff_Boltz_kT0.1_Lin0.1_M1_new_extant_adjacencies_with_scaff_adj)
###      5- OUTPUT file containing the list of new adjacencies predicted by DeCOSTAR not present in chromosome map (AGP file) 
###         (RESULTS/7sauria/ADJ_graph/SVG/CTG/GM_unfiltered/GM_allchr_CTG.svg)
###
###   OUTPUT:
###      - File with new adjacencies predicted by DeCoSTAR not present in chromosome map 
###
###   Name: get_newADJ_not_in_chrMAP.py           Author: Yoann Anselmetti
###   Creation date: 2016/06/13                   Last modification: 2018/06/28
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


def find_newadj_path(gene,dict_newadj,list_gene,dict_geneCTG,list_out):
   # Analyze adjacenies of gene
   if gene not in list_gene:
      list_newadj=dict_newadj[gene]
      for newadj in list_newadj:
         if not newadj in list_out:
            list_out.append(newadj)
            if gene==newadj.g1:
               find_newadj_path(newadj.g2,dict_newadj,list_gene,dict_geneCTG,list_out)
            elif gene==newadj.g2:
               find_newadj_path(newadj.g1,dict_newadj,list_gene,dict_geneCTG,list_out)

      # Get gene localized to the other CTG end and get new adj of this gene
      gene2=""
      ctg=dict_geneCTG[gene]
      if gene==dict_CTG[ctg].g1:
         gene2=dict_CTG[ctg].g2
      elif gene==dict_CTG[ctg].g2:
         gene2=dict_CTG[ctg].g1

      if gene2 not in list_gene and gene2!=gene:
         list_newadj2=dict_newadj[gene]
         for newadj2 in list_newadj2:
            if not newadj2 in list_out:
               list_out.append(newadj2)
               if gene2==newadj2.g1:
                  find_newadj_path(newadj2.g2,dict_newadj,list_gene,dict_geneCTG,list_out)
               elif gene2==newadj2.g2:
                  find_newadj_path(newadj2.g1,dict_newadj,list_gene,dict_geneCTG,list_out)



################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   AGP_file=argv[1]
   CTG_file=argv[2]
   INPUT_species=argv[3]
   NewADJ_file=argv[4]
   OUTPUT_file=argv[5]

   OUTPUT_DIR=path.dirname(path.realpath(OUTPUT_file))
   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)

   CTG=namedtuple("CTG",["g1","g2"])
   ADJ=namedtuple("ADJ",["g1","g2"])

   NEWADJ=namedtuple("ADJ",["ctg1","ctg2","oriC1","oriC2","distC","gf1","gf2","g1","g2","oriG1","oriG2","distG","vscore","dscore","links","supp"])


#######
### GET GENE ANNOTATION AND STORE IT IN dict_gene
#######
   dict_CTG=dict()
   dict_geneCTG=dict()
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
            dict_geneCTG[g1]=ctg
            dict_geneCTG[g2]=ctg



#######
### BROWSE CHR MAP TO WRITE CTG ADJACENCIES GRAPH
#######
   dict_gmADJ=dict()
   # Browse new extant adjacencies file and store adj of selected species in "list_ADJ"
   map_file=open(AGP_file,'r')
   chromosome=""
   gene=""
   i=0
   bool_first=True
   listCTG_chrMAP=list()
   list_CTH_chrMAP_noGENE=list()
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

         if status!="N":
            if ctg in dict_CTG:
               listCTG_chrMAP.append(ctg)
               if chromosome!=chrom:
                  if(orientation=="+"):
                     gene=dict_CTG[ctg].g2
                  elif(orientation=="-"):
                     gene=dict_CTG[ctg].g1
                  chromosome=chrom
                  g1=dict_CTG[ctg].g1
                  g2=dict_CTG[ctg].g2
                  if not g1 in dict_gmADJ:
                     dict_gmADJ[g1]=list()
                  dict_gmADJ[g1].append(ADJ(g1,g2))
                  if not g2 in dict_gmADJ:
                     dict_gmADJ[g2]=list()
                  dict_gmADJ[g2].append(ADJ(g1,g2))
               else:
                  g1=gene
                  g2,g3="",""
                  if(orientation=="+"):
                     g2=dict_CTG[ctg].g1
                     g3=dict_CTG[ctg].g2
                     gene=dict_CTG[ctg].g2
                  elif(orientation=="-"):
                     g2=dict_CTG[ctg].g2
                     g3=dict_CTG[ctg].g1
                     gene=dict_CTG[ctg].g1

                  if not g1 in dict_gmADJ:
                     dict_gmADJ[g1]=list()
                  dict_gmADJ[g1].append(ADJ(g1,g2))
                  if not g2 in dict_gmADJ:
                     dict_gmADJ[g2]=list()
                  dict_gmADJ[g2].append(ADJ(g1,g2))
                  # dict_gmADJ[g2].append(ADJ(g2,g3))
                  if not g3 in dict_gmADJ:
                     dict_gmADJ[g3]=list()
                  dict_gmADJ[g3].append(ADJ(g2,g3))

            else:
               list_CTH_chrMAP_noGENE.append(ctg)
               print "!!!WARNING!!! => The scaffold "+ctg+" is not present in CTG_file "+CTG_file
      else:
         exit("!!! ERROR, line:\n\t"+line+" of file "+NewADJ_file+" is incorrectly written !!!")
   map_file.close()


##################################################################
### BROWSE NEW EXTANT ADJACENCIES FILE PROPOSED BY ARt-DeClone ###
##################################################################
   dict_newADJ=dict()
   list_newCTG=list()
   # Browse new extant adjacencies file and store adj of selected species in "list_ADJ"
   newadj_file=open(NewADJ_file,'r')
   for line in newadj_file:
      species=""
      ctg1=""
      ctg2=""
      oriC1=""
      oriC2=""
      distC=""
      gf1=""
      gf2=""
      g1=""
      g2=""
      oriG1=""
      oriG2=""
      distG=""
      vscore="?"
      dscore="?"
      links="?"
      support="?"

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
         if ctg1 in listCTG_chrMAP and ctg2 not in listCTG_chrMAP:
            list_newCTG.append(ctg2)
         if ctg2 in listCTG_chrMAP and ctg1 not in listCTG_chrMAP:
            list_newCTG.append(ctg1)
         if not g1 in dict_newADJ:
            dict_newADJ[g1]=list()
         dict_newADJ[g1].append(NEWADJ(ctg1,ctg2,oriC1,oriC2,distC,gf1,gf2,g1,g2,oriG1,oriG2,distG,vscore,dscore,links,support))
         if not g2 in dict_newADJ:
            dict_newADJ[g2]=list()
         dict_newADJ[g2].append(NEWADJ(ctg1,ctg2,oriC1,oriC2,distC,gf1,gf2,g1,g2,oriG1,oriG2,distG,vscore,dscore,links,support))

   newadj_file.close()


   #######################################################################################################################
   ### Browse gene list of GM and determine if there is new adajcency proposed for this gene that is not present in GM ###
   #######################################################################################################################
   list_outADJ=list()
   for gene in dict_gmADJ:
      # print gene
      if gene in dict_newADJ:
         # print "\t=> "+gene
         for newADJ in dict_newADJ[gene]:
            # print "\t\tnewADJ => ",
            # print newADJ
            adj_exists=False
            for gmADJ in dict_gmADJ[gene]:
               # print "\t\tgmADJ => ",
               # print gmADJ
               if (gmADJ.g1==newADJ.g1 and gmADJ.g1==newADJ.g1) or (gmADJ.g1==newADJ.g2 and gmADJ.g1==newADJ.g2):
                  adj_exists=True
                  break
            if not adj_exists:
               if not newADJ in list_outADJ:
                  # print "!!!=> ",
                  # print newADJ
                  list_outADJ.append(newADJ)
                  if gene==newADJ.g1:
                     find_newadj_path(newADJ.g2,dict_newADJ,dict_gmADJ.keys(),dict_geneCTG,list_outADJ)
                  elif gene==newADJ.g2:
                     find_newadj_path(newADJ.g1,dict_newADJ,dict_gmADJ.keys(),dict_geneCTG,list_outADJ)

   #####################################################################################
   ### Print new adjacencies connected to Genetic map but not present in Genetic Map ###
   #####################################################################################
   output=open(OUTPUT_file,"w")
   output.write("ctg1\tctg2\torientation_ctg1\torientation_ctg2\tdist_ctg1-ctg2\tgene_family1\tgene_family2\tgene1\tgene2\torientation_gene1\torientation_gene2\tdist_gene1-gene2\tlinks_variation_score\tlinks_dispersity_score\t#links\tDeClone_support\n")
   for adj in list_outADJ:
      # print adj.ctg1+"\t"+adj.ctg2+"\t"+adj.oriC1+"\t"+adj.oriC2+"\t"+adj.distC+"\t"+adj.gf1+"\t"+adj.gf2+"\t"+adj.g1+"\t"+adj.g2+"\t"+adj.oriG1+"\t"+adj.oriG2+"\t"+adj.distG+"\t"+adj.vscore+"\t"+adj.dscore+"\t"+adj.links+"\t"+adj.supp
      output.write(adj.ctg1+"\t"+adj.ctg2+"\t"+adj.oriC1+"\t"+adj.oriC2+"\t"+adj.distC+"\t"+adj.gf1+"\t"+adj.gf2+"\t"+adj.g1+"\t"+adj.g2+"\t"+adj.oriG1+"\t"+adj.oriG2+"\t"+adj.distG+"\t"+adj.vscore+"\t"+adj.dscore+"\t"+adj.links+"\t"+adj.supp+"\n")
   output.close()

   print "List of CTG that are linked to the chromosome MAP by DeCoSTAR:"
   for ctg in list_newCTG:
      print ctg

   print "\nList of CTG WITH gene in INPUT chromosome MAP:"
   for ctg in listCTG_chrMAP:
      print ctg

   print "\nList of CTG WITHOUT gene in INPUT chromosome MAP:"
   for ctg in list_CTH_chrMAP_noGENE:
      print ctg


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))