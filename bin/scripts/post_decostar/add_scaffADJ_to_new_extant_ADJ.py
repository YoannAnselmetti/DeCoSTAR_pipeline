#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create New adjacencies file of DeCoSTAR with CTG infos and SCAFF infos
###
###   INPUT:
###      1- INPUT new adjacencies file proposed by ARt-DeCoSTAR with CTG infos
###         (9passeriformes_dataset/results/decostar/ADseq+scaff_Boltz_kT0.1/DeCoSTAR_9passeriformes_ADseq+scaff_Boltz_kT0.1_Lin0.1_M1_new_extant_adjacencies)
###      2- INPUT file containing scaffolding adj proposed by SCAFF
###         (9passeriformes_dataset/data/data_DeCoSTAR/scaff_BESST_DeCoSTAR)
###      3- OUTPUT file path for New adjacencies file with contig ID
###         (9passeriformes_dataset/results/decostar/ADseq+scaff_Boltz_kT0.1/DeCoSTAR_9passeriformes_ADseq+scaff_Boltz_kT0.1_Lin0.1_M1_new_extant_adjacencies_with_scaff_adj)
###
###   OUTPUT:  (run in few sec)
###      - New extant adjacencies file predicted by DeCoSTAR with CTG/SCAFF infos
###
###   Name: 05-add_scaffADJ_to_new_extant_ADJ.py  Author: Yoann Anselmetti
###   Creation date: 2016/03/07                   Last modification: 2018/06/26
###

from sys import argv, stdout
from re import search
from os import close, path
from datetime import datetime
from collections import namedtuple   #New in version 2.6
import subprocess


def rev_ori(ori):
   if ori=="-":
      return "+"
   elif ori=="+":
      return "-"
   elif ori=="?":
      return "?"
   else:
      exit("ERROR, orientation: \""+ori+"\" is incorrect, it should be \"+\" or \"-\"!!!")



def sameADJ(stored,scaff):
   if (stored.g1==scaff.g1 and stored.g2==scaff.g2):
      if (stored.ori1==scaff.ori1 and stored.ori2==scaff.ori2):
         return True
      else:
         return False
   elif (stored.g1==scaff.g2 and stored.g2==scaff.g1):
      if (stored.ori1==rev_ori(scaff.ori2) and stored.ori2==rev_ori(scaff.ori1)):
         return True
      else:
         return False
   else:
      return False 


################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   NewADJ_file=argv[1]
   SCAFF_file=argv[2]
   OUTPUT_file=argv[3]

   shortADJ=namedtuple("shortADJ",["g1","g2","ori1","ori2"])
   ADJ=namedtuple("ADJ",["spe","ctg1","ctg2","oriC1","oriC2","gf1","gf2","g1","g2","oriG1","oriG2","sup"])

##################
### BROWSE NEW EXTANT ADJACENCIES FILE PROPOSED BY DeCoSTAR
##################
   print "1/ Store new extant adjacencies predicted by DeCoSTAR... ",
   dict_ADJ=dict()
   adj_in_DeCoSTAR=0
   # Browse new extant adjacencies file and store adj of selected species in "list_ADJ"
   newadj_file=open(NewADJ_file,'r')
   for line in newadj_file:
      r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r:
         species=r.group(1)
         ctg1=r.group(2)
         ctg2=r.group(3)
         oriC1=r.group(4)
         oriC2=r.group(5)
         gf1=r.group(6)
         gf2=r.group(7)
         g1=r.group(8)
         g2=r.group(9)
         oriG1=r.group(10)
         oriG2=r.group(11)
         support=r.group(12)

         if species!="#species":
            adj_in_DeCoSTAR+=1
            adj=ADJ(species,ctg1,ctg2,oriC1,oriC2,gf1,gf2,g1,g2,oriG1,oriG2,support)
            if not species in dict_ADJ:
               dict_ADJ[species]=list()
            dict_ADJ[species].append(adj)
   newadj_file.close()
   print "DONE\n"



##################
### BROWSE NEW ADJ FILE PROPOSED BY DeCoSTAR
##################
   print "2/ Write new adjacencies predicted by DeCoSTAR with SCAFF infos in file "+OUTPUT_file+"... ",
   stdout.flush()
   dict_spe_newAdj={}
   scaff_adj_nb=0
   common_adj=0
   list_scaff_not_in_DeCoSTAR=list()
   scaff_file=open(SCAFF_file,'r')
   output_file=open(OUTPUT_file,'w')
   for line in scaff_file:
      r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r:
         species=r.group(1)
         ctg1=r.group(2)
         ctg2=r.group(3)
         oriC1=r.group(4)
         oriC2=r.group(5)
         distC=r.group(6)
         gf1=r.group(7)
         gf2=r.group(8)
         g1=r.group(9)
         g2=r.group(10)
         oriG1=r.group(11)
         oriG2=r.group(12)
         distG=r.group(13)
         vscore=r.group(14)
         dscore=r.group(15)
         links=r.group(16)

         scaffADJ=shortADJ(g1,g2,oriG1,oriG2)

         if species!="#species":
            bool_SCAFF=False
            scaff_adj_nb+=1
            if species in dict_ADJ:
               for adj in dict_ADJ[species]:
                  storedADJ=shortADJ(adj.g1,adj.g2,adj.oriG1,adj.oriG2)
                  if sameADJ(storedADJ,scaffADJ):
                     common_adj+=1
                     output_file.write(species+"\t"+ctg1+"\t"+ctg2+"\t"+oriC1+"\t"+oriC2+"\t"+distC+"\t"+gf1+"\t"+gf2+"\t"+g1+"\t"+g2+"\t"+oriG1+"\t"+oriG2+"\t"+distG+"\t"+vscore+"\t"+dscore+"\t"+links+"\t"+adj.sup+"\n")
                     dict_ADJ[species].remove(adj)
                     bool_SCAFF=True
                     break
            if not bool_SCAFF:
               list_scaff_not_in_DeCoSTAR.append(line)
      else:
         exit("\n!!! ERROR, format of line:\n\t"+line+"is incorrect!!!")
   scaff_file.close()

   DeCoSTAR_not_in_SCAFF=0
   for species in dict_ADJ:
      DeCoSTAR_not_in_SCAFF+=len(dict_ADJ[species])
      for adj in dict_ADJ[species]:
         output_file.write(adj.spe+"\t"+adj.ctg1+"\t"+adj.ctg2+"\t"+adj.oriC1+"\t"+adj.oriC2+"\t?\t"+adj.gf1+"\t"+adj.gf2+"\t"+adj.g1+"\t"+adj.g2+"\t"+adj.oriG1+"\t"+adj.oriG2+"\t?\tN/A\tN/A\tN/A\t"+adj.sup+"\n")
   output_file.close()
   print "DONE\n"

   print "3/ Statictis on common adjacencies between BESST and DeCoSTAR:\n"
   print "\tThere are "+str(common_adj)+" scaffolding ADJ proposed by BESST AND DeCoSTAR"
   print "\tThere are "+str(len(list_scaff_not_in_DeCoSTAR))+"/"+str(scaff_adj_nb)+" scaffolding ADJ proposed by BESST that are not predicted by DeCoSTAR"
   print "\tThere are "+str(DeCoSTAR_not_in_SCAFF)+"/"+str(adj_in_DeCoSTAR)+" ADJ predicted by DeCoSTAR that are not proposed by BESST"

##################
### Sort OUTPUT new adjacencies file
##################
   command_line="sort "+OUTPUT_file+" -o "+OUTPUT_file+"; sed -i 1i'#species\tctg1\tctg2\torientation_ctg1\torientation_ctg2\tdist_ctg1-ctg2\tgene_family1\tgene_family2\tgene1\tgene2\torientation_gene1\torientation_gene2\tdist_gene1-gene2\tvscore\tdscore\tlinks\tsupport\n' "+OUTPUT_file
   subprocess.call(command_line,shell=True)

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
