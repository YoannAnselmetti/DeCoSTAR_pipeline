#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:                                                                 
###       Create scaffolding adj file for all species from BESST score      
###       files                                                             
###                                                                         
###   INPUT:                                                                
###      1- BESST directory                                                 
###         (data/DATA_SEQ/SCAFFOLDING/BESST-2.2.6/Bowtie2_ALL/TRIMMOMATIC3/ALL)
###      2- Maximum distance between linked contigs                         
###         (Ex: 1000000000)                                                
###      3- Links number threshold in scaffolding adjacencies               
###         (Ex: 3)                                                         
###      4- OUTPUT file                                                     
###         (data/data_DeCoSTAR/scaff_BESST_ALL_3_TRIMMOMATIC3)
###                                                                         
###   OUTPUT:                                                               
###      - scaffolding adj file for all species to create instance for      
###        ART-DeCo_seq. 1st step, OUTPUT file need to be refined by:       
###           => 05b-scaff_adj_file_final.py                                
###                                                                         
###   Name: create_scaff_adj_prefile.py   Author: Yoann Anselmetti      
###   Creation date: 2015/12/02           Last modification: 2017/11/09
###


from sys import argv
from re import search
from os import close, listdir, path, makedirs
from datetime import datetime
from collections import namedtuple   #New in version 2.6
import errno


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise


def store_ADJ(BESST_dir,species,dict_spe_edge_scaff):
   for SRX in sorted(listdir(BESST_dir+"/"+species)):
      if SRX!="BESST_output" and SRX!="newADJ_BESST_complete":
         for SRR in sorted(listdir(BESST_dir+"/"+species+"/"+SRX)):
            ctg_scaff_graph_file=""
            # If BAM merged
            if SRR=="BESST_output":
               ctg_scaff_graph_file=BESST_dir+"/"+species+"/"+SRX+"/BESST_output/score_file_pass_1.tsv"
            # If BAM NOT merged
            else:
               ctg_scaff_graph_file=BESST_dir+"/"+species+"/"+SRX+"/"+SRR+"/BESST_output/score_file_pass_1.tsv"

            # print ctg_scaff_graph_file
            scaff_graph=open(ctg_scaff_graph_file,'r')
            for line in scaff_graph:
               r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\n\t]*)\n$',line)
               if r:
                  ctg1=r.group(1)
                  ori1=r.group(2)
                  ctg2=r.group(3)
                  ori2=r.group(4)
                  gap=r.group(5)
                  var_score=r.group(6)
                  disp_score=r.group(7)
                  links_nb=r.group(8)

                  if ctg1!="scf1/ctg1":
                     if int(links_nb)>links_min and float(gap)<=gap_max:
                        score1='{0:.12f}'.format(float(var_score))
                        score2='{0:.12f}'.format(float(disp_score))

                        adj=ADJ(species,ctg1,ctg2,ori1,ori2)
                        edge=EDGE(species,ctg1,ctg2,ori1,ori2,float(gap),score1,score2,int(links_nb))

                        ##########
                        ### WARNING!!! => SUPPOSE THAT BESST ALWAYS PRINT CTG PAIR IN THE SAME ORDER
                        ##########
                        # If the current edge is already in dict_spe_edge_scaff.
                        if adj in dict_spe_edge_scaff[species]:
                           if dict_spe_edge_scaff[species][adj].vscore<score1:
                              dict_spe_edge_scaff[species][adj]=edge
                           elif dict_spe_edge_scaff[species][adj].vscore==score1:
                              # If vscore are equals take the adj with the best dscore
                              if dict_spe_edge_scaff[species][adj].dscore<score2:
                                 dict_spe_edge_scaff[species][adj]=edge
                              elif dict_spe_edge_scaff[species][adj].dscore==score2:
                                 # If dscore are equals take the adaj with the higher number of links
                                 if dict_spe_edge_scaff[species][adj].link<int(links_nb):
                                    dict_spe_edge_scaff[species][adj]=edge
                        else:
                           dict_spe_edge_scaff[species][adj]=edge
               else:
                  exit("ERROR, the line:\n\t"+line+"\nis incorrectly written in file "+ctg_scaff_graph_file)
            scaff_graph.close()
   return dict_spe_edge_scaff


################
###   MAIN 
################
if __name__ == '__main__':

   start_time = datetime.now()

   BESST_dir="data/DATA_SEQ/SCAFFOLDING/BESST-2.2.6/Bowtie2_ALL/TRIMMOMATIC3/ALL"
   gap_max=1000000000000
   links_min=3
   OUTPUT_file="data/data_DeCoSTAR/scaff_BESST_ALL_3_TRIMMOMATIC3"

   OUTPUT_DIR=path.dirname(OUTPUT_file)

   # To be sure than directory have no "/" to the end of the path
   OUTPUT_DIR=path.normpath(OUTPUT_DIR)

   # Create OUTPUT_DIR if not existing
   if not path.exists(OUTPUT_DIR):
      mkdir_p(OUTPUT_DIR)


   # STRUCTURE for gene edge and adjacency (EDGE==ADJ)
   EDGE=namedtuple("EDGE",["spe","ctg1","ctg2","ori1","ori2","gap","vscore","dscore","link"])
   ADJ=namedtuple("ADJ",["spe","ctg1","ctg2","ori1","ori2"])


###################################################################################################################
### BROWSE BESST DIRECTORY TO GET SCAFFOLDING ADJACENCIES PROPOSED BY BESST AND STORE IT IN dict_spe_edge_scaff ###
###################################################################################################################
   dict_spe_edge_scaff={}
   for species in sorted(listdir(BESST_dir)):
      print species
      dict_spe_edge_scaff[species]=dict()
      dict_spe_edge_scaff = store_ADJ(BESST_dir,species,dict_spe_edge_scaff)


#######################################################
### WRITE OUTPUT_FILE => PRE SCAFFOLDING GRAPH FILE ###
#######################################################
   output=open(OUTPUT_file,"w")
   output.write("#species\tctg1\tctg2\torientation_ctg1\torientation_ctg2\tctg1-ctg2_dist\tvscore\tdscore\t#links\n")
   for spe in sorted(dict_spe_edge_scaff):
      # print spe
      for adj in sorted(dict_spe_edge_scaff[spe]):
         # print "\t",
         edge=dict_spe_edge_scaff[spe][adj]
         output.write(edge.spe+"\t"+edge.ctg1+"\t"+edge.ctg2+"\t"+edge.ori1+"\t"+edge.ori2+"\t"+str(edge.gap)+"\t"+str(edge.vscore)+"\t"+str(edge.dscore)+"\t"+str(edge.link)+"\n")
   output.close()


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
