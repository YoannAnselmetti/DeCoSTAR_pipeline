#! /usr/bin/env python2
# -*- coding: utf-8 -*-
###
###   Goal:
###       Create scaffolding adj file for all species from BESST score files
###
###   INPUT:
###      1- BESST directory
###         (data/DATA_SEQ/SCAFFOLDING/BESST-2.2.6/Bowtie2_k100/TRIM/)
###      2- CTG file
###         (data/data_DeCoSTAR/CTG_file)
###      3- Maximum distance between linked contigs
###         (Ex: 1000000000)
###      4- Links number threshold in scaffolding adjacencies
###         (Ex: 4)
###      5- OUTPUT file
###         (data/data_DeCoSTAR/scaff_BESST_ALL_4_TRIM)
###
###   OUTPUT:
###      - Scaffolding adjacencies gene file for all species
###
###   Name: create_scaff_adj_file.py   Author: Yoann Anselmetti
###   Creation date: 2015/12/02               Last modification: 2018/07/26
###

from sys import argv, stdout
from re import search
from os import close, listdir, path, makedirs
from datetime import datetime
from collections import namedtuple   #New in version 2.6
from glob import glob
import errno



def getDIR(file_path):
   return file_path.rsplit("/",1)[0]



def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else:
         raise



def rev_ori(ori):
   if ori=="-":
      return "+"
   elif ori=="+":
      return "-"
   elif ori=="?":
      return "?"
   else:
      exit("ERROR, orientation: \""+ori+"\" is incorrect, it should be \"+\" or \"-\"!!!")



def store_CTG(CTG_file):
   CTG_format="#species\tctg\tctg_size\tctg_gene_nb\t5'_gene_family\t5'_gene\torientation_5'_gene\tstart_5'_gene\t3'_gene_family\t3'_gene\torientation_3'_gene\tend_3'_gene\n"
   dict_spe_ctg=dict()
   # When get a contigs pairs (edge scaffolding link) => Get genes that are linked by scaffolding graph with the distance
   contig_file=open(CTG_file,'r')
   for line in contig_file:
      r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r:
         spe=r.group(1)
         contig=r.group(2)
         contig_size=r.group(3)
         contig_geneNb=r.group(4)
         GF1=r.group(5)
         g1=r.group(6)
         oriG1=r.group(7)
         start_g1=r.group(8)
         GF2=r.group(9)
         g2=r.group(10)
         oriG2=r.group(11)
         stop_g2=r.group(12)

         if spe!="#species":
            if contig_size!="?": # Cause contig with "size==?"" are not present in genome assemblies
               ctg=CTG(spe,contig,int(contig_size),GF1,g1,oriG1,int(start_g1),GF2,g2,oriG2,int(stop_g2))
               if not spe in dict_spe_ctg:
                  dict_spe_ctg[spe]=dict()
               dict_spe_ctg[spe][contig]=ctg
            # else:
            #    print "\t=> Contig "+contig+" is not present in FASTA file assembly of species "+spe
      else:
         exit("ERROR, line "+line+" of file "+CTG_file+" is incorrectly written!!!\nIt should match with the following format:\n"+CTG_format)
   contig_file.close()

   return dict_spe_ctg



def best_adj(stored_edge,current_edge):
   current_score=(current_edge.vscore+current_edge.dscore)/2.0
   stored_score=(stored_edge.vscore+stored_edge.dscore)/2.0
   # If current score > stored score, replace stored edge by current edge
   if current_score>stored_score:
      return True
   # If current score == stored score, keep the edge with the higher number of links 
   elif current_score==stored_score:
      if current_edge.link>stored_edge.link:
         return True
      else:
         return False



def read_and_store_scaff_ADj(ctg_scaff_graph_file,dict_spe_edge_scaff):
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

            # print line

            if ctg1!="scf1/ctg1":
               # Parse ctg1, ori1, ctg2 and ori2 if they are composed of several contigs
               if ";" in ctg1:
                  ctg1=ctg1.rsplit(";",1)[1]
                  ori1=ori1.rsplit(";",1)[1]
               ctg2=ctg2.split(";",1)[0]
               ori2=ori2.split(";",1)[0]

               if int(links_nb)>links_min and float(gap)<=gap_max:
                  score1=float('{0:.12f}'.format(float(var_score)))
                  score2=float('{0:.12f}'.format(float(disp_score)))

                  adj=ADJ(species,ctg1,ctg2,ori1,ori2)
                  edge=EDGE(species,ctg1,ctg2,ori1,ori2,float(gap),score1,score2,int(links_nb))

                  rev_adj=ADJ(species,ctg2,ctg1,rev_ori(ori2),rev_ori(ori1))
                  rev_edge=EDGE(species,ctg2,ctg1,rev_ori(ori2),rev_ori(ori1),float(gap),score1,score2,int(links_nb))

                  # If the current adj is in "dict_spe_edge_scaff"
                  if adj in dict_spe_edge_scaff[species]:
                     if best_adj(dict_spe_edge_scaff[species][adj],edge):
                        dict_spe_edge_scaff[species][adj]=edge
                  else:
                     # If the current adj is in "dict_spe_edge_scaff" in the reverse order
                     if rev_adj in dict_spe_edge_scaff[species]:
                        if verbose:
                           print "\tAdjacency:\n",edge,"\nis present in forward and reverse orientation"
                        if best_adj(dict_spe_edge_scaff[species][rev_adj],rev_edge):
                           dict_spe_edge_scaff[species][rev_adj]=rev_edge
                     # If current adj is not present in "dict_spe_edge_scaff", store it in "dict_spe_edge_scaff"
                     else:
                        dict_spe_edge_scaff[species][adj]=edge
         else:
            exit("ERROR, the line:\n\t"+line+"\nis incorrectly written in file "+ctg_scaff_graph_file)
      scaff_graph.close()



def store_ADJ(BESST_dir,species,dict_spe_edge_scaff):
   score_files_nb=len(glob(BESST_dir+"/"+species+"/BESST_output/score_file_pass_*.tsv"))
   if score_files_nb:
      i=1
      while i<=score_files_nb:
         ctg_scaff_graph_file=BESST_dir+"/"+species+"/BESST_output/score_file_pass_"+str(i)+".tsv"
         print "\t"+ctg_scaff_graph_file
         read_and_store_scaff_ADj(ctg_scaff_graph_file,dict_spe_edge_scaff)
         i+=1
   else:
      for SRX in listdir(BESST_dir+"/"+species):
         ctg_scaff_graph_file=BESST_dir+"/"+species+"/"+SRX+"/BESST_output/score_file_pass_1.tsv"
         print "\t"+ctg_scaff_graph_file
         read_and_store_scaff_ADj(ctg_scaff_graph_file,dict_spe_edge_scaff)

   return dict_spe_edge_scaff



def get_gene_infos(ctg_order,ori,dist,CTG,CTG_file):
   GF_gene=""
   gene=""
   ori_gene=""
   # If gene involved in scaffolding adjacency is in 5' position on contig CTG
   if ((ori=="+" and ctg_order=="second") or (ori=="-" and ctg_order=="first")):
      GF_gene=CTG.gf1
      gene=CTG.g1
      ori_gene=CTG.ori1
      dist+=CTG.start
   # If gene involved in scaffolding adjacency is in 3' position on contig CTG
   elif ((ori=="-" and ctg_order=="second") or (ori=="+" and ctg_order=="first")):
      GF_gene=CTG.gf2
      gene=CTG.g2
      ori_gene=CTG.ori2
      dist+=CTG.size-CTG.end
   else:
      exit("ERROR on gene orientation in contigs extremities file \""+CTG_file+"\" (Should be \"+\" or \"-\")")

   # Change gene orientation if CTG orientation is "-" (ori=="-")
   if (ori=="-"):
      oriG=ori_gene
      ori_gene=rev_ori(oriG)

   return GF_gene,gene,ori_gene,dist





###########
### MAIN 
###########
if __name__ == '__main__':

   start_time = datetime.now()

   BESST_dir=argv[1]
   CTG_file=argv[2]
   gap_max=int(argv[3])
   links_min=int(argv[4])
   OUTPUT_file=argv[5]

   verbose=False

   # Create OUTPOUT_DIR for OUTPUT_file
   OUTPUT_DIR=getDIR(OUTPUT_file)
   mkdir_p(OUTPUT_DIR)

   # STRUCTURE for gene edge and adjacency (EDGE==ADJ)
   EDGE=namedtuple("EDGE",["spe","ctg1","ctg2","ori1","ori2","gap","vscore","dscore","link"])
   ADJ=namedtuple("ADJ",["spe","ctg1","ctg2","ori1","ori2"])
   CTG=namedtuple("CTG",["spe","size"])
   CTG=namedtuple("CTG",["spe","ctg","size","gf1","g1","ori1","start","gf2","g2","ori2","end"])



###########
### STORE CTG INFOS OF "CTG_file" IN "dict_spe_ctg"
###########
   print "1/ Store infos contained in CTG file \""+CTG_file+"\"...",
   stdout.flush()
   dict_spe_ctg=store_CTG(CTG_file)
   print "DONE"

   # for spe in dict_spe_ctg:
   #    for ctg in dict_spe_ctg[spe]:
   #       print ctg+":",dict_spe_ctg[spe][ctg]



###########
### BROWSE "BESST_dir" TO GET SCAFFOLDING ADJACENCIES PROPOSED BY BESST AND STORE IT IN "dict_edge_scaff"
###########
   print "2/ Get, filter and store scaffolding adjacencies predicted by BESST present in directory \""+BESST_dir+"\":"
   Nb_scaff_edge_tot=0
   Nb_scaff_edge_kept=0
   dict_spe_edge_scaff=dict()
   for species in sorted(listdir(BESST_dir)):
      print species+":"
      dict_spe_edge_scaff[species]=dict()
      dict_spe_edge_scaff=store_ADJ(BESST_dir,species,dict_spe_edge_scaff)

   # for spe in dict_spe_edge_scaff:
   #    for adj in dict_spe_edge_scaff[spe]:
   #       print adj,"=>",dict_spe_edge_scaff[spe][adj]



###########
### ADD CTG/GENE INFOS ON SCAFFOLDING ADJACENCIES AND WRITE THEM IN "OUTPUT_file"
###########
   print "3/ Add gene/contig infos to scaffolding adjacencies and write them in OUTPUT file \""+OUTPUT_file+"\"... ",
   stdout.flush()
   scaff_gene_file=open(OUTPUT_file,"w")
   scaff_gene_file.write("#species\tctg1\tctg2\torientation_ctg1\torientation_ctg2\tctg1-ctg2_dist\tgene1_family\tgene2_family\tgene1\tgene2\torientation_gene1\torientation_gene2\tgene1-gene2_dist\tvscore\tdscore\t#links\n")
   for spe in sorted(dict_spe_edge_scaff):
      for adj in sorted(dict_spe_edge_scaff[spe]):
         edge=dict_spe_edge_scaff[spe][adj]
         species=edge.spe
         ctg1=edge.ctg1
         ctg2=edge.ctg2
         oriC1=edge.ori1
         oriC2=edge.ori2
         gap=edge.gap
         vscore=edge.vscore
         dscore=edge.dscore
         link=edge.link

         dist=gap
         Nb_scaff_edge_tot+=1
         bool_ctg1,bool_ctg2=False,False
         GF_gene1,gene1,ori_gene1,GF_gene2,gene2,ori_gene2="","","","","",""
         # If ctg1 is present in CTG_file: Get information on gene involved in the linked between ctg1 and ctg2
         if ctg1 in dict_spe_ctg[species]:
            bool_ctg1=True
            CTG1=dict_spe_ctg[species][ctg1]
            # Get infos for gene1
            gf1,g1,oriG1,dist = get_gene_infos("first",oriC1,dist,CTG1,CTG_file)

         # If ctg2 is present in CTG_file: Get information on gene involved in the linked between ctg1 and ctg2
         if ctg2 in dict_spe_ctg[species]:
            bool_ctg2=True
            CTG2=dict_spe_ctg[species][ctg2]
            # Get infos for gene2
            gf2,g2,oriG2,dist = get_gene_infos("second",oriC2,dist,CTG2,CTG_file)

         # If the 2 contigs are present in file CTG_file: Print scaffolding adj in OUTPUT_scaff_file
         if bool_ctg1 and bool_ctg2:
            scaff_gene_file.write(species+"\t"+ctg1+"\t"+ctg2+"\t"+oriC1+"\t"+oriC2+"\t"+str(gap)+"\t"+gf1+"\t"+gf2+"\t"+g1+"\t"+g2+"\t"+oriG1+"\t"+oriG2+"\t"+str(dist)+"\t"+str(vscore)+"\t"+str(dscore)+"\t"+str(link)+"\n")

   scaff_gene_file.close() 
   print "DONE"



   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
