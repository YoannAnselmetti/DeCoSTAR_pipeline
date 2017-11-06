#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:                                                                    
###      Create scaffolding graph file for creation of instances for AD/ADseq
###
###   INPUT:                                                                   
###      1- Pre-scaffolding adjacencies file                                   
###         (data/data_DeCoSTAR/scaff_BESST_ALL_3_TRIMMOMATIC3)      
###      2- CTG file                                                           
###         (data/data_DeCoSTAR/CTG_file)                            
###      3- add gene to CTG scaffolds file                                     
###         (data/data_DeCoSTAR/scaff_BESST_DeCoSTAR)
###
###   OUTPUT:   (RUN in ~5min)                                                 
###      - Create scaffolding adj file for creation of Instances for ADseq     
###        => Take only 1 solution for each couple of gene adjacencies         
###           (Same if different orientations possible)                        
###                                                                            
###   Name: 14b-create_scaff_adj_file_final.py    Author: Yoann Anselmetti     
###   Creation date: 2015/11/06                   Last modification: 2017/03/15
###

from sys import argv, stdout
from re import search
from os import close
from collections import namedtuple   #New in version 2.6
from datetime import datetime



def get_gene_infos(ctg_order,ori,dist,CTG,contigEXT_file,preSCAFF_file):
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
      # Get correct orientation of gene2
   else:
      exit("ERROR on gene orientation in contigs extremities file "+contigEXT_file+" (Should be \"+\" or \"-\")")

   # Change gene orientation if CTG is in orientation "-"
   if (ori=="-"):
      oriG=ori_gene
      if (oriG=="-"):
         ori_gene="+"
      elif (oriG=="+"):
         ori_gene="-"
      elif (oriG=="?"):
         ori_gene="?"
      else:
         exit("ERROR on gene orientation in scaffolding graph file "+preSCAFF_file+" (Should be \"+\" or \"-\")")


   return GF_gene,gene,ori_gene,dist



if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   preSCAFF_file="data/data_DeCoSTAR/scaff_BESST_ALL_3_TRIMMOMATIC3"
   contigEXT_file="data/data_DeCoSTAR/CTG_file"
   OUTPUT_scaff_file="data/data_DeCoSTAR/scaff_BESST_DeCoSTAR"

   # STRUCTURE for gene edge and adjacency (EDGE==ADJ)
   ADJ=namedtuple("ADJ",["spe","g1","g2"])
   CTG=namedtuple("CTG",["spe","ctg","size","gf1","g1","ori1","start","gf2","g2","ori2","end"])
   EDGE=namedtuple("EDGE",["spe","ctg1","ctg2","oriC1","oriC2","gap","gf1","gf2","g1","g2","oriG1","oriG2","dist","vscore","dscore","link"])

######################################
### STORE CTG INFOS IN dict_spe_ID_ctg ###
######################################
   print "1/ Store CTG infos.. ",
   dict_spe_ID_ctg=dict()
   # When get a contigs pairs (edge scaffolding link) => Get genes that are linked by scaffolding graph with the distance
   with open(contigEXT_file,'r') as contig_file:
      for line in contig_file:
         spe=""
         contig=""
         contig_size=""
         GF1=""
         g1=""
         oriG1=""
         start_g1=""
         GF2=""
         g2=""
         oriG2=""
         stop_g2=""
         r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         r2=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         if r:
            spe=r.group(1)
            contig=r.group(2)
            contig_size=r.group(3)
            GF1=r.group(4)
            g1=r.group(5)
            oriG1=r.group(6)
            start_g1=r.group(7)
            GF2=r.group(8)
            g2=r.group(9)
            oriG2=r.group(10)
            stop_g2=r.group(11)
         elif r2:
            spe=r2.group(1)
            contig=r2.group(2)
            contig_size=r2.group(3)
            contig_geneNb=r2.group(4)
            GF1=r2.group(5)
            g1=r2.group(6)
            oriG1=r2.group(7)
            start_g1=r2.group(8)
            GF2=r2.group(9)
            g2=r2.group(10)
            oriG2=r2.group(11)
            stop_g2=r2.group(12)
         else:
            exit("ERROR in line "+line+" of file "+contigEXT_file+" !!!")

         if spe!="#species":
            if contig_size=="?":
               print "\n\t=> Contig "+contig+" is not present in FASTA file assembly of species "+spe
            else:
               if not spe in dict_spe_ID_ctg:
                  dict_spe_ID_ctg[spe]=dict()
               ctg=CTG(spe,contig,int(contig_size),GF1,g1,oriG1,int(start_g1),GF2,g2,oriG2,int(stop_g2))
               dict_spe_ID_ctg[spe][contig]=ctg

   contig_file.close()
   print "DONE"

   # for spe in dict_spe_ID_ctg:
   #    print spe+":"
   #       for ctg in dict_spe_ID_ctg[spe]:
   #          print "\t"+ctg+" :",
   #          print dict_spe_ID_ctg[spe][ctg]


###############################################################################################
### CREATE SCAFFOLDING GRAPH FILE WITH GENES INVOLVED IN SCAFFFOLDING EDGES INSTEAD CONTIGS ###
###############################################################################################
   print "2/ Store all scaffolding gene adjacencies present in file "+preSCAFF_file+"... ",
   stdout.flush()
   dict_spe_edge_scaff=dict()
   Nb_scaff_edge_tot=0
   Nb_scaff_edge_kept=0
   # Browse edge from scaffolding graph file and store them in "dict_spe_edge_scaff" (1 scaffolding graph file for all species)
   with open(preSCAFF_file,'r') as pre_scaff_file:
      for line in pre_scaff_file:
         r1=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         if r1:
            species=r1.group(1)
            ctg1=r1.group(2)
            ctg2=r1.group(3)
            ori1=r1.group(4)
            ori2=r1.group(5)
            gap=r1.group(6)
            vscore=r1.group(7)
            dscore=r1.group(8)
            links_nb=r1.group(9)

            if species!="#species":
               Nb_scaff_edge_tot+=1
               # distance between gene1 and gene2 => If gap is a negative distance == size of the overlap
               # -> distance between gene1 and gene2 = gap + d(3'gene1-ov) + d(ov-5'gene2)
               dist=float(gap)

               bool_ctg1=False
               bool_ctg2=False

               GF_gene1=""
               gene1=""
               ori_gene1=""
               GF_gene2=""
               gene2=""
               ori_gene2=""
               # If ctg1 is present in CTG_file: Get information on gene involved in the linked between ctg1 and ctg2
               if ctg1 in dict_spe_ID_ctg[species]:
                  bool_ctg1=True
                  CONTIG1=dict_spe_ID_ctg[species][ctg1]
                  # Get infos for gene1
                  GF_gene1,gene1,ori_gene1,dist = get_gene_infos("first",ori1,dist,CONTIG1,contigEXT_file,preSCAFF_file)

               # If ctg2 is present in contigEXT_file: Get information on gene involved in the linked between ctg1 and ctg2
               if ctg2 in dict_spe_ID_ctg[species]:
                  bool_ctg2=True
                  CONTIG2=dict_spe_ID_ctg[species][ctg2]
                  # Get infos for gene2
                  GF_gene2,gene2,ori_gene2,dist = get_gene_infos("second",ori2,dist,CONTIG2,contigEXT_file,preSCAFF_file)

               score1='{0:.12f}'.format(float(vscore))
               score2='{0:.12f}'.format(float(dscore))


               # If the 2 contigs are present in file contigEXT_file: Print scaffolding adj in OUTPUT_scaff_file
               if bool_ctg1 and bool_ctg2:
                  edge=EDGE(species,ctg1,ctg2,ori1,ori2,float(gap),GF_gene1,GF_gene2,gene1,gene2,ori_gene1,ori_gene2,float(dist),score1,score2,int(links_nb))
                  bool_adj=False
                  # If not species in "dict_spe_edge_scaff"
                  if not species in dict_spe_edge_scaff:
                     dict_spe_edge_scaff[species]=dict()
                     adj=ADJ(species,gene1,gene2)
                     dict_spe_edge_scaff[species][adj]=edge

                  else:
                     for adj in dict_spe_edge_scaff[species]:
                        # If the gene adajcences is ever present, then keep the adj with the besst vscore
                        if (gene1==adj.g1 and gene2==adj.g2) or (gene1==adj.g2 and gene2==adj.g1):
                           bool_adj=True
                           if dict_spe_edge_scaff[species][adj].vscore<score1:
                              dict_spe_edge_scaff[species][adj]=edge
                           elif dict_spe_edge_scaff[species][adj].vscore==score1:
                              # If vscore are equals take the adj with the best dscore
                              if dict_spe_edge_scaff[species][adj].dscore<score2:
                                 dict_spe_edge_scaff[species][adj]=edge
                              elif dict_spe_edge_scaff[species][adj].dscore==score2:
                                 # If dscore are equals take the adj with the higher number of links
                                 if dict_spe_edge_scaff[species][adj].link<int(links_nb):
                                    dict_spe_edge_scaff[species][adj]=edge
                           break
                     if not bool_adj:
                        adj=ADJ(species,gene1,gene2)
                        dict_spe_edge_scaff[species][adj]=edge
         else:
            exit("ERROR in line:\n\t"+line+"of file "+preSCAFF_file+" !!!")
   pre_scaff_file.close()
   print "DONE"


#########################################################
### WRITE OUTPUT_FILE => FINAL SCAFFOLDING GRAPH FILE ###
#########################################################
   print "3/ Write final scaffolding gene adjacencies file... ",
   stdout.flush()
   scaff_gene_file=open(OUTPUT_scaff_file,"w")
   scaff_gene_file.write("#species\tctg1\tctg2\torientation_ctg1\torientation_ctg2\tctg1-ctg2_dist\tgene1_family\tgene2_family\tgene1\tgene2\torientation_gene1\torientation_gene2\tgene1-gene2_dist\tvscore\tdscore\t#links\n")
   for species in sorted(dict_spe_edge_scaff):
      print species
      for adj in sorted(dict_spe_edge_scaff[species]):
         scaff_gene_file.write(dict_spe_edge_scaff[species][adj].spe+"\t"+dict_spe_edge_scaff[species][adj].ctg1+"\t"+dict_spe_edge_scaff[species][adj].ctg2+"\t"+dict_spe_edge_scaff[species][adj].oriC1+"\t"+dict_spe_edge_scaff[species][adj].oriC2+"\t"+str(dict_spe_edge_scaff[species][adj].gap)+"\t"+dict_spe_edge_scaff[species][adj].gf1+"\t"+dict_spe_edge_scaff[species][adj].gf2+"\t"+dict_spe_edge_scaff[species][adj].g1+"\t"+dict_spe_edge_scaff[species][adj].g2+"\t"+dict_spe_edge_scaff[species][adj].oriG1+"\t"+dict_spe_edge_scaff[species][adj].oriG2+"\t"+str(dict_spe_edge_scaff[species][adj].dist)+"\t"+str(dict_spe_edge_scaff[species][adj].vscore)+"\t"+str(dict_spe_edge_scaff[species][adj].dscore)+"\t"+str(dict_spe_edge_scaff[species][adj].link)+"\n")
   scaff_gene_file.close()
   print "DONE"

   edges_kept=0
   for species in dict_spe_edge_scaff:
      edges_kept+=len(dict_spe_edge_scaff[species])
   print "\n=> "+str(edges_kept)+"/"+str(Nb_scaff_edge_tot)+" edges have been kept (Others edges have at least one of the 2 linked contigs that have no gene in gene trees and have been erased OR linked several times the same couple fo genes (Keep the ADJ with best score))"

end_time = datetime.now()
print('Duration: {}'.format(end_time - start_time))
