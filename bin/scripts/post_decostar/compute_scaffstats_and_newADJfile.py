#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Compute scaffolding statistics and generate new extant adjacencies file 
###
###   INPUT:
###      1- INPUT FASTA files directory
###         (18Anopheles_dataset/data/INPUT_DATA/FASTA/SCAFF)
###      2- CTG file
###         (18Anopheles_dataset/data/data_DeCoSTAR/CTG_file)
###      3- File containing the kept adjacencies after linearization 
###         (18Anopheles_dataset/results/decostar/Xtopo/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_Lin0.1_M1_kept)
###      4- File containing the discarded adjacencies after linearization 
###         (18Anopheles_dataset/results/decostar/Xtopo/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_Lin0.1_M1_disc)
###      5- Separator
###         (@)
###      6- OUTPUT new ADJ file
###         (18Anopheles_dataset/results/decostar/Xtopo/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_Lin0.1_M1_new_extant_adjacencies)
###      7- OUTPUT scaffolding stats directory
###         (18Anopheles_dataset/results/scaffolding_stats/Xtopo/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_Lin0.1_M1)
###
###   OUTPUT:
###      - Produce adjacencies file containing new extant adjacencies proposed by DeCo*
###      - Produce file and graphics summarizing scaffolding statistics of DeCo*
###
###   Name: compute_scaffstats_and_newADJfile.py     Author: Yoann Anselmetti
###   Creation date: 2016/07/18                      Last modification: 2019/01/16
###

from sys import argv, stdout
from datetime import datetime
from re import search
from os import close, path, makedirs, listdir
from collections import namedtuple   #New in version 2.6
from copy import deepcopy
import errno
from Bio import SeqIO
import numpy as np 
from matplotlib import pyplot as plt



def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else:
         raise



def Nx_Lx(x,list_size,assembly_size):
   sum_size=0
   Nx=0
   Lx=0
   for size in sorted(list_size,reverse=True):
      # print gene
      sum_size+=size
      Lx+=1
      if sum_size>=float(assembly_size/(100.0/float(x))):
         Nx=size
         break
   return Nx,Lx



def oriCTG(init_ori,new_ori):
   if (init_ori=="-" or init_ori=="+") and (new_ori=="-" or new_ori=="+"):
      if init_ori==new_ori:
         return "+"
      else:
         return "-"
   else:
      exit("ERROR, gene orientation is incorrect. Should be \"+\" or \"-\" and not \""+init_ori+" or "+new_ori+"\" !!!")



def get_CTGsize(FASTA_dir):
   dict_CTGsize=dict()
   for FASTA in sorted(listdir(FASTA_dir)):
      species=""
      r=search("^([^\.]*)\..*$",FASTA)
      if r:
         species=r.group(1)
         dict_CTGsize[species]=dict()
      else:
         exit("!!! ERROR, FASTA file name: "+FASTA+" is incorrectly written (Must be: ${species_name}.fa) !!!")
      FASTA_FILE=FASTA_dir+"/"+FASTA
      if verbose>0:
         print "\tStore sequences size of FASTA file \""+FASTA_FILE+"\""
      # Browse FASTA file of current species to get list of scaffolds 
      fasta_file=open(FASTA_dir+"/"+FASTA)
      for sequence in SeqIO.parse(fasta_file,"fasta"):
         ID=sequence.id
         size=len(sequence.seq)
         dict_CTGsize[species][ID]=size
      fasta_file.close()
   return dict_CTGsize



def store_CTG(CTG_file,dict_CTGsize):
   # dict_spe_ctg_size,dict_spe_gene,dict_spe_ctg_gene=dict(),dict(),dict()
   dict_spe_ctg_size,dict_spe_gene=dict(),dict()
   CTG_format="#species\tctg\tctg_size\tctg_gene_nb\t5'_gene_family\t5'_gene\torientation_5'_gene\tstart_5'_gene\t3'_gene_family\t3'_gene\torientation_3'_gene\tend_3'_gene\n"
   # When get a contigs pairs (edge scaffolding link) => Get genes pair linked by scaffolding with the distance
   contig_file=open(CTG_file,'r')
   for line in contig_file:
      r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r:
         species=r.group(1)
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

         if species!="#species":
            if contig_size=="?":
               print "!!! WARNING: Size of contig/scaffold "+contig+" is unknown !!! -> Size set to \"10,000 bp\""
               contig_size="10000"

            # Test if CTGsize(FASTA file)==CTGsize(CTG file)
            if int(contig_size)!=dict_CTGsize[species][contig]:
               exit("!!! ERROR, size of contig \""+contig+"\" is different between FASTA file ("+contig_size+") and CTG file ("+dict_CTGsize[species][contig]+") !!!")

            # Store CTG size (bp and #gene) for each species
            if not species in dict_spe_ctg_size:
               dict_spe_ctg_size[species]=dict()
            dict_spe_ctg_size[species][contig]=SIZE(int(contig_size),int(contig_geneNb))

            # # Association CTG and gene localised to extremities for each species
            # if not species in dict_spe_ctg_gene:
            #    dict_spe_ctg_gene[species]=dict()
            # dict_spe_ctg_gene[species][contig]=PAIR(g1,g2)

            # Store gene infos by species 
            if not species in dict_spe_gene:
               dict_spe_gene[species]=dict()
            gene1=GENE(species,contig,GF1,g1,oriG1)
            gene2=GENE(species,contig,GF2,g2,oriG2)
            dict_spe_gene[species][g1]=gene1
            dict_spe_gene[species][g2]=gene2

            # Remove current CTG from the dict_CTGsize
            if contig in dict_CTGsize[species]:
               dict_CTGsize[species].pop(contig,None)
      else:
         exit("!!! ERROR, line "+line+" of file "+CTG_file+" is incorrectly written!!!\nIt should match with the following format:\n"+CTG_format+" !!!")
   contig_file.close()

   # Add scaffolds without gene present in gene families (used as input of DeCoSTAR)
   for species in dict_CTGsize:
      for ctg in dict_CTGsize[species]:
         size=dict_CTGsize[species][ctg]
         dict_spe_ctg_size[species][ctg]=SIZE(size,0)

   # return dict_spe_ctg_size,dict_spe_gene,dict_spe_ctg_gene
   return dict_spe_ctg_size,dict_spe_gene



def mergeCTG(species,ctg1,ctg2,dict_newSCAFF_ID,dict_newSCAFF,dict_spe_ctg_size,scaff_id):
   # Get new ID of ctg1 and ctg2 if they have been previously changed
   CTG1,CTG2="",""
   new_ctg1,new_ctg2=False,False
   list_CTG=list()
   if species in dict_newSCAFF_ID:
      if ctg1 in dict_newSCAFF_ID[species]:
         CTG1=dict_newSCAFF_ID[species][ctg1]
         list_CTG+=dict_newSCAFF[species][CTG1]
         new_ctg1=True
      else:
         CTG1=ctg1
         list_CTG.append(CTG1)
      if ctg2 in dict_newSCAFF_ID[species]:
         CTG2=dict_newSCAFF_ID[species][ctg2]
         list_CTG+=dict_newSCAFF[species][CTG2]
         new_ctg2=True
      else:
         CTG2=ctg2
         list_CTG.append(CTG2)
   else:
      CTG1=ctg1
      CTG2=ctg2
      list_CTG.append(CTG1)
      list_CTG.append(CTG2)

   # Get size (bp & #gene) for CTG1 and CTG2 
   size1=dict_spe_ctg_size[species][CTG1].bp
   size2=dict_spe_ctg_size[species][CTG2].bp
   geneNb1=dict_spe_ctg_size[species][CTG1].gene
   geneNb2=dict_spe_ctg_size[species][CTG2].gene
   size=size1+size2
   geneNb=geneNb1+geneNb2

   # Remove CTG1 and CTG2 from dict_spe_ctg_size
   dict_spe_ctg_size[species].pop(CTG1,None)
   dict_spe_ctg_size[species].pop(CTG2,None)

   # Set the ID of the new scaffold
   newSCAFF_ID=""
   if new_ctg1:
      newSCAFF_ID=CTG1
   elif new_ctg2:
      newSCAFF_ID=CTG2
   else:
      scaff_id,newSCAFF_ID=get_SCAFF_ID(scaff_id)

   ##########
   ### Merge CTG1 and CTG2 in scaffold "newSCAFF_ID"
   ##########
   # Set the list of the scaffold "newSCAFF_ID"
   if not species in dict_newSCAFF:
      dict_newSCAFF[species]=dict()
   else:
      if CTG1 in dict_newSCAFF[species]:
         dict_newSCAFF[species].pop(CTG1,None)
      if CTG2 in dict_newSCAFF[species]:
         dict_newSCAFF[species].pop(CTG2,None)   
   dict_newSCAFF[species][newSCAFF_ID]=list_CTG

   # Set the "newSCAFF_ID" of CTG in "list_CTG"
   if not species in dict_newSCAFF_ID:
      dict_newSCAFF_ID[species]=dict()
   for ctg in list_CTG:
      dict_newSCAFF_ID[species][ctg]=newSCAFF_ID

   # Set the size of the scaffold "newSCAFF_ID"
   dict_spe_ctg_size[species][newSCAFF_ID]=SIZE(size,geneNb)

   return scaff_id,dict_spe_ctg_size



def get_SCAFF_ID(i):
   SCAFF_ID=""
   i+=1
   if i>=1000000 and i<=9999999:
      SCAFF_ID="SCAFF"+str(i)
   elif i>=100000:
      SCAFF_ID="SCAFF0"+str(i)
   elif i>=10000:
      SCAFF_ID="SCAFF00"+str(i)
   elif i>=1000:
      SCAFF_ID="SCAFF000"+str(i)
   elif i>=100:
      SCAFF_ID="SCAFF0000"+str(i)
   elif i>=10:
      SCAFF_ID="SCAFF00000"+str(i)
   elif i>=1:
      SCAFF_ID="SCAFF000000"+str(i)
   else:
      exit("ERROR too much scaffolds!!! -> (>9999999)")
   return i,SCAFF_ID



def print_scaffolding_stats(scaff_stats_DIR,x_list,dict_spe_ctg_size_INIT,dict_spe_ctg_size_FIN,dict_spe_new_adj_tot,dict_spe_new_adj_scaff):
   output_stats=open(scaff_stats_DIR+"/scaffolding_stats_file.tab","w")
   output_stats.write("#species\tassemblySIZE(bp)\tassemblySIZE(#gene)\t#newADJ(#scaffADJ)\t#scaffolds(initial)\t#scaffolds(final)\tminSIZE(initial)\tminSIZE(final)\tmaxSIZE(initial)\tmaxSIZE(final)")
   for x in x_list:
      output_stats.write("\tN"+str(x)+"_bp(initial)\tN"+str(x)+"_bp(final)\tN"+str(x)+"_bp(fold_improvment)\tL"+str(x)+"_bp(initial)\tL"+str(x)+"_bp(final)\tL"+str(x)+"_bp(fold_improvment)\tN"+str(x)+"_#gene(initial)\tN"+str(x)+"_#gene(final)\tN"+str(x)+"_#gene(fold_improvment)\tL"+str(x)+"_#gene(initial)\tL"+str(x)+"_#gene(final)\tL"+str(x)+"_#gene(fold_improvment)")
   output_stats.write("\n")

   dict_scaffolding_stats=dict()
   list_size_bp_allspe,list_size_geneNb_allspe,list_size_bp_allspe_INIT,list_size_geneNb_allspe_INIT=list(),list(),list(),list()
   assembly_size_bp_allspe,assembly_size_gene_allspe=0,0   
   nb_scaff_tot,nb_scaff_tot_INIT,new_adj_allspe,new_scaff_adj_allspe=0,0,0,0
   minSCAFF_allspe,maxSCAFF_allspe,minSCAFF_allspe_INIT,maxSCAFF_allspe_INIT=100000000000000000000000,0,100000000000000000000000,0
   for spe in sorted(dict_spe_ctg_size_FIN):
      assembly_size_bp,assembly_size_gene=0,0
      ##########
      ### STATS IN INITIAL GENOME ASSEMBLIES (BEFORE DECOSTAR)
      ##########
      list_size_bp_INIT,list_size_geneNb_INIT=list(),list()
      minSCAFF_INIT,maxSCAFF_INIT=100000000000000000000000,0
      for size in dict_spe_ctg_size_INIT[spe].values():
         # print size
         # Get value of the longest scaffold
         if size.bp>maxSCAFF_allspe_INIT:
            maxSCAFF_INIT=size.bp
            maxSCAFF_allspe_INIT=size.bp
         elif size.bp>maxSCAFF_INIT:
            maxSCAFF_INIT=size.bp
         # Get value of the shortest scaffold
         if size.bp<minSCAFF_allspe_INIT:
            minSCAFF_INIT=size.bp
            minSCAFF_allspe_INIT=size.bp
         elif size.bp<minSCAFF_INIT:
            minSCAFF_INIT=size.bp
         # Store values of the scaffold length (in bp and #gene(->present in gene families))
         list_size_bp_INIT.append(size.bp)
         list_size_bp_allspe_INIT.append(size.bp)
         list_size_geneNb_INIT.append(size.gene)
         list_size_geneNb_allspe_INIT.append(size.gene)
         # Get values of the genome assembly (in bp and #gene(->present in gene families))
         assembly_size_bp+=size.bp
         assembly_size_gene+=size.gene
      assembly_size_bp_allspe+=assembly_size_bp
      assembly_size_gene_allspe+=assembly_size_gene
      scaffolds_nb_INIT=len(dict_spe_ctg_size_INIT[spe].keys())
      nb_scaff_tot_INIT+=scaffolds_nb_INIT

      ##########
      ### STATS IN INITIAL GENOME ASSEMBLIES (BEFORE DECOSTAR)
      ##########
      list_size_bp,list_size_geneNb=list(),list()
      minSCAFF,maxSCAFF=100000000000000000000000,0
      for size in dict_spe_ctg_size_FIN[spe].values():
         # print size
         # Get value of the longest scaffold
         if size.bp>maxSCAFF_allspe:
            maxSCAFF=size.bp
            maxSCAFF_allspe=size.bp
         elif size.bp>maxSCAFF:
            maxSCAFF=size.bp
         # Get value of the shortest scaffold
         if size.bp<minSCAFF_allspe:
            minSCAFF=size.bp
            minSCAFF_allspe=size.bp
         elif size.bp<minSCAFF:
            minSCAFF=size.bp
         # Store values of the scaffold length (in bp and #gene(->present in gene families))
         list_size_bp.append(size.bp)
         list_size_bp_allspe.append(size.bp)
         list_size_geneNb.append(size.gene)
         list_size_geneNb_allspe.append(size.gene)

      scaffolds_nb=len(dict_spe_ctg_size_FIN[spe].keys())
      nb_scaff_tot+=scaffolds_nb

      new_adj=0
      new_scaff_adj=0
      if spe in dict_spe_new_adj_tot:
         new_adj=dict_spe_new_adj_tot[spe]      
      if spe in dict_spe_new_adj_scaff:
         new_scaff_adj=dict_spe_new_adj_scaff[spe]

      new_adj_allspe+=new_adj
      new_scaff_adj_allspe+=new_scaff_adj

      print "\nFor species "+spe+":"
      print "\t- Assembly size (bp): "+str(assembly_size_bp)
      print "\t- Assembly size (#gene): "+str(assembly_size_gene)
      print ""
      print "\tBEFORE DeCoSTAR:"
      print "\t\t- #scaffolds: "+str(scaffolds_nb_INIT)
      print "\t\t- min(size): "+str(minSCAFF_INIT),"bp"
      print "\t\t- max(size): "+str(maxSCAFF_INIT),"bp"
      print "\tAFTER DeCoSTAR:"
      print "\t\t- #scaffolds: "+str(scaffolds_nb)
      print "\t\t- min(size): "+str(minSCAFF),"bp"
      print "\t\t- max(size): "+str(maxSCAFF),"bp"
      print "\t\t- #new adjacencies: "+str(new_adj)+" ("+str(new_scaff_adj)+" are scaff. adj.)"
      print ""
      output_stats.write(spe+"\t"+str(assembly_size_bp)+"\t"+str(assembly_size_gene)+"\t"+str(new_adj)+"("+str(new_scaff_adj)+")\t"+str(scaffolds_nb_INIT)+"\t"+str(scaffolds_nb)+"\t"+str(minSCAFF_INIT)+"\t"+str(minSCAFF)+"\t"+str(maxSCAFF_INIT)+"\t"+str(maxSCAFF))

      # Loop on x for Nx and Lx stats
      dict_scaffolding_stats[spe]=dict()
      for x in x_list:
         print "\tN"+str(x)+"/L"+str(x)+" statistics:"
         Nx_bp_INIT,Lx_bp_INIT=Nx_Lx(x,list_size_bp_INIT,assembly_size_bp)
         Nx_gene_INIT,Lx_gene_INIT=Nx_Lx(x,list_size_geneNb_INIT,assembly_size_gene)
         print "\t\tBEFORE DeCoSTAR:"
         print "\t\t\t- N"+str(x)+" (bp): "+str(Nx_bp_INIT)
         print "\t\t\t- L"+str(x)+" (bp): "+str(Lx_bp_INIT)
         print "\t\t\t- N"+str(x)+" (#gene): "+str(Nx_gene_INIT)
         print "\t\t\t- L"+str(x)+" (#gene): "+str(Lx_gene_INIT)
         print ""
         Nx_bp_FIN,Lx_bp_FIN=Nx_Lx(x,list_size_bp,assembly_size_bp)
         Nx_gene_FIN,Lx_gene_FIN=Nx_Lx(x,list_size_geneNb,assembly_size_gene)
         Nx_bp_fold=round(float(Nx_bp_FIN)/float(Nx_bp_INIT),2)
         Lx_bp_fold=round(float(Lx_bp_INIT)/float(Lx_bp_FIN),2)
         Nx_gene_fold=round(float(Nx_gene_FIN)/float(Nx_gene_INIT),2)
         Lx_gene_fold=round(float(Lx_gene_INIT)/float(Lx_gene_FIN),2)
         print "\t\tAFTER DeCoSTAR:"
         print "\t\t\t- N"+str(x)+" (bp): "+str(Nx_bp_FIN)
         print "\t\t\t\t+ fold improvment: "+str(Nx_bp_fold)
         print "\t\t\t- L"+str(x)+" (bp): "+str(Lx_bp_FIN)
         print "\t\t\t\t+ fold improvment: "+str(Lx_bp_fold)
         print "\t\t\t- N"+str(x)+" (#gene): "+str(Nx_gene_FIN)
         print "\t\t\t\t+ fold improvment: "+str(Nx_gene_fold)
         print "\t\t\t- L"+str(x)+" (#gene): "+str(Lx_gene_FIN)
         print "\t\t\t\t+ fold improvment: "+str(Lx_gene_fold)
         print ""
         output_stats.write("\t"+str(Nx_bp_INIT)+"\t"+str(Nx_bp_FIN)+"\t"+str(Nx_bp_fold)+"\t"+str(Lx_bp_INIT)+"\t"+str(Lx_bp_FIN)+"\t"+str(Lx_bp_fold)+"\t"+str(Nx_gene_INIT)+"\t"+str(Nx_gene_FIN)+"\t"+str(Nx_gene_fold)+"\t"+str(Lx_gene_INIT)+"\t"+str(Lx_gene_FIN)+"\t"+str(Lx_gene_fold))
         dict_scaffolding_stats[spe][x]=SCAFFOLDING_STATS(scaffolds_nb_INIT,scaffolds_nb,Nx_bp_INIT,Nx_gene_INIT,Lx_bp_INIT,Lx_gene_INIT,Nx_bp_FIN,Nx_gene_FIN,Lx_bp_FIN,Lx_gene_FIN)
      output_stats.write("\n")

   print "\nFor ALL species:"
   print "\t- Assembly size (bp): "+str(assembly_size_bp_allspe)
   print "\t- Assembly size (#gene): "+str(assembly_size_gene_allspe)
   print ""
   print "\tBEFORE DeCoSTAR:"
   print "\t\t- #scaffolds: "+str(nb_scaff_tot_INIT)
   print "\t\t- min(size): "+str(minSCAFF_allspe_INIT),"bp"
   print "\t\t- max(size: "+str(maxSCAFF_allspe_INIT),"bp"
   print "\tAFTER DeCoSTAR:"
   print "\t\t- #scaffolds: "+str(nb_scaff_tot)
   print "\t\t- min(size): "+str(minSCAFF_allspe),"bp"
   print "\t\t- max(size): "+str(maxSCAFF_allspe),"bp"
   print "\t\t- #new adjacencies: "+str(new_adj_allspe)+ "("+str(new_scaff_adj_allspe)+" are scaff. adj.)"
   print ""
   output_stats.write("All species\t"+str(assembly_size_bp_allspe)+"\t"+str(assembly_size_gene_allspe)+"\t"+str(new_adj_allspe)+"("+str(new_scaff_adj_allspe)+")\t"+str(nb_scaff_tot_INIT)+"\t"+str(nb_scaff_tot)+"\t"+str(minSCAFF_allspe_INIT)+"\t"+str(minSCAFF_allspe)+"\t"+str(maxSCAFF_allspe_INIT)+"\t"+str(maxSCAFF_allspe))
   # Loop on x for Nx and Lx stats
   # dict_scaffolding_stats["ALL species"]=dict()
   for x in x_list:
      print "\tN"+str(x)+"/L"+str(x)+" statistics:"
      Nx_bp_allspe_INIT,Lx_bp_allspe_INIT=Nx_Lx(x,list_size_bp_allspe_INIT,assembly_size_bp_allspe)
      Nx_gene_allspe_INIT,Lx_gene_allspe_INIT=Nx_Lx(x,list_size_geneNb_allspe_INIT,assembly_size_gene_allspe)
      print "\t\tBEFORE DeCoSTAR:"
      print "\t\t\t- N"+str(x)+" (bp): "+str(Nx_bp_allspe_INIT)
      print "\t\t\t- L"+str(x)+" (bp): "+str(Lx_bp_allspe_INIT)
      print "\t\t\t- N"+str(x)+" (#gene): "+str(Nx_gene_allspe_INIT)
      print "\t\t\t- L"+str(x)+" (#gene): "+str(Lx_gene_allspe_INIT)
      print ""
      Nx_bp_allspe_FIN,Lx_bp_allspe_FIN=Nx_Lx(x,list_size_bp_allspe,assembly_size_bp_allspe)
      Nx_gene_allspe_FIN,Lx_gene_allspe_FIN=Nx_Lx(x,list_size_geneNb_allspe,assembly_size_gene_allspe)
      Nx_bp_allspe_fold=round(float(Nx_bp_allspe_FIN)/float(Nx_bp_allspe_INIT),2)
      Lx_bp_allspe_fold=round(float(Lx_bp_allspe_INIT)/float(Lx_bp_allspe_FIN),2)
      Nx_gene_allspe_fold=round(float(Nx_gene_allspe_FIN)/float(Nx_gene_allspe_INIT),2)
      Lx_gene_allspe_fold=round(float(Lx_gene_allspe_INIT)/float(Lx_gene_allspe_FIN),2)
      print "\t\tAFTER DeCoSTAR:"
      print "\t\t\t- N"+str(x)+" (bp): "+str(Nx_bp_allspe_FIN)
      print "\t\t\t\t+ fold improvment: "+str(Nx_bp_allspe_fold)
      print "\t\t\t- L"+str(x)+" (bp): "+str(Lx_bp_allspe_FIN)
      print "\t\t\t\t+ fold improvment: "+str(Lx_bp_allspe_fold)
      print "\t\t\t- N"+str(x)+" (#gene): "+str(Nx_gene_allspe_FIN)
      print "\t\t\t\t+ fold improvment: "+str(Nx_gene_allspe_fold)
      print "\t\t\t- L"+str(x)+" (#gene): "+str(Lx_gene_allspe_FIN)
      print "\t\t\t\t+ fold improvment: "+str(Lx_gene_allspe_fold)
      output_stats.write("\t"+str(Nx_bp_allspe_INIT)+"\t"+str(Nx_bp_allspe_FIN)+"\t"+str(Nx_bp_allspe_fold)+"\t"+str(Lx_bp_allspe_INIT)+"\t"+str(Lx_bp_allspe_FIN)+"\t"+str(Lx_bp_allspe_fold)+"\t"+str(Nx_gene_allspe_INIT)+"\t"+str(Nx_gene_allspe_FIN)+"\t"+str(Nx_gene_allspe_fold)+"\t"+str(Lx_gene_allspe_INIT)+"\t"+str(Lx_gene_allspe_FIN)+"\t"+str(Lx_gene_allspe_fold))
      # dict_scaffolding_stats["ALL species"][x]=SCAFFOLDING_STATS(Nx_bp_allspe_INIT,Nx_gene_allspe_INIT,Lx_bp_allspe_INIT,Lx_gene_allspe_INIT,Nx_bp_allspe,Nx_gene_allspe,Lx_bp_allspe,Lx_gene_allspe)
   output_stats.write("\n")

   plot_scaffolding_stats(x_list,dict_scaffolding_stats,scaff_stats_DIR)



def plot_scaffolding_stats(x_list,dict_scaffolding_stats,scaff_stats_DIR):
   for x in x_list:
      list_species=list()
      list_scaffolds_nb_INIT,list_scaffolds_nb_FIN=list(),list()
      list_Nx_bp_INIT,list_Nx_bp_FIN,list_Nx_gene_INIT,list_Nx_gene_FIN=list(),list(),list(),list()
      list_Lx_bp_INIT,list_Lx_bp_FIN,list_Lx_gene_INIT,list_Lx_gene_FIN=list(),list(),list(),list()
      for species in sorted(dict_scaffolding_stats):
         short_species=species[0]+"."+species.split("_")[1]
         list_species.append(short_species)
         scaff_stats=dict_scaffolding_stats[species][x]
         list_scaffolds_nb_INIT.append(scaff_stats.nb_scaffolds_INIT)
         list_scaffolds_nb_FIN.append(scaff_stats.nb_scaffolds_FIN)
         list_Nx_bp_INIT.append(scaff_stats.Nx_bp_INIT)
         list_Nx_bp_FIN.append(scaff_stats.Nx_bp_FIN)
         list_Nx_gene_INIT.append(scaff_stats.Nx_gene_INIT)
         list_Nx_gene_FIN.append(scaff_stats.Nx_gene_FIN)
         list_Lx_bp_INIT.append(scaff_stats.Lx_bp_INIT)
         list_Lx_bp_FIN.append(scaff_stats.Lx_bp_FIN)
         list_Lx_gene_INIT.append(scaff_stats.Lx_gene_INIT)
         list_Lx_gene_FIN.append(scaff_stats.Lx_gene_FIN)
      # GRAPH to compare INITIAL and FINAl stats
      graph_INITIAL_vs_FINAL_lower(list_species,list_scaffolds_nb_INIT,list_scaffolds_nb_FIN,scaff_stats_DIR,"#scaffolds")
      graph_INITIAL_vs_FINAL_upper(list_species,list_Nx_bp_INIT,list_Nx_bp_FIN,scaff_stats_DIR,"N"+str(x)+"(bp)")
      graph_INITIAL_vs_FINAL_upper(list_species,list_Nx_gene_INIT,list_Nx_gene_FIN,scaff_stats_DIR,"N"+str(x)+"(#gene)")
      graph_INITIAL_vs_FINAL_lower(list_species,list_Lx_bp_INIT,list_Lx_bp_FIN,scaff_stats_DIR,"L"+str(x)+"(bp)")
      graph_INITIAL_vs_FINAL_lower(list_species,list_Lx_gene_INIT,list_Lx_gene_FIN,scaff_stats_DIR,"L"+str(x)+"(#gene)")

      graph_INITIAL_vs_FINAL_Nx_nbSCAFF(list_species,list_Nx_bp_INIT,list_Nx_bp_FIN,list_scaffolds_nb_INIT,list_scaffolds_nb_FIN,scaff_stats_DIR,"N"+str(x)+"(bp)","#scaffolds")



def graph_INITIAL_vs_FINAL_upper(list_species,list_INITIAL,list_FINAL,OUTPUT,TAG):
   i=0
   while i<len(list_species):
      plt.plot((list_INITIAL[i],list_FINAL[i]),(list_species[i],list_species[i]),color="darkgreen",linestyle="-",zorder=1)
      i+=1
   plt.scatter(list_INITIAL,list_species,color="darkcyan",marker="s",zorder=2)
   plt.scatter(list_FINAL,list_species,color="darkgreen",marker=">",zorder=2)
   plt.xlabel(TAG)
   plt.yticks(fontstyle="italic",fontweight="bold",rotation=15)
   fig_name=OUTPUT+"/"+TAG+".pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()



def graph_INITIAL_vs_FINAL_lower(list_species,list_INITIAL,list_FINAL,OUTPUT,TAG):
   i=0
   while i<len(list_species):
      plt.plot((list_INITIAL[i],list_FINAL[i]),(list_species[i],list_species[i]),color="saddlebrown",linestyle="-",zorder=1)
      i+=1
   plt.scatter(list_INITIAL,list_species,color="chocolate",marker="s",zorder=2)
   plt.scatter(list_FINAL,list_species,color="saddlebrown",marker="<",zorder=2)
   plt.xlabel(TAG)
   plt.yticks(fontstyle="italic",fontweight="bold",rotation=15)
   fig_name=OUTPUT+"/"+TAG+".pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()



def graph_INITIAL_vs_FINAL_Nx_nbSCAFF(list_species,list_INITIAL_Nx,list_FINAL_Nx,list_INITIAL_nbSCAFF,list_FINAL_nbSCAFF,OUTPUT,Nx_TAG,nbSCAFF_TAG):
   width=10.0
   path=20.0
   fig,ax1a=plt.subplots()

   list_X_ax1=np.arange(0.0-(width/2.0),float(len(list_species)*path)-(width/2.0),path)
   list_X_ax2=np.arange(0.0+(width/2.0),float(len(list_species)*path)+(width/2.0),path)

   ax1a.tick_params(axis="y",labelrotation=30)
   ax1a.set_xlabel(Nx_TAG,fontweight="bold",color="darkgreen")
   ax1a.scatter(list_INITIAL_Nx,list_X_ax1,color="darkcyan",marker="s",zorder=2)
   ax1a.scatter(list_FINAL_Nx,list_X_ax1,color="darkgreen",marker=">",zorder=2)

   ax2a=ax1a.twiny()  # instantiate a second axes that shares the same y-axis
   ax2a.set_xlabel(nbSCAFF_TAG,fontweight="bold",color="saddlebrown")
   ax2a.scatter(list_INITIAL_nbSCAFF,list_X_ax2,color="chocolate",marker="s",zorder=2)
   ax2a.scatter(list_FINAL_nbSCAFF,list_X_ax2,color="saddlebrown",marker="<",zorder=2)

   i=0
   while i<len(list_species):
      ax1a.plot((list_INITIAL_Nx[i],list_FINAL_Nx[i]),(list_X_ax1[i],list_X_ax1[i]),color="darkgreen",linestyle="-",zorder=1)
      ax2a.plot((list_INITIAL_nbSCAFF[i],list_FINAL_nbSCAFF[i]),(list_X_ax2[i],list_X_ax2[i]),color="saddlebrown",linestyle="-",zorder=1)
      i+=1

   plt.yticks(np.arange(0,len(list_species)*path,path),list_species,rotation=15,fontstyle="italic",fontweight="bold")

   fig_name=OUTPUT+"/"+Nx_TAG+"_"+nbSCAFF_TAG+".pdf"
   fig.tight_layout()
   fig.savefig(fig_name,format='pdf')
   plt.close(fig)



################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   FASTA_dir=argv[1]
   CTG_file=argv[2]
   kept_newAdj_file=open(argv[3],'r')
   disc_newAdj_file=open(argv[4],'r')
   sep=argv[5]
   OUTPUT_file=argv[6]
   scaff_stats_DIR=argv[7]

   x_list=[50,90]
   verbose=1

   SIZE=namedtuple("SIZE",["bp","gene"])
   PAIR=namedtuple("PAIR",["g1","g2"])
   GENE=namedtuple("GENE",["spe","ctg","gf","gene","ori"])
   NEWADJ=namedtuple("NEWADJ",["spe","ctg1","ctg2","oriC1","oriC2","gf1","gf2","g1","g2","oriG1","oriG2","score"])
   SCAFFOLDING_STATS=namedtuple("SCAFFOLDING_STATS",["nb_scaffolds_INIT","nb_scaffolds_FIN","Nx_bp_INIT","Nx_gene_INIT","Lx_bp_INIT","Lx_gene_INIT","Nx_bp_FIN","Nx_gene_FIN","Lx_bp_FIN","Lx_gene_FIN"])

   mkdir_p(scaff_stats_DIR)

######################
### STORE SEQUENCES SIZE FROM FASTA FILES
######################
   print "\n1/ Store sequences size of initial genome assemblies:"
   dict_CTGsize=get_CTGsize(FASTA_dir)


######################
### STORE CTG INFOS
######################
   print "\n2/ Store CTG infos...",
   # dict_spe_ctg_size_INIT,dict_spe_gene=dict(),dict()
   dict_spe_ctg_size_INIT,dict_spe_gene=store_CTG(CTG_file,dict_CTGsize)
   dict_CTGsize.clear()
   print "DONE"


######################
### MERGE CTG/SCAFF LINKED WITH NEW ADJACENCIES PREDICTED BY DECOSTAR
######################
   print "\n3/ Merge sequences with new adjacencies predicted by DeCo*:"
   dict_spe_new_adj_tot,dict_spe_new_adj_scaff,dict_spe_newADJ,dict_newSCAFF_ID,dict_newSCAFF=dict(),dict(),dict(),dict(),dict()
   dict_spe_ctg_size_FIN=deepcopy(dict_spe_ctg_size_INIT)
   scaff_id=0
   # Read new adjacencies file
   for adj in kept_newAdj_file:
      species,gene1,gene2,ori1,ori2,prior_score,post_score,scaff=adj.split()
      species_name=gene1.split(sep)[0]
      if species_name in dict_spe_gene:
         g1=gene1.split(sep)[1]
         g2=gene2.split(sep)[1]
         ########
         ### HUGE APPROXIMATION: THESE ADJACENCIES CAN BE OBSERVED ADJACENCIES IF WE DECIDE TO SET THE PRE-SCORE<1.0 !!!
         ########
         # If (prior_score!=1), then ADJ is a new one
         if float(prior_score)!=1.0:
            # print species_name
            # print gene1
            if not species_name in dict_spe_new_adj_tot:
               dict_spe_new_adj_tot[species_name]=0
            dict_spe_new_adj_tot[species_name]+=1

            if float(prior_score)!=0.0:
               if not species_name in dict_spe_new_adj_scaff:
                  dict_spe_new_adj_scaff[species_name]=0
               dict_spe_new_adj_scaff[species_name]+=1

            # Get CTG infos
            ctg1=dict_spe_gene[species_name][g1].ctg
            ctg2=dict_spe_gene[species_name][g2].ctg
            oriG1=dict_spe_gene[species_name][g1].ori
            oriG2=dict_spe_gene[species_name][g2].ori
            oriC1=oriCTG(ori1,oriG1)
            oriC2=oriCTG(ori2,oriG2)
            gf1=dict_spe_gene[species_name][g1].gf
            gf2=dict_spe_gene[species_name][g2].gf
            g1=dict_spe_gene[species_name][g1].gene
            g2=dict_spe_gene[species_name][g2].gene

            # Merge ctg1 and ctg2
            scaff_id,dict_spe_ctg_size_FIN=mergeCTG(species_name,ctg1,ctg2,dict_newSCAFF_ID,dict_newSCAFF,dict_spe_ctg_size_FIN,scaff_id)

            # Store new adjacencies
            new_adj=NEWADJ(species_name,ctg1,ctg2,oriC1,oriC2,gf1,gf2,g1,g2,ori1,ori2,post_score)
            if not species_name in dict_spe_newADJ:
               dict_spe_newADJ[species_name]=list()
            dict_spe_newADJ[species_name].append(new_adj)

   kept_newAdj_file.close()
   dict_newSCAFF_ID.clear()
   dict_newSCAFF.clear()


######################
### WRITE NEW ADJACENCIES FILE
######################
   print "\n4/ Write new adjacencies file \""+OUTPUT_file+"\":"
   # Write new adjacencies file
   output_file=open(OUTPUT_file,'w')
   output_file.write("#species\tcontig1\tcontig2\torientation_contig1\torientation_contig2\tgene_family1\tgene_family2\tgene1\tgene2\torientation_gene1\torientation_gene2\tposterior_score\n")
   for species in sorted(dict_spe_newADJ):
      for adj in sorted(dict_spe_newADJ[species]):
         output_file.write(adj.spe+"\t"+adj.ctg1+"\t"+adj.ctg2+"\t"+adj.oriC1+"\t"+adj.oriC2+"\t"+adj.gf1+"\t"+adj.gf2+"\t"+adj.g1+"\t"+adj.g2+"\t"+adj.oriG1+"\t"+adj.oriG2+"\t"+adj.score+"\n")
   output_file.close()
   dict_spe_newADJ.clear()


# ##########
# ### SPLIT CTG WHERE ADJ HAVE BEEN DISCARDED BY LINEARIZATION
# ##########
#    print "\n5/ Split contigs where initial adjacencies have been discarded by linearization:"
#    # Read discarded adjacencies file
#    for disc_adj in disc_newAdj_file:
#       species_ID,gene1,gene2,ori1,ori2,pre_score,post_score,linear_step=disc_adj.split()
#       if pre_score==post_score=="1.0":




######################
### COMPUTE GENOME SCAFFOLDING STATISTICS
######################
   print "\n6/ Compute genome scaffolding statistics:"
   print_scaffolding_stats(scaff_stats_DIR,x_list,dict_spe_ctg_size_INIT,dict_spe_ctg_size_FIN,dict_spe_new_adj_tot,dict_spe_new_adj_scaff)
   # for species in dict_spe_ctg_size_INIT:
   #    print "\t\tSpecies",species,":"
   #    for contig in dict_spe_ctg_size_INIT[species]:
   #       print "\t\t\t"+contig+":"
   #       print "\t\t\t\tINITIAL -> bp:",dict_spe_ctg_size_INIT[species][contig].bp,"| #gene",dict_spe_ctg_size_INIT[species][contig].gene
   #       if contig in dict_spe_ctg_size_FIN:
   #          print "\t\t\t\tFINAL -> bp:",dict_spe_ctg_size_FIN[species][contig].bp,"| #gene",dict_spe_ctg_size_FIN[species][contig].gene


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))