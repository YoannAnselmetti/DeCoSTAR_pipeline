#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Compare annotation gene files between original and after mapping of genes on genome assembly computed by minia
###
###   INPUT:
###      1- CTG file
###         (18Anopheles_dataset/data/data_DeCoSTAR/CTG_file)
###      2- File containing the kept adjacencies after linearization 
###         (18Anopheles_dataset/results/decostar/Xtopo/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_Lin0.1_M1_kept)
###      3- Separator
###         (@)
###      3- Nx value (ex: 50 for N50 stats)
###         (50)
###      4- OUTPUT new ADJ file
###         (18Anopheles_dataset/results/decostar/Xtopo/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_Lin0.1_M1_new_extant_adjacencies)
###
###   OUTPUT:
###      - Give adjacencies that are in adjacencies file proposed by DeCo* (param 1) and are inconsistent with adajcencies in reference genome assembly (param 3)
###
###   Name: compute_scaffstats_and_newADJfile.py     Author: Yoann Anselmetti
###   Creation date: 2016/07/18                      Last modification: 2018/06/22
###

from sys import argv, stdout
from datetime import datetime
from re import search
from os import close, path, makedirs
from collections import namedtuple   #New in version 2.6
from copy import deepcopy
import errno


def mkdir_p(dir_path):
   output_dir=path.dirname(path.abspath(dir_path))
   try:
      makedirs(output_dir)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(output_dir):
         pass
      else: raise


def Nx(x,list_size,assembly_size):
   sum_size=0
   Nx=0
   for size in sorted(list_size,reverse=True):
      # print gene
      sum_size+=size
      if sum_size>=float(assembly_size/(100.0/float(x))):
         Nx=size
         break
   return Nx


def oriCTG(init_ori,new_ori):
   if (init_ori=="-" or init_ori=="+") and (new_ori=="-" or new_ori=="+"):
      if init_ori==new_ori:
         return "+"
      else:
         return "-"
   else:
      exit("ERROR, gene orientation is incorrect. Should be \"+\" or \"-\" and not \""+init_ori+" or "+new_ori+"\" !!!")


def store_CTG(CTG_file):
   dict_spe_ctg_size,dict_spe_gene,dict_spe_gene_ctg,dict_spe_ctg_gene,dict_spe_assembly_size_bp,dict_spe_assembly_size_gene=dict(),dict(),dict(),dict(),dict(),dict()
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
            if contig_size=="?":
               contig_size="10000"

            # Compute assembly size (bp) for each species
            if not spe in dict_spe_assembly_size_bp:
               dict_spe_assembly_size_bp[spe]=0
            dict_spe_assembly_size_bp[spe]+=int(contig_size)

            # Compute assembly size (#gene) for each species
            if not spe in dict_spe_assembly_size_gene:
               dict_spe_assembly_size_gene[spe]=0
            dict_spe_assembly_size_gene[spe]+=int(contig_geneNb)

            # Association CTG and gene localised to extremities for each species
            if not spe in dict_spe_ctg_gene:
               dict_spe_ctg_gene[spe]=dict()
            dict_spe_ctg_gene[spe][contig]=PAIR(g1,g2)

            # Store CTG size (bp and #gene) for each species
            if not spe in dict_spe_ctg_size:
               dict_spe_ctg_size[spe]=dict()
            dict_spe_ctg_size[spe][contig]=SIZE(int(contig_size),int(contig_geneNb))

            # Store association gene-CTG for each species
            if not spe in dict_spe_gene_ctg:
               dict_spe_gene_ctg[spe]=dict()
            dict_spe_gene_ctg[spe][g1]=contig
            dict_spe_gene_ctg[spe][g2]=contig

            if not spe in dict_spe_gene:
               dict_spe_gene[spe]=dict()

            gene1=GENE(spe,contig,GF1,g1,oriG1)
            gene2=GENE(spe,contig,GF2,g2,oriG2)
            dict_spe_gene[spe][g1]=gene1
            dict_spe_gene[spe][g2]=gene2

      else:
         exit("ERROR, line "+line+" of file "+CTG_file+" is incorrectly written!!!\nIt should match with the following format:\n"+CTG_format)

   contig_file.close()
   return dict_spe_ctg_size,dict_spe_gene,dict_spe_gene_ctg,dict_spe_ctg_gene,dict_spe_assembly_size_bp,dict_spe_assembly_size_gene




################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   CTG_file=argv[1]
   newAdj_file=open(argv[2],'r')
   sep=argv[3]
   x=argv[4]
   OUTPUT_file=argv[5]

   verbose=False

   SIZE=namedtuple("SIZE",["bp","gene"])
   PAIR=namedtuple("PAIR",["g1","g2"])
   GENE=namedtuple("GENE",["spe","ctg","gf","gene","ori"])
   NEWADJ=namedtuple("NEWADJ",["spe","ctg1","ctg2","oriC1","oriC2","gf1","gf2","g1","g2","oriG1","oriG2","score"])


######################################
### STORE CTG INFOS IN dict_ID_ctg ###
######################################
   print "1/ Store CTG infos...",
   dict_spe_ctg_size,dict_spe_gene,dict_spe_gene_ctg,dict_spe_ctg_gene,dict_spe_assembly_size_bp,dict_spe_assembly_size_gene=dict(),dict(),dict(),dict(),dict(),dict()
   dict_spe_ctg_size,dict_spe_gene,dict_spe_gene_ctg,dict_spe_ctg_gene,dict_spe_assembly_size_bp,dict_spe_assembly_size_gene=store_CTG(CTG_file)
   print "DONE"



######################
### Analysis of new extant adjacencies predicted by DeCo*
######################
   dict_spe_new_adj_tot,dict_spe_new_adj_scaff,dict_spe_newADJ=dict(),dict(),dict()

   dict_spe_ctg_size_INIT=deepcopy(dict_spe_ctg_size)
   dict_spe_ctg_gene_INIT=deepcopy(dict_spe_ctg_gene)
   # dict_spe_gene_ctg_INIT=dict_spe_gene_ctg.copy()

   for adj in newAdj_file:
      r=search('^([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*)[ \t]([^ \t\n]*)\n$', adj)
      if r:
         # print adj
         species=r.group(1)
         gene1=r.group(2)
         gene2=r.group(3)
         ori1=r.group(4)
         ori2=r.group(5)
         prior_score=r.group(6)
         post_score=r.group(7)
         scaff=r.group(8)


         species_name=gene1.split(sep)[0]
         if species_name in dict_spe_gene:
            g1=gene1.split(sep)[1]
            g2=gene2.split(sep)[1]
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

               # Get infos on contig of genes "g1" and "g2" ("ctg1" and "ctg2")
               ctg1=dict_spe_gene_ctg[species_name][g1]
               ctg2=dict_spe_gene_ctg[species_name][g2]
               size1=dict_spe_ctg_size[species_name][ctg1].bp
               size2=dict_spe_ctg_size[species_name][ctg2].bp
               geneNb1=dict_spe_ctg_size[species_name][ctg1].gene
               geneNb2=dict_spe_ctg_size[species_name][ctg2].gene

               # Sum sizes of ctg1 and ctg2
               size=size1+size2
               geneNb=geneNb1+geneNb2

               # Merge ctg1 and ctg2 in ctg1
               dict_spe_ctg_size[species_name].pop(ctg2,None)
               dict_spe_ctg_size[species_name][ctg1]=SIZE(size,geneNb)

               # Modify the genes to the extremities of "ctg1" ("ctg1" and "ctg2" merging)
               g1_ctg1=dict_spe_ctg_gene[species_name][ctg1].g1
               g2_ctg1=dict_spe_ctg_gene[species_name][ctg1].g2
               g1_ctg2=dict_spe_ctg_gene[species_name][ctg2].g1
               g2_ctg2=dict_spe_ctg_gene[species_name][ctg2].g2

               G1,G2="",""
               if g1==g1_ctg1:
                  G1=g2_ctg1
               else:
                  G1=g1_ctg1

               if g2==g1_ctg2:
                  G2=g2_ctg2
               else:
                  G2=g1_ctg2


               if verbose:
                  print "BEFORE: "
                  print "\t-CTG1:",
                  print dict_spe_ctg_gene[species_name][ctg1]
                  print "\t-CTG2:",
                  print dict_spe_ctg_gene[species_name][ctg2]

               dict_spe_ctg_gene[species_name].pop(ctg2,None)
               dict_spe_ctg_gene[species_name][ctg1]=PAIR(G1,G2)

               if verbose:
                  print "AFTER:",
                  print dict_spe_ctg_gene[species_name][ctg1]


               # THIS DOESN'T WORK FOR GENE THAT ARE ALONE ON CONTIGS!!!
               # dict_spe_gene_ctg[species_name].pop(g1,None)
               # dict_spe_gene_ctg[species_name].pop(g2,None)
               dict_spe_gene_ctg[species_name][G2]=ctg1


               # To print new adajcencies predict by the method
               oriG1=dict_spe_gene[species_name][g1].ori
               oriG2=dict_spe_gene[species_name][g2].ori
               oriC1=oriCTG(ori1,oriG1)
               oriC2=oriCTG(ori2,oriG2)
               ctg1=dict_spe_gene[species_name][g1].ctg
               ctg2=dict_spe_gene[species_name][g2].ctg
               gf1=dict_spe_gene[species_name][g1].gf
               gf2=dict_spe_gene[species_name][g2].gf
               g1=dict_spe_gene[species_name][g1].gene
               g2=dict_spe_gene[species_name][g2].gene

               
               new_adj=NEWADJ(species_name,ctg1,ctg2,oriC1,oriC2,gf1,gf2,g1,g2,ori1,ori2,post_score)
               if not species_name in dict_spe_newADJ:
                  dict_spe_newADJ[species_name]=list()
               dict_spe_newADJ[species_name].append(new_adj)



   newAdj_file.close()

   # Write new adjacencies file
   output_file=open(OUTPUT_file,'w')
   output_file.write("#species\tcontig1\tcontig2\torientation_contig1\torientation_contig2\tgene_family1\tgene_family2\tgene1\tgene2\torientation_gene1\torientation_gene2\tposterior_score\n")
   for species in sorted(dict_spe_newADJ):
      for adj in sorted(dict_spe_newADJ[species]):
         output_file.write(adj.spe+"\t"+adj.ctg1+"\t"+adj.ctg2+"\t"+adj.oriC1+"\t"+adj.oriC2+"\t"+adj.gf1+"\t"+adj.gf2+"\t"+adj.g1+"\t"+adj.g2+"\t"+adj.oriG1+"\t"+adj.oriG2+"\t"+adj.score+"\n")
   output_file.close()


   list_size_bp_allspe=list()
   list_size_gene_allspe=list()
   assembly_size_bp_allspe=0
   assembly_size_gene_allspe=0
   nb_scaff_tot=0
   new_adj_allspe=0
   new_scaff_adj_allspe=0
   for spe in sorted(dict_spe_ctg_size):
      list_size_bp,list_size_geneNb=list(),list()
      for size in dict_spe_ctg_size[spe].values():
         # print size
         list_size_bp.append(size.bp)
         list_size_bp_allspe.append(size.bp)
         list_size_geneNb.append(size.gene)
         list_size_gene_allspe.append(size.gene)

      assembly_size_bp=dict_spe_assembly_size_bp[spe]
      assembly_size_bp_allspe+=assembly_size_bp
      assembly_size_gene=dict_spe_assembly_size_gene[spe]
      assembly_size_gene_allspe+=assembly_size_gene

      Nx_bp=Nx(x,list_size_bp,assembly_size_bp)
      Nx_gene=Nx(x,list_size_geneNb,assembly_size_gene)

      scaffolds_nb=len(dict_spe_ctg_gene[spe].keys())
      nb_scaff_tot+=scaffolds_nb

      new_adj=0
      new_scaff_adj=0
      if spe in dict_spe_new_adj_tot:
         new_adj=dict_spe_new_adj_tot[spe]      
      if spe in dict_spe_new_adj_scaff:
         new_scaff_adj=dict_spe_new_adj_scaff[spe]

      new_adj_allspe+=new_adj
      new_scaff_adj_allspe+=new_scaff_adj

      print "For species "+spe+":"
      print "\t- Assembly size (bp): "+str(assembly_size_bp)
      print "\t- Assembly size (#gene): "+str(assembly_size_gene)
      print ""
      print "\tBEFORE DeCoSTAR:"
      print "\t\t- #scaffolds: "+str(len(dict_spe_ctg_gene_INIT[spe].keys()))
      list_size_bp,list_size_geneNb=list(),list()
      for size in dict_spe_ctg_size_INIT[spe].values():
         list_size_bp.append(size.bp)
         list_size_geneNb.append(size.gene)
      print "\t\t- N"+str(x)+" (bp): "+str(Nx(x,list_size_bp,assembly_size_bp))
      print "\t\t- N"+str(x)+" (#gene): "+str(Nx(x,list_size_geneNb,assembly_size_gene))
      print ""
      print "\tAFTER DeCoSTAR:"
      print "\t\t- #scaffolds: "+str(scaffolds_nb)
      print "\t\t- #new adjacencies: "+str(new_adj)+" ("+str(new_scaff_adj)+" are scaff. adj.)"
      print "\t\t- N"+str(x)+" (bp): "+str(Nx_bp)
      print "\t\t- N"+str(x)+" (#gene): "+str(Nx_gene)
      print "\n"


   Nx_bp_allspe=Nx(x,list_size_bp_allspe,assembly_size_bp_allspe)
   Nx_gene_allspe=Nx(x,list_size_gene_allspe,assembly_size_gene_allspe)
   print "For ALL species:"
   print "\t- Assembly size (bp): "+str(assembly_size_bp_allspe)
   print "\t- Assembly size (#gene): "+str(assembly_size_gene_allspe)
   print ""
   print "\tBEFORE DeCoSTAR:"
   nb_scaff_tot_INIT=0
   for spe in dict_spe_ctg_gene_INIT:
      nb_scaff_tot_INIT+=len(dict_spe_ctg_gene_INIT[spe].keys())
   print "\t\t- #scaffolds: "+str(nb_scaff_tot_INIT)
   list_size_bp,list_size_geneNb=list(),list()
   for spe in dict_spe_ctg_size_INIT:
      for size in dict_spe_ctg_size_INIT[spe].values():
         list_size_bp.append(size.bp)
         list_size_geneNb.append(size.gene)
   print "\t\t- N"+str(x)+" (bp): "+str(Nx(x,list_size_bp,assembly_size_bp_allspe))
   print "\t\t- N"+str(x)+" (#gene): "+str(Nx(x,list_size_geneNb,assembly_size_gene_allspe))
   print ""
   print "\tAFTER DeCoSTAR:"
   print "\t\t- #scaffolds: "+str(nb_scaff_tot)
   print "\t\t- #new adjacencies: "+str(new_adj_allspe)+ "("+str(new_scaff_adj_allspe)+" are scaff. adj.)"
   print "\t\t- N"+str(x)+" (bp): "+str(Nx_bp_allspe)
   print "\t\t- N"+str(x)+" (#gene): "+str(Nx_gene_allspe)


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))