#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Compare annotation gene files between original and after mapping of genes on genome assembly computed by minia
###
###   INPUT:
###      1- CTG file
###         (data/data_DeCoSTAR/CTG_file)
###      2- DeCo* file results after linearization 
###         (results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_0.1_M1_kept)
###         (results/decostar/Xtopo-scaff/DeCoSTAR_Anopheles_Xtopo-scaff_Boltz_0.1_0.1_M1_kept)
###         (results/decostar/Xtopo_RAW/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1_0.1_M1_kept)
###         (results/decostar/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1_0.1_M1_kept)
###      3- Nx value (ex: 50 for N50 stats)
###         (50)
###      4- OUTPUT new ADJ file
###         (results/decostar/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_0.1_M1_new_adjacencies)
###         (results/decostar/DeCoSTAR_Anopheles_Xtopo-scaff_Boltz_0.1_0.1_M1_new_adjacencies)
###         (results/decostar/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1_0.1_M1_new_adjacencies)
###         (results/decostar/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1_0.1_M1_new_adjacencies)
###
###   OUTPUT:
###      - Give adjacencies that are in adjacencies file proposed by DeCo* (param 1) and are inconsistent with adajcencies in reference genome assembly (param 3)
###
###   Name: compute_scaffstats_and_newADJfile.py     Author: Yoann Anselmetti
###   Creation date: 2016/07/18                      Last modification: 2017/11/03
###

from sys import argv, stdout
from datetime import datetime
from re import search
from os import close, path, makedirs
from collections import namedtuple   #New in version 2.6
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


################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   contig_file=open(argv[1],'r')
   newAdj_file=open(argv[2],'r')
   x=argv[3]
   OUTPUT_file=argv[4]


   SIZE=namedtuple("SIZE",["bp","gene"])
   PAIR=namedtuple("PAIR",["g1","g2"])
   GENE=namedtuple("GENE",["spe","ctg","gf","gene","ori"])
   NEWADJ=namedtuple("NEWADJ",["spe","ctg1","ctg2","oriC1","oriC2","gf1","gf2","g1","g2","oriG1","oriG2","score"])


######################################
### STORE CTG INFOS IN dict_ID_ctg ###
######################################
   print "1/ Store CTG infos.. ",
   dict_spe_ctg_size=dict()
   dict_spe_gene=dict()
   dict_spe_gene_ctg=dict()
   dict_spe_ctg_gene=dict()
   dict_spe_assembly_size_bp=dict()
   dict_spe_assembly_size_gene=dict()
   # When get a contigs pairs (edge scaffolding link) => Get genes that are linked by scaffolding graph with the distance
   for line in contig_file:
      r2=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r2:
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
         exit("ERROR in line "+line+" of file "+argv[1]+" !!!")

         
   contig_file.close()
   print "DONE"


   # for spe in sorted(dict_spe_ctg_size):
   #    print spe
   #    for ctg in dict_spe_ctg_size[spe]:
   #       print "\t"+ctg+":\t",
   #       print dict_spe_ctg_size[spe][ctg]


######################
### Analysis of new extant adjacencies predicted by DeCo*
######################
   dict_spe_new_adj_tot=dict()
   dict_spe_new_adj_scaff=dict()
   dict_spe_newADJ=dict()
   for adj in newAdj_file:
      r=search('^([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*) ([^ ]*)[ \t]([^ \t\n]*)\n$', adj)
      if r:
         species=r.group(1)
         gene1=r.group(2)
         gene2=r.group(3)
         ori1=r.group(4)
         ori2=r.group(5)
         prior_score=r.group(6)
         post_score=r.group(7)
         scaff=r.group(8)


         species_name=gene1.split("@")[0]
         if species_name in dict_spe_gene:
            g1=gene1.split("@")[1]
            g2=gene2.split("@")[1]
            # If prior score!=1, then ADJ is a scaffolding one
            if float(prior_score)!=1.0:
               if not species_name in dict_spe_new_adj_tot:
                  dict_spe_new_adj_tot[species_name]=0
               dict_spe_new_adj_tot[species_name]+=1

               if float(prior_score)!=0.0:
                  if not species_name in dict_spe_new_adj_scaff:
                     dict_spe_new_adj_scaff[species_name]=0
                  dict_spe_new_adj_scaff[species_name]+=1


               ctg1=dict_spe_gene_ctg[species_name][g1]
               ctg2=dict_spe_gene_ctg[species_name][g2]
               size1=dict_spe_ctg_size[species_name][ctg1].bp
               size2=dict_spe_ctg_size[species_name][ctg2].bp
               geneNb1=dict_spe_ctg_size[species_name][ctg1].gene
               geneNb2=dict_spe_ctg_size[species_name][ctg2].gene

               size=size1+size2
               geneNb=geneNb1+geneNb2

               # Merge ctg1 and ctg2 in ctg1
               dict_spe_ctg_size[species_name].pop(ctg2,None)
               dict_spe_ctg_size[species_name][ctg1]=SIZE(size,geneNb)

               g1_ctg1=dict_spe_ctg_gene[species_name][ctg1].g1
               g2_ctg1=dict_spe_ctg_gene[species_name][ctg1].g2
               g1_ctg2=dict_spe_ctg_gene[species_name][ctg2].g1
               g2_ctg2=dict_spe_ctg_gene[species_name][ctg2].g2


               G1=""
               G2=""
               if g1==g1_ctg1:
                  G1=g2_ctg1
               else:
                  G1=g1_ctg1

               if g2==g1_ctg2:
                  G2=g2_ctg2
               else:
                  G2=g1_ctg2

               # print "BEFORE: "
               # print "\t-CTG1:",
               # print dict_spe_ctg_gene[species_name][ctg1]
               # print "\t-CTG2:",
               # print dict_spe_ctg_gene[species_name][ctg2]

               dict_spe_ctg_gene[species_name].pop(ctg2,None)
               dict_spe_ctg_gene[species_name][ctg1]=PAIR(G1,G2)

               # print "AFTER:",
               # print dict_spe_ctg_gene[species_name][ctg1]


               # THIS DOESN'T WORK FOR GNE THAT ARE ALONGE ON CONTIGS!!!
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

               
               new_adj=NEWADJ(species_name,ctg1,ctg2,oriC1,oriC2,gf1,gf2,g1,g2,oriG1,oriG2,post_score)
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
      list_size_bp=list()
      list_size_geneNb=list()
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
      print "\t- #scaffolds: "+str(scaffolds_nb)
      print "\t- #new adjacencies: "+str(new_adj)+" ("+str(new_scaff_adj)+" are scaff. adj.)"
      print "\t- Assembly size (bp): "+str(assembly_size_bp)
      print "\t- Assembly size (#gene): "+str(assembly_size_gene)
      print "\t- N"+str(x)+" (bp): "+str(Nx_bp)
      print "\t- N"+str(x)+" (#gene): "+str(Nx_gene)
      print ""


   Nx_bp_allspe=Nx(x,list_size_bp_allspe,assembly_size_bp_allspe)
   Nx_gene_allspe=Nx(x,list_size_gene_allspe,assembly_size_gene_allspe)

   print "For ALL species:"
   print "\t- #scaffolds: "+str(nb_scaff_tot)
   print "\t- #new adjacencies: "+str(new_adj_allspe)+ "("+str(new_scaff_adj_allspe)+" are scaff. adj.)"
   print "\t- Assembly size (bp): "+str(assembly_size_bp_allspe)
   print "\t- Assembly size (#gene): "+str(assembly_size_gene_allspe)
   print "\t- N"+str(x)+" (bp): "+str(Nx_bp_allspe)
   print "\t- N"+str(x)+" (#gene): "+str(Nx_gene_allspe)
   print ""




   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))