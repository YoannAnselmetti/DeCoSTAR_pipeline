#! /usr/bin/env python2
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create new genome assemblies file with adjacencies predicted with DeCOSTAR (in FASTA format file)
###
###   INPUT:
###      1- New adjacencies file
###         (27avian_dataset/results/decostar/ADseq+scaff_Boltz_kT0.1/DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Linear_0.1_M1_new_adjacencies_with_scaff)
###      2- Genome assemblies directory
###         (27avian_dataset/data/INPUT_DATA/FASTA/SCAFF)
###      3- Directory containing FASTA files with 
###         (27avian_dataset/results/FASTA/SCAFF/ADseq+scaff_Boltz_kT0.1_Lin0.1_M1_)
###
###   OUTPUT:
###      - New genome assemblies FASTA files including new adjacencies predicted by DeCoSTAR
###
###   Name: create_FASTA_from_new_adjacencies.py    Author: Yoann Anselmetti     
###   Creation date: 2018/06/11                     Last modification: 2018/08/28
###

from sys import argv
from re import search
from os import close, path, listdir, makedirs
from datetime import datetime
import errno
from collections import namedtuple   #New in version 2.6

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC



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


def rev_ctg_order(listCTG):
   new_listCTG=list()
   for i in range(1,len(listCTG)+1):
      rev_ctg=CTG(listCTG[-i].id,rev_ori(listCTG[-i].ori))
      new_listCTG.append(rev_ctg)
   return new_listCTG


def mergeCTG(species,ctg1,ctg2,ori_ctg1,ori_ctg2,dict_newCTG_ID,dict_newCTG):
   # Get new ID of ctg1 and ctg2 if they have been previously changed
   new_ctg1,new_ctg2="",""
   if species in dict_newCTG_ID:
      if ctg1 in dict_newCTG_ID[species]:
         new_ctg1=dict_newCTG_ID[species][ctg1]
      if ctg2 in dict_newCTG_ID[species]:
         new_ctg2=dict_newCTG_ID[species][ctg2]

   # If ctg1 has been previously changed
   new_ctg1_order=list()
   new_ctg2_order=list()
   if new_ctg1:
      first_ctg1=dict_newCTG[species][new_ctg1][0]
      last_ctg1=dict_newCTG[species][new_ctg1][-1]
      # Get the correct order of contigs in new_ctg1
      if ctg1==first_ctg1.id:
         new_ctg1_order=rev_ctg_order(dict_newCTG[species][new_ctg1])
      elif ctg1==last_ctg1.id:
         new_ctg1_order=dict_newCTG[species][new_ctg1]
      else:
         exit("ERROR, contig \""+ctg1+"\" should be to the start or the end of the new contig \""+new_ctg1+"\" !!!")
      dict_newCTG[species].pop(new_ctg1,None)

   # If ctg2 has been previously changed
   if new_ctg2:        
      first_ctg2=dict_newCTG[species][new_ctg2][0]
      last_ctg2=dict_newCTG[species][new_ctg2][-1]
      # Get the correct order of contigs in new_27avian_dataset/data/INPUT_DATA/27avian_species.t2
      if ctg2==last_ctg2.id:
         new_ctg2_order=rev_ctg_order(dict_newCTG[species][new_ctg2])
      elif ctg2==first_ctg2.id:
         new_ctg2_order=dict_newCTG[species][new_ctg2]
      else:
         exit("ERROR, contig \""+ctg2+"\" should be to the start or the end of the new contig \""+new_ctg2+"\" !!!")
      dict_newCTG[species].pop(new_ctg2,None)

   # If ctg1 has NOT been previously changed
   if not new_ctg1_order:
      new_ctg1_order.append(CTG(ctg1,ori_ctg1))
   # If ctg2 has NOT been previously changed
   if not new_ctg2_order:
      new_ctg2_order.append(CTG(ctg2,ori_ctg2))

   # Concatenate the two CTG order to define the new SCAFF
   new_ctg_order=new_ctg1_order+new_ctg2_order

   listCTG=list()
   new_CTG_ID=""
   for elem in new_ctg_order:
      ctg=elem.id
      ori=elem.ori
      listCTG.append(ctg)
      if not new_CTG_ID:
         new_CTG_ID+=ctg+"("+ori+")"
      else:
         new_CTG_ID+=":"+ctg+"("+ori+")"

   if not species in dict_newCTG:
      dict_newCTG[species]=dict()
   dict_newCTG[species][new_CTG_ID]=new_ctg_order

   if not species in dict_newCTG_ID:
      dict_newCTG_ID[species]=dict()
   for ctg in listCTG:
      dict_newCTG_ID[species][ctg]=new_CTG_ID






################
###   MAIN 
################
if __name__ == '__main__':

   start_time = datetime.now()

   CTG=namedtuple("CTG",["id","ori"])

   # Recovery of input parameters
   newAdj_file=argv[1]
   INPUT_SCAFF_dir=argv[2]
   OUTPUT_SCAFF=argv[3]

   default_gap_size=100

   verbose=1

   OUTPUT_DIR=path.dirname(path.realpath(OUTPUT_SCAFF))
   # Create OUTPUT_DIR if not existing
   mkdir_p(OUTPUT_DIR)


   print "1/ Merge scaffolds linked by linearized new adjacencies"
   dict_newCTG,dict_newCTG_ID,dict_distCTG=dict(),dict(),dict()
   input_file=open(newAdj_file,"r")
   for new_adj in input_file:
      r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",new_adj)
      if r:
         species=r.group(1)
         ctg1=r.group(2)
         ctg2=r.group(3)
         ori_ctg1=r.group(4)
         ori_ctg2=r.group(5)
         dist_ctg=r.group(6)
         gf1=r.group(7)
         gf2=r.group(8)
         g1=r.group(9)
         g2=r.group(10)
         ori_g1=r.group(11)
         ori_g2=r.group(12)
         dist_gene=r.group(13)
         vscore=r.group(14)
         dscore=r.group(15)
         links=r.group(16)
         support=r.group(17)

         if species!="#species":
            dict_distCTG[(ctg1,ctg2)]=dist_ctg
            mergeCTG(species,ctg1,ctg2,ori_ctg1,ori_ctg2,dict_newCTG_ID,dict_newCTG)
   input_file.close()


   # Print the new CTG ID after merging initial contigs/scaffoldsin new scaffolds 
   if verbose>1:
      # Print new linked scaffolds / species 
      for species in sorted(dict_newCTG):
         print "\n"+species+":"
         for ctg in dict_newCTG[species]:
            print "\t"+ctg
            if verbose>2:
               for elem in dict_newCTG[species][ctg]:
                  print "\t\t",
                  print elem


   print "\n2/ Write new FASTA files after scaffolding with linearized new adjacencies predicted by DeCoSTAR"
   for species in sorted(dict_newCTG):
      if verbose>0:
         print "\t"+species
      output_scaff=open(OUTPUT_SCAFF+species+".fasta","w")

      dict_SCAFF=dict()
      FASTA_FILE=""
      for FASTA in sorted(listdir(INPUT_SCAFF_dir)):
         i=0
         spe=""
         r=search("^([^\.]*)\..*$",FASTA)
         if r:
            spe=r.group(1)
         else:
            exit("!!! ERROR, FASTA file name: "+FASTA+" is incorrectly written !!!")
         if spe==species:
            FASTA_FILE=INPUT_SCAFF_dir+"/"+FASTA
            # Browse FASTA file of current species to get list of scaffolds 
            fasta_file=open(INPUT_SCAFF_dir+"/"+FASTA)
            for sequence in SeqIO.parse(fasta_file,"fasta"):
               scaff=sequence.id
               seq=sequence.seq
               dict_SCAFF[scaff]=seq
            fasta_file.close()
            break

      for scaff in dict_newCTG[species]:
         output_scaff.write("> "+scaff+"\n")
         listCTG=dict_newCTG[species][scaff]
         seq=""
         ID=""
         nb_N=0
         for ctg in listCTG:
            # Write "N" gap between the 2 contigs scaffold
            if seq:
               if (ID,ctg.id) in dict_distCTG:
                  if dict_distCTG[(ID,ctg.id)]=="?":
                     nb_N=default_gap_size
                  else:
                     nb_N=int(float(dict_distCTG[(ID,ctg.id)]))
                     # print nb_N
                     if nb_N<0:
                        if verbose>1:
                           print "NEGATIVE distance between contigs "+ID+" and "+ctg.id
                        nb_N=default_gap_size
               elif (ctg.id,ID) in dict_distCTG:
                  if dict_distCTG[(ctg.id,ID)]=="?":
                     nb_N=default_gap_size
                  else:
                     nb_N=int(float(dict_distCTG[(ctg.id,ID)]))
                     # print nb_N
                     if nb_N<0:
                        if verbose>1:
                           print "NEGATIVE distance between contigs "+ctg.id+" and "+ID
                        nb_N=default_gap_size
               else:
                  exit("ERROR: CTG adjacency ("+ID+"-"+ctg.id+") is not present in DeCoSTAR predicted adjacencies file: "+newAdj_file)
               i=0
               while i<nb_N:
                  output_scaff.write("N")
                  i+=1

            # Write sequence of current contig
            ID=ctg.id
            ori=ctg.ori
            seq=dict_SCAFF[ID]
            del dict_SCAFF[ID]

            if ori=="+":
               output_scaff.write(str(seq))
            elif ori=="-":
               # Write the REVERSE of the sequence (NOT THE REVERSE COMPLEMENT)
               for i in range(1,len(str(seq))+1):
                  output_scaff.write(str(seq)[-i])
            else:
               exit("ERROR, CTG orientation should be equal to \"-\" or \"+\" and not "+ori+" !!!")
         output_scaff.write("\n")
      # Write remaining scaffolds in the FASTA file
      for scaff in dict_SCAFF:
         output_scaff.write("> "+scaff+"\n")
         output_scaff.write(str(dict_SCAFF[scaff])+"\n")
      output_scaff.close()


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time-start_time))
