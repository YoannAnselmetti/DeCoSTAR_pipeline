#! /usr/bin/env python2
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create scaffolds file in FASTA format file
###
###   INPUT:
###      1- New adjacencies file
###         (27avian_dataset/results/decostar/ADseq+scaff_Boltz_kT0.1/DeCoSTAR_9passeriformes_ADseq+scaff_Boltz_kT0.1_Linear_0.1_M1_new_adjacencies)
###      2- Genome assemblies directory
###         (27avian_dataset/data/INPUT_DATA/FASTA/SCAFF)
###      3- Directory containing FASTA files with 
###         (27avian_dataset/results/FASTA/SCAFF/ADseq+scaff_Boltz_kT0.1_Lin0.1_M1)
###
###   OUTPUT:
###      - Contig file with genes located to extremities
###
###   Name: create_FASTA_from_new_adj.py    Author: Yoann Anselmetti     
###   Creation date: 2018/06/11             Last modification: 2018/06/22
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
   # Get new ID of ctg1 and ctg2 if they heva been previously changed
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

   verbose=True

   OUTPUT_DIR=path.dirname(path.realpath(OUTPUT_SCAFF))
   # Create OUTPUT_DIR if not existing
   mkdir_p(OUTPUT_DIR)


   dict_spe_SCAFF=dict()
   FASTA_list=listdir(INPUT_SCAFF_dir)
   print "1/ Storing SCAFF sequences of species:"
   # Browse list of Genome assemblies to FASTA file format
   for FASTA in sorted(FASTA_list):
      i=0
      species=""
      r=search("^([^\.]*)\..*$",FASTA)
      if r:
         species=r.group(1)
         dict_spe_SCAFF[species]=dict()
      else:
         exit("!!! ERROR, FASTA file name: "+FASTA+" is incorrectly written !!!")

      print "\t"+species

      # Browse FASTA file of current species to get list of scaffolds 
      fasta_file=open(INPUT_SCAFF_dir+"/"+FASTA)
      for sequence in SeqIO.parse(fasta_file,"fasta"):
         scaff=sequence.id
         seq=sequence.seq
         dict_spe_SCAFF[species][scaff]=seq
      fasta_file.close()


   print "2/ Merge scaffolds linked by linearized new adjacencies:"
   dict_newCTG,dict_newCTG_ID=dict(),dict()
   input_file=open(newAdj_file,"r")
   for new_adj in input_file:
      r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",new_adj)
      if r:
         species=r.group(1)
         ctg1=r.group(2)
         ctg2=r.group(3)
         ori_ctg1=r.group(4)
         ori_ctg2=r.group(5)
         gf1=r.group(6)
         gf2=r.group(7)
         g1=r.group(8)
         g2=r.group(9)
         ori_g1=r.group(10)
         ori_g2=r.group(11)
         support=r.group(12)

         if species!="#species":
            mergeCTG(species,ctg1,ctg2,ori_ctg1,ori_ctg2,dict_newCTG_ID,dict_newCTG)


   # Print new linked scaffolds / species 
   for species in sorted(dict_newCTG):
      print "\n"+species+":"
      for ctg in dict_newCTG[species]:
         print "\t"+ctg
         # for elem in dict_newCTG[species][ctg]:
         #    print "\t\t",
         #    print elem


   print "\n3/ Write new FASTA files after scaffolding with linearized new adjacencies predicted by DeCoSTAR:"
   for species in sorted(dict_newCTG):
      output_scaff=open(OUTPUT_SCAFF+species+".fasta","w")
      for scaff in dict_newCTG[species]:
         output_scaff.write("> "+scaff+"\n")
         listCTG=dict_newCTG[species][scaff]
         seq=""
         for ctg in listCTG:
            if seq:
               output_scaff.write("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
            ID=ctg.id
            ori=ctg.ori
            seq=dict_spe_SCAFF[species][ID]
            if ori=="+":
               output_scaff.write(str(seq))
            elif ori=="-":
               # Write the REVERSE of the sequence (NOT THE REVERSE COMPLEMENT)
               for i in range(1,len(str(seq))+1):
                  output_scaff.write(str(seq)[-i])
            else:
               exit("ERROR, CTG orientation should be equal to \"-\" or \"+\" and not "+ori+" !!!")
         output_scaff.write("\n")
      output_scaff.close()


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
