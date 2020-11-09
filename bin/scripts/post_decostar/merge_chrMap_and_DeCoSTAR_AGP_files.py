#! /usr/bin/env python3
# -*- coding: utf-8 -*-
###
###   Goal:
###      Merge DeCoSTAR and chromosome map AGP files
###
###   INPUT:
###      1- DeCoSTAR AGP file
###         (27avian_dataset/results/AGP/DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Lin0.1_M2_Zosterops_borbonicus.agp)
###      2- Chromosome map AGP file (can be corrected with DeCoSTAR predicitions)
###         (27avian_dataset/results/Zosterops_borbonicus_ZeFi-ZoBo+DeCoSTAR_final_version.agp)
###      3- OUTPUT AGP file
###         (27avian_dataset/results/AGP/Zosterops_borbonicus_ZeFi-ZoBo+DeCoSTAR_final_version.agp)
###
###   OUTPUT:
###      - AGP file merging DeCoSTAR and Chromosome map adjacencies
###
###   Name: merge_chrMap_and_DeCoSTAR_AGP_files.py    Author: Yoann Anselmetti     
###   Creation date: 2018/11/13                       Last modification: 2020/11/04
###

from sys import argv
from re import search
from os import close, path, listdir, makedirs
from datetime import datetime
import errno
from collections import namedtuple   #New in version 2.6
from Bio import SeqIO


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else:
         raise


def parse_DeCoSTAR_AGP_file(AGP_file,dict_newSCAFF,dict_newSCAFF_ID,dict_distCTG,dict_CTG):
   input_file=open(AGP_file,"r")
   gap_in_progress=False
   stored_obj,stored_ctg,stored_gapSize="","",""
   for line in input_file:
      obj,objBeg,objEnd,partNumb,compoType,compoId_gapLength,compoBeg_gapType,compoEnd_link,orientation_linkEvidence=line.split()
      # If line corresponds to a sequence (contig, scaffold, chromosome, ...)
      if compoType in ['A','D','F','G','O','P','W']:
         dict_CTG[compoId_gapLength]=int(compoEnd_link)
         if gap_in_progress:
            gap_in_progress=False
            dict_distCTG[(stored_ctg,compoId_gapLength)]=stored_gapSize
         else:
            if stored_obj==obj:
               dict_distCTG[(stored_ctg,compoId_gapLength)]="?"
         if not obj in dict_newSCAFF:
            dict_newSCAFF[obj]=list()
         dict_newSCAFF[obj].append(CTG(compoId_gapLength,orientation_linkEvidence))
         dict_newSCAFF_ID[compoId_gapLength]=obj
         stored_ctg=compoId_gapLength
         stored_obj=obj
      # If line corresponds to a gap
      elif compoType in ['N','U']:
         gap_in_progress=True
         if compoType=='N':
            stored_gapSize=compoId_gapLength
         else:
            stored_gapSize="?"
      else:
         exit("ERROR, column 5 should be equal to ['A','D','F','G','O','P','W','N','U'] and not '"+compoType+"' (cf. AGP file format: https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/#FORMAT)")


def parse_chrMAP_AGP_file(AGP_file,dict_newSCAFF,dict_newSCAFF_ID,dict_distCTG):
   input_file=open(AGP_file,"r")
   gap_in_progress=False
   stored_obj,stored_ctg,stored_gapSize="","",""
   for line in input_file:
      # print "line:",line
      obj,objBeg,objEnd,partNumb,compoType,compoId_gapLength,compoBeg_gapType,compoEnd_link,orientation_linkEvidence=line.split()
      # If line corresponds to a sequence (contig, scaffold, chromosome, ...)
      if compoType in ['A','D','F','G','O','P','W']:
         if gap_in_progress:
            gap_in_progress=False
            if (not (stored_ctg,compoId_gapLength) in dict_distCTG) and (not (compoId_gapLength,stored_ctg) in dict_distCTG):
               dict_distCTG[(stored_ctg,compoId_gapLength)]=stored_gapSize
         else:
            if stored_obj==obj:
               if (not (stored_ctg,compoId_gapLength) in dict_distCTG) and (not (compoId_gapLength,stored_ctg) in dict_distCTG):
                  dict_distCTG[(stored_ctg,compoId_gapLength)]="?"
         if compoId_gapLength in dict_newSCAFF_ID:
            new_scaff=dict_newSCAFF_ID[compoId_gapLength]
            if new_scaff in dict_newSCAFF:
               del dict_newSCAFF[new_scaff]
         if not obj in dict_newSCAFF:
            dict_newSCAFF[obj]=list()
         dict_newSCAFF[obj].append(CTG(compoId_gapLength,orientation_linkEvidence))
         dict_newSCAFF_ID[compoId_gapLength]=obj
         stored_ctg=compoId_gapLength
         stored_obj=obj
      # If line corresponds to a gap
      elif compoType in ['N','U']:
         gap_in_progress=True
         if compoType=='N':
            stored_gapSize=compoId_gapLength
         else:
            stored_gapSize="?"
      else:
         exit("ERROR, column 5 should be equal to ['A','D','F','G','O','P','W','N','U'] and not '"+compoType+"' (cf. AGP file format: https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/#FORMAT)")


################
###   MAIN 
################
if __name__ == '__main__':

   start_time = datetime.now()

   CTG=namedtuple("CTG",["id","ori"])

   # Recovery of input parameters
   AGP_decostar_file=argv[1]
   AGP_chrMAP_file=argv[2]
   output_AGP_file=argv[3]

   default_gap_size=100
   verbose=1

   # # Create OUTPUT_DIR if not existing
   # mkdir_p(OUTPUT_SCAFF)


   # Store the DeCoSTAR AGP file
   dict_newSCAFF,dict_newSCAFF_ID,dict_distCTG,dict_CTG=dict(),dict(),dict(),dict()
   parse_DeCoSTAR_AGP_file(AGP_decostar_file,dict_newSCAFF,dict_newSCAFF_ID,dict_distCTG,dict_CTG)
   if verbose>1:
      print("\nDeCoSTAR:")
      print("\tsize(dict_newSCAFF):",len(dict_newSCAFF))
      print("\tsize(dict_newSCAFF_ID):",len(dict_newSCAFF_ID))
      print("\tsize(dict_distCTG):",len(dict_distCTG))
      print("\tsize(dict_CTG):",len(dict_CTG))


   parse_chrMAP_AGP_file(AGP_chrMAP_file,dict_newSCAFF,dict_newSCAFF_ID,dict_distCTG)
   if verbose>1:
      print("\nDeCoSTAR+Chromosome Map:")
      print("\tsize(dict_newSCAFF):",len(dict_newSCAFF))
      print("\tsize(dict_newSCAFF_ID):",len(dict_newSCAFF_ID))
      print("\tsize(dict_distCTG):",len(dict_distCTG))
      print("\tsize(dict_CTG):",len(dict_CTG))


   output_agp=open(output_AGP_file,"w")
   for scaff in sorted(dict_newSCAFF):
      listCTG=dict_newSCAFF[scaff]
      ID=""
      gap_size=0
      posSTART=1
      posEND=0
      bool_default=False
      for ctg in listCTG:
         ctg_size=dict_CTG[ctg.id]
         del dict_CTG[ctg.id]
         # Get size of the gap between the 2 scaffolded contigs (ID and ctg.id)
         if ID:
            if (ID,ctg.id) in dict_distCTG:
               if dict_distCTG[(ID,ctg.id)]=="?":
                  bool_default=True
                  gap_size=default_gap_size
               else:
                  gap_size=int(float(dict_distCTG[(ID,ctg.id)]))
                  # print gap_size
                  if gap_size<0:
                     if verbose>1:
                        print("\t!!!WARNING!!! => NEGATIVE distance between contigs "+ID+" and "+ctg.id)
            elif (ctg.id,ID) in dict_distCTG:
               if dict_distCTG[(ctg.id,ID)]=="?":
                  bool_default=True
                  gap_size=default_gap_size
               else:
                  gap_size=int(float(dict_distCTG[(ctg.id,ID)]))
                  # print gap_size
                  if gap_size<0:
                     if verbose>1:
                        print("\t!!!WARNING!!! => NEGATIVE distance between contigs "+ctg.id+" and "+ID)
            else:
               exit("ERROR: CTG adjacency ("+ID+"-"+ctg.id+") is not present in DeCoSTAR predicted adjacencies file: "+AGP_decostar_file)

            # Write gap in AGP file
            posEND+=gap_size
            if bool_default:
               output_agp.write(scaff+"\t"+str(posSTART)+"\t"+str(posEND)+"\t.\tU\t"+str(gap_size)+"\tscaffold\tno\tna\n")
            else:
               output_agp.write(scaff+"\t"+str(posSTART)+"\t"+str(posEND)+"\t.\tN\t"+str(gap_size)+"\tscaffold\tyes\tpaired-ends\n")
            bool_default=False
            posSTART=posEND+1
         posEND+=ctg_size

         # Write current contig in the AGP file
         ID=ctg.id
         ori=ctg.ori
         output_agp.write(scaff+"\t"+str(posSTART)+"\t"+str(posEND)+"\t.\tW\t"+ID+"\t1\t"+str(ctg_size)+"\t"+ori+"\n")
         posSTART=posEND+1

   # Write remaining contigs (not scaffolded) in the AGP file
   for ctg in sorted(dict_CTG):
      ctg_size=dict_CTG[ctg]
      output_agp.write(ctg+"\t1\t"+str(ctg_size)+"\t.\tW\t"+ctg+"\t1\t"+str(ctg_size)+"\t+\n")
   output_agp.close()


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time-start_time))
