#! /usr/bin/env python3
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create new genome assemblies file with adjacencies predicted with DeCoSTAR (in FASTA format file)
###
###   INPUT:
###      1- AGP file
###         (27avian_dataset/results/AGP/)
###      2- prefix file
###         (DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Lin0.1_M2_)
###      3- Genome assemblies directory
###         (27avian_dataset/data/INPUT_DATA/FASTA/SCAFF)
###      4- Directory containing FASTA files containing new scaffolds 
###         (27avian_dataset/results/FASTA/SCAFF/DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Lin0.1_M2_)
###
###   OUTPUT:
###      - New genome assembly FASTA files including new adjacencies predicted by DeCoSTAR
###
###   Name: create_newFASTA_from_AGP.py    Author: Yoann Anselmetti
###   Creation date: 2018/06/11            Last modification: 2020/11/04
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


def parse_AGP_file(AGP_file,species):
   input_file=open(AGP_file,"r")
   gap_in_progress=False
   stored_obj,stored_ctg,stored_gapSize="","",""
   for line in input_file:
      obj,objBeg,objEnd,partNumb,compoType,compoId_gapLength,compoBeg_gapType,compoEnd_link,orientation_linkEvidence=line.split()
      # If line corresponds to a sequence (contig, scaffold, chromosome, ...)
      if compoType in ['A','D','F','G','O','P','W']:
         if gap_in_progress:
            gap_in_progress=False
            dict_distCTG[(stored_ctg,compoId_gapLength)]=stored_gapSize
         else:
            if stored_obj==obj:
               dict_distCTG[(stored_ctg,compoId_gapLength)]="?"
         if not obj in dict_newSCAFF:
            dict_newSCAFF[obj]=list()
         dict_newSCAFF[obj].append(CTG(compoId_gapLength,orientation_linkEvidence))
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
   AGP_dir=argv[1]
   file_prefix=argv[2]
   INPUT_SCAFF_dir=argv[3]
   OUTPUT_SCAFF=argv[4]

   write_unscaffolded_ctg=True
   default_gap_size=100
   verbose=1

   # Create OUTPUT_DIR if not existing
   mkdir_p(OUTPUT_SCAFF)

   for AGP_file in sorted(listdir(AGP_dir)):
      if ".agp" in AGP_file:
         species=AGP_file.split(file_prefix)[1].split(".agp")[0]
         print("\nFor species "+species+":")

      print("\t1/ Parse AGP file to get CTG order and gaps sizes")
      dict_newSCAFF,dict_distCTG=dict(),dict()
      parse_AGP_file(AGP_dir+"/"+AGP_file,species)


      # Print the new CTG ID after merging initial contigs/scaffolds in new scaffolds 
      if verbose>1:
         if verbose>2:
            print("\n\t1bis/ Print association between new scaffolds ID and old scaffolds/contigs ID:")
         else:
            print("\n\t1bis/ Print new scaffolds ID:")
         # Print new linked scaffolds
         for scaff in sorted(dict_newSCAFF):
            print("\t\t"+ctg)
            if verbose>2:
               for elem in dict_newSCAFF[scaff]:
                  print("\t\t\t", end=' ')
                  print(elem)

      # Print distance between CTG pairs linked by DeCoSTAR
      if verbose>2:
         print("\n\t1ter/ Print distance between CTG pairs linked by DeCoSTAR:")
         for distCTG in sorted(dict_distCTG):
            print(distCTG,"\t",dict_distCTG[distCTG])


      print("\n\t2/ Write new FASTA files after scaffolding with linearized new adjacencies predicted by DeCoSTAR")
      output_scaff=open(OUTPUT_SCAFF+"/"+file_prefix+species+".fasta","w")
      dict_SCAFF=dict()
      FASTA_FILE=""
      for FASTA in sorted(listdir(INPUT_SCAFF_dir)):
         i=0
         spe=""
         r=search("^([^\.]*)\..*$",FASTA)
         if r:
            spe=r.group(1)
         else:
            exit("!!! ERROR, FASTA file name: "+FASTA+" is incorrectly written (Must be: ${species_name}.fa) !!!")
         if spe==species:
            FASTA_FILE=INPUT_SCAFF_dir+"/"+FASTA
            print("\t\tStore contigs/scaffolds of FASTA file \""+FASTA_FILE+"\"")
            # Browse FASTA file of current species to get list of scaffolds 
            fasta_file=open(INPUT_SCAFF_dir+"/"+FASTA)
            for sequence in SeqIO.parse(fasta_file,"fasta"):
               scaff=sequence.id
               seq=sequence.seq
               dict_SCAFF[scaff]=seq
            fasta_file.close()
            break

      for scaff in dict_newSCAFF:
         output_scaff.write("> "+scaff+"\n")
         listCTG=dict_newSCAFF[scaff]
         seq=""
         ID=""
         gap_size=0
         posSTART=1
         posEND=0
         bool_default=False
         for ctg in listCTG:
            # Write "N" gap between the 2 contigs scaffold
            if seq:
               if (ID,ctg.id) in dict_distCTG:
                  if dict_distCTG[(ID,ctg.id)]=="?":
                     bool_default=True
                     gap_size=default_gap_size
                  else:
                     gap_size=int(float(dict_distCTG[(ID,ctg.id)]))
                     # print gap_size
                     if gap_size<0:
                        if verbose>1:
                           print("\t\tNEGATIVE distance between contigs "+ID+" and "+ctg.id)
               elif (ctg.id,ID) in dict_distCTG:
                  if dict_distCTG[(ctg.id,ID)]=="?":
                     bool_default=True
                     gap_size=default_gap_size
                  else:
                     gap_size=int(float(dict_distCTG[(ctg.id,ID)]))
                     # print gap_size
                     if gap_size<0:
                        if verbose>1:
                           print("\t\tNEGATIVE distance between contigs "+ctg.id+" and "+ID)
               else:
                  exit("ERROR: CTG adjacency ("+ID+"-"+ctg.id+") is not present in DeCoSTAR predicted adjacencies file: "+AGP_dir+"/"+AGP_file)

               nb_N=0
               if gap_size<0:
                  nb_N=default_gap_size
               else:
                  nb_N=gap_size
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

      # Write remaining contigs (not scaffolded) in the FASTA and AGP files
      for scaff in dict_SCAFF:
         output_scaff.write("> "+scaff+"\n")
         output_scaff.write(str(dict_SCAFF[scaff])+"\n")
      output_scaff.close()

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time-start_time))
