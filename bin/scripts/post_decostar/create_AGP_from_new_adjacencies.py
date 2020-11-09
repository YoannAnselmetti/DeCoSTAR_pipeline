#! /usr/bin/env python3
# -*- coding: utf-8 -*-
###
###   Goal:
###      Produce AGP files summarizing scaffolding results of DeCoSTAR
###
###   INPUT:
###      1- New adjacencies file
###         (27avian_dataset/results/decostar/ADseq+scaff_Boltz_kT0.1/DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Linear_0.1_M_new_adjacencies_with_scaff)
###      2- Genome assemblies directory
###         (27avian_dataset/data/INPUT_DATA/FASTA/SCAFF)
###      3- prefix file 
###         (DeCoSTAR_27avian_ADseq+scaff_Boltz_kT0.1_Lin0.1_M2_)
###      4- AGP files path/prefix 
###         (27avian_dataset/results/AGP/SCAFF/)
###
###   OUTPUT:
###      - AGP files summarizing scaffolding results of DeCoSTAR
###
###   Name: create_AGP_from_new_adjacencies.py    Author: Yoann Anselmetti     
###   Creation date: 2018/06/11                   Last modification: 2020/11/04
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


def mergeCTG(species,ctg1,ctg2,ori_ctg1,ori_ctg2,dict_newSCAFF_ID,dict_newSCAFF,scaffold_ID_with_ctg_order,scaff_id):
   # Get new ID of ctg1 and ctg2 if they have been previously changed
   new_ctg1,new_ctg2="",""
   if species in dict_newSCAFF_ID:
      if ctg1 in dict_newSCAFF_ID[species]:
         new_ctg1=dict_newSCAFF_ID[species][ctg1]
      if ctg2 in dict_newSCAFF_ID[species]:
         new_ctg2=dict_newSCAFF_ID[species][ctg2]

   # If ctg1 has been previously changed
   new_ctg1_order=list()
   new_ctg2_order=list()
   if new_ctg1:
      first_ctg1=dict_newSCAFF[species][new_ctg1][0]
      last_ctg1=dict_newSCAFF[species][new_ctg1][-1]
      # Get the correct order of contigs in new_ctg1
      if ctg1==first_ctg1.id:
         new_ctg1_order=rev_ctg_order(dict_newSCAFF[species][new_ctg1])
      elif ctg1==last_ctg1.id:
         new_ctg1_order=dict_newSCAFF[species][new_ctg1]
      else:
         exit("ERROR, contig \""+ctg1+"\" should be to the start or the end of the new contig \""+new_ctg1+"\" !!!")
      # Remove the previous version of "new_ctg1"
      dict_newSCAFF[species].pop(new_ctg1,None)

   # If ctg2 has been previously changed
   if new_ctg2:        
      first_ctg2=dict_newSCAFF[species][new_ctg2][0]
      last_ctg2=dict_newSCAFF[species][new_ctg2][-1]
      # Get the correct order of contigs in new_27avian_dataset/data/INPUT_DATA/27avian_species.t2
      if ctg2==last_ctg2.id:
         new_ctg2_order=rev_ctg_order(dict_newSCAFF[species][new_ctg2])
      elif ctg2==first_ctg2.id:
         new_ctg2_order=dict_newSCAFF[species][new_ctg2]
      else:
         exit("ERROR, contig \""+ctg2+"\" should be to the start or the end of the new contig \""+new_ctg2+"\" !!!")
      # Remove the previous version of "new_ctg2"
      dict_newSCAFF[species].pop(new_ctg2,None)

   # If ctg1 has NOT been previously changed
   if not new_ctg1_order:
      new_ctg1_order.append(CTG(ctg1,ori_ctg1))
   # If ctg2 has NOT been previously changed
   if not new_ctg2_order:
      new_ctg2_order.append(CTG(ctg2,ori_ctg2))

   # Concatenate the two CTG order to define the new SCAFF
   new_ctg_order=new_ctg1_order+new_ctg2_order


   # Set the ID of the new scaffold
   new_CTG_ID=""
   ### If we set the scaffold ID with the ctg order list
   if scaffold_ID_with_ctg_order:
      for elem in new_ctg_order:
         ctg=elem.id
         ori=elem.ori
         short_ctg=ctg.split(sep)[0]
         if not new_CTG_ID:
            new_CTG_ID+=short_ctg+"("+ori+")"
         else:
            new_CTG_ID+=":"+short_ctg+"("+ori+")"
   ### If we set the scaffold ID with the format "SCAFFXXXXXX"
   else:
      if new_ctg1:
         new_CTG_ID=new_ctg1
      elif new_ctg2:
         new_CTG_ID=new_ctg2
      else:
         scaff_id,new_CTG_ID=get_SCAFF_ID(scaff_id)
   
   # 
   if not species in dict_newSCAFF:
      dict_newSCAFF[species]=dict()
   dict_newSCAFF[species][new_CTG_ID]=new_ctg_order

   if not species in dict_newSCAFF_ID:
      dict_newSCAFF_ID[species]=dict()
   for elem in new_ctg_order:
      ctg=elem.id
      dict_newSCAFF_ID[species][ctg]=new_CTG_ID

   return scaff_id



################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   CTG=namedtuple("CTG",["id","ori"])

   # Recovery of input parameters
   newAdj_file=argv[1]
   FASTA_dir=argv[2]
   file_prefix=argv[3]
   OUTPUT_AGP=argv[4]

   scaffold_ID_with_ctg_order=False
   write_unscaffolded_ctg=True
   default_gap_size=100
   verbose=1
   sep="_size"

   # Create OUTPUT_DIR if not existing
   mkdir_p(OUTPUT_AGP)


   print("1/ Merge scaffolds linked by linearized new adjacencies")
   dict_newSCAFF,dict_newSCAFF_ID,dict_distCTG=dict(),dict(),dict()
   input_file=open(newAdj_file,"r")
   scaff_id=0
   stored_species=""
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
            if stored_species!=species:
               scaff_id=0
               stored_species=species
            dict_distCTG[(ctg1,ctg2)]=dist_ctg
            scaff_id=mergeCTG(species,ctg1,ctg2,ori_ctg1,ori_ctg2,dict_newSCAFF_ID,dict_newSCAFF,scaffold_ID_with_ctg_order,scaff_id)
   input_file.close()
   dict_newSCAFF_ID.clear()


   # Print the new CTG ID after merging initial contigs/scaffoldsin new scaffolds 
   if verbose>1:
      if verbose>2:
         print("\n1bis/ Print association between new scaffolds ID and old scaffolds/contigs ID / species:")
      else:
         print("\n1bis/ Print new scaffolds ID / species:")
      # Print new linked scaffolds / species 
      for species in sorted(dict_newSCAFF):
         print("\n"+species+":")
         for ctg in dict_newSCAFF[species]:
            print("\t"+ctg)
            if verbose>2:
               for elem in dict_newSCAFF[species][ctg]:
                  print("\t\t", end=' ')
                  print(elem)

   # Print distance between CTG pairs linked by DeCoSTAR
   if verbose>2:
      print("\n\t1ter/ Print distance between CTG pairs linked by DeCoSTAR:")
      for distCTG in sorted(dict_distCTG):
         print(distCTG,"\t",dict_distCTG[distCTG])


   print("\n2/ Write AGP files after scaffolding with linearized new adjacencies predicted by DeCoSTAR")
   for species in sorted(dict_newSCAFF):
      if verbose>0:
         print("\t"+species)
      output_agp=open(OUTPUT_AGP+"/"+file_prefix+species+".agp","w")

      dict_CTG=dict()
      FASTA_FILE=""
      for FASTA in sorted(listdir(FASTA_dir)):
         i=0
         spe=""
         r=search("^([^\.]*)\..*$",FASTA)
         if r:
            spe=r.group(1)
         else:
            exit("!!! ERROR, FASTA file name: "+FASTA+" is incorrectly written (Must be: ${species_name}.fa) !!!")
         if spe==species:
            FASTA_FILE=FASTA_dir+"/"+FASTA
            # Browse FASTA file of current species to get list of scaffolds 
            fasta_file=open(FASTA_dir+"/"+FASTA)
            for sequence in SeqIO.parse(fasta_file,"fasta"):
               ctg=sequence.id
               size=len(sequence.seq)
               dict_CTG[ctg]=size
            fasta_file.close()
            break

      for scaff in sorted(dict_newSCAFF[species]):
         listCTG=dict_newSCAFF[species][scaff]
         ID=""
         gap_size=0
         posSTART=1
         posEND=0
         bool_default=False
         stored_CTG="" 
         for ctg in listCTG:
            ### Allow to discard duplicate CTG (corresponding to a circularizing adjacency predicted by DeCoSTAR)
            if ctg.id!=stored_CTG:
               stored_CTG=ctg.id
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
                     exit("ERROR: CTG adjacency ("+ID+"-"+ctg.id+") is not present in DeCoSTAR predicted adjacencies file: "+newAdj_file)

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
      if write_unscaffolded_ctg:
         for ctg in sorted(dict_CTG):
            ctg_size=dict_CTG[ctg]
            output_agp.write(ctg+"\t1\t"+str(ctg_size)+"\t.\tW\t"+ctg+"\t1\t"+str(ctg_size)+"\t+\n")
      output_agp.close()


   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time-start_time))
