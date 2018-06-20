#! /usr/bin/env python2
# -*- coding: utf-8 -*-
###                                                                         
###   Goal:                                                                 
###      Create scaffs file with gene located to the extremities of scaffs
###                                                                         
###   INPUT:                                                                
###      1- Annotation gene file                                            
###         (data/data_DeCoSTAR/GENE_file)                        
###      2- Genome Assembly FASTA files directory                           
###         (data/INPUT_DATA/FASTA/SCAFF)                 
###      3- Contigs OUTPUT file giving info on genes to extremities         
###         (data/data_DeCoSTAR/SCAFF_file)                         
###                                                                         
###   OUTPUT:                                                               
###      - Contig file with genes located to extremities                    
###                                                                         
###   Name: create_SCAFF_file.py              Author: Yoann Anselmetti     
###   Creation date: 2016/02/15             Last modification: 2018/04/24
###                                                                         

from sys import argv
from re import search
from os import close, path, listdir, makedirs
from datetime import datetime
from Bio import SeqIO
import errno


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else:
         raise


################
###   MAIN 
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   annotG_file=argv[1]
   SCAFF_dir=argv[2]
   SCAFF_file=argv[3]

   verbose=True

   OUTPUT_DIR=path.dirname(path.realpath(SCAFF_file))
   # Create OUTPUT_DIR if not existing
   mkdir_p(OUTPUT_DIR)

   # To be sure than directory have no "/" to the end of the path
   SCAFF_dir=path.normpath(SCAFF_dir)


#################################################################
### Store list of scaff size/species in dict_spe_listSCAFFsize ###
#################################################################
   dict_spe_SCAFF=dict()
   FASTA_list=listdir(SCAFF_dir)
   print "1/ Storing SCAFF size of species:"
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
      fasta_file=open(SCAFF_dir+"/"+FASTA)
      for sequence in SeqIO.parse(fasta_file,"fasta"):
         scaff=sequence.id
         size=len(sequence.seq)

         if scaff in dict_spe_SCAFF[species]:
               exit("!!! ERROR, there are two scaffolds that have the same ID ("+scaff+") !!! => Did you uncompress scaffolds assemblies in the directory: \""+SCAFF_dir+"\"?")

         dict_spe_SCAFF[species][scaff]=size
         i+=1

      fasta_file.close()
      if i==0:
         exit("ERROR, there is no scaff for this species "+species+"!!! => Check if the genome assemblies are uncompressed.")


   print "\n2/ Write CTG file \""+SCAFF_file+"\":"
####################################################################
### BROWSE ANNOTATION GENE FILE TO WRITE GENE TO SCAFF EXTREMITIES ###
####################################################################
   input_file=open(annotG_file,'r')
   output_file=open(SCAFF_file,'w')
   output_file.write("#species\tctg\tctg_size\tctg_gene_nb\t5'_gene_family\t5'_gene\torientation_5'_gene\tstart_5'_gene\t3'_gene_family\t3'_gene\torientation_3'_gene\tend_3'_gene\n")
   gf1=""
   gene1=""
   ori1=""
   start1=""
   str_spe=""
   str_ctg=""
   str_gf=""
   str_gene=""
   str_ori=""
   str_stop=""
   gene_nb=0
   for line in input_file:
      r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
      if r:
         cur_spe=r.group(1)
         cur_ctg=r.group(2)
         cur_gf=r.group(3)
         cur_gene=r.group(4)
         cur_ori=r.group(5)
         cur_start=r.group(6)
         cur_stop=r.group(7)
         exon_nb=r.group(8)
         exon_pos=r.group(9)

         if cur_spe!="#species":
            if str_spe=="":
               # print cur_spe
               gf1=cur_gf
               gene1=cur_gene
               ori1=cur_ori
               start1=cur_start

            else:
               if (str_spe!=cur_spe) or (str_ctg!=cur_ctg):
                  if str_ctg in dict_spe_SCAFF[str_spe]:
                     output_file.write(str_spe+"\t"+str_ctg+"\t"+str(dict_spe_SCAFF[str_spe][str_ctg])+"\t"+str(gene_nb)+"\t"+gf1+"\t"+gene1+"\t"+ori1+"\t"+start1+"\t"+str_gf+"\t"+str_gene+"\t"+str_ori+"\t"+str_stop+"\n")
                  # For scaff that are not present in SCAFF directory (Genome assemblies)
                  else:
                     if verbose:
                        print "\tContig "+str_ctg+" is not present in FASTA file assembly of species "+str_spe
                     output_file.write(str_spe+"\t"+str_ctg+"\t?\t"+str(gene_nb)+"\t"+gf1+"\t"+gene1+"\t"+ori1+"\t"+start1+"\t"+str_gf+"\t"+str_gene+"\t"+str_ori+"\t"+str_stop+"\n")

                  gene_nb=0
                  gf1=cur_gf
                  gene1=cur_gene
                  ori1=cur_ori
                  start1=cur_start

            str_spe=cur_spe
            str_ctg=cur_ctg
            str_gf=cur_gf
            str_gene=cur_gene
            str_ori=cur_ori
            str_stop=cur_stop
            gene_nb+=1

   input_file.close()

   # Write last scaff in SCAFF_file
   if str_ctg in dict_spe_SCAFF[str_spe]:
      output_file.write(str_spe+"\t"+str_ctg+"\t"+str(dict_spe_SCAFF[str_spe][str_ctg])+"\t"+str(gene_nb)+"\t"+gf1+"\t"+gene1+"\t"+ori1+"\t"+start1+"\t"+str_gf+"\t"+str_gene+"\t"+str_ori+"\t"+str_stop+"\n")
   # For scaff that are not present in SCAFF directory (Genome assemblies)
   else:
      output_file.write(str_spe+"\t"+str_ctg+"\t?\t"+str(gene_nb)+"\t"+gf1+"\t"+gene1+"\t"+ori1+"\t"+start1+"\t"+str_gf+"\t"+str_gene+"\t"+str_ori+"\t"+str_stop+"\n")
   output_file.close()

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
