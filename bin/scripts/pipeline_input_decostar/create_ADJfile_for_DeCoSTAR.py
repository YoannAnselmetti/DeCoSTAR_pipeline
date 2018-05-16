#! /usr/bin/env python2
# -*- coding: utf-8 -*-
###
###   Goal:                                                                    
###      Create Adjacency files for DeCo* from GENE file and scaffolding       
###      adjacencies file                                                      
###                                                                            
###   INPUT:                                                                   
###      1- INPUT file with gene position for all species                      
###         (data/data_DeCoSTAR/GENE_file)                           
###      2- INPUT scaffolding adjacencies file                                           
###         (data/data_DeCoSTAR/scaff_BESST_DeCoSTAR)   
###      3- OUTPUT file path for Adjacency file                                
###         (data/data_DeCoSTAR/decostar/adjacencies.txt)
###      4- OUTPUT file path for Adjacency file without scaffolding                               
###         (data/data_DeCoSTAR/decostar/adjacencies-scaff.txt)
###      5- Character separator between species name and gene ID               
###         (@)                                                                
###                                                                            
###   OUTPUT:   (RUN in <1sec)                                                 
###      - INPUT files for ARt-DeCo_seq (Gene & Adjacency files)               
###                                                                            
###   Name: create_ADJfile_for_DeCoSTAR.py     Author: Yoann Anselmetti     
###   Creation date: 2016/04/01                Last modification: 2017/11/09
###

from sys import argv, stdout
from datetime import datetime
from re import search
from os import close, path, makedirs
import errno


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else:
         raise


if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   GENE_file=argv[1]
   SCAFF_file=argv[2]
   OUTPUT_Adj=argv[3]
   OUTPUT_Adj_no_scaff=argv[4]
   separator=argv[5]

   
   # Create OUTPUT_DIR if not existing
   OUTPUT_DIR=path.dirname(path.realpath(OUTPUT_Adj))
   mkdir_p(OUTPUT_DIR)

   gene_file=open(GENE_file,'r')
   scaff_file=open(SCAFF_file,'r')
   output_adj=open(OUTPUT_Adj,'w')
   output_adj_no_scaff=open(OUTPUT_Adj_no_scaff,'w')

   species=""
   contig=""
   gene_stored=""
   ori_stored=""
   stop_stored=""
   gene2=""

   for gene in gene_file:
      r=search('^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$', gene)
      if r:
         name_species=r.group(1)
         contig_ID=r.group(2)
         gf_ID=r.group(3)
         gene_ID=r.group(4)
         orientation=r.group(5)
         start=r.group(6)
         stop=r.group(7)
         exon_nb=r.group(8)
         exon_pos=r.group(9)

         if name_species!="#species":
            if name_species==species and contig_ID==contig:
               dist=int(start)-int(stop_stored)
               # output_adj.write(name_species+separator+gene_stored+"\t"+name_species+separator+gene_ID+"\t"+ori_stored+"\t"+orientation+"\t"+str(dist)+"\n")
               output_adj.write(name_species+separator+gene_stored+"\t"+name_species+separator+gene_ID+"\t"+ori_stored+"\t"+orientation+"\n")
               output_adj_no_scaff.write(name_species+separator+gene_stored+"\t"+name_species+separator+gene_ID+"\t"+ori_stored+"\t"+orientation+"\n")
               gene_stored=gene_ID
               ori_stored=orientation
               stop_stored=stop
            else:
               species=name_species
               contig=contig_ID
               gene_stored=gene_ID
               ori_stored=orientation
               stop_stored=stop
      else:
           exit("!!! ERROR line: "+gene+" is bad written !!!")
   gene_file.close()
   output_adj_no_scaff.close()


   for scaff_adj in scaff_file:
      r=search('^([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\t]*)[\t]*([^\n]*)\n$', scaff_adj)
      if r:
         species=r.group(1)
         ctg1=r.group(2)
         ctg2=r.group(3)
         ori_CTG1=r.group(4)
         ori_CTG2=r.group(5)
         gap=r.group(6)
         gf1=r.group(7)
         gf2=r.group(8)
         g1=r.group(9)
         g2=r.group(10)
         ori1=r.group(11)
         ori2=r.group(12)
         dist=r.group(13)
         link_score=r.group(14)
         disp_score=r.group(15)
         links_nb=r.group(16)

         if species!="#species":
            score=str((float(link_score)+float(disp_score))/2.0)
            # output_adj.write(species+separator+g1+"\t"+species+separator+g2+"\t"+ori1+"\t"+ori2+"\t"+dist+"\t"+score+"\n")
            output_adj.write(species+separator+g1+"\t"+species+separator+g2+"\t"+ori1+"\t"+ori2+"\t"+score+"\n")
      else:
           exit("!!! ERROR line: "+scaff_adj+" is bad written !!!")
   scaff_file.close()
   output_adj.close()

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
