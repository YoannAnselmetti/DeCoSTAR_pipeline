#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###
###   Goal:
###      Transform GFF files in GENE files and get statistics graphics on exons.
###
###   INPUT:
###      1 - INPUT directory containing sorted GFF files
###         (data/GFF_to_GENE_files/sort_GFF)
###      2 - OUTPUT directory where GENE file will be stored
###         (data/GFF_to_GENE_files/GENE_file)
###      3 - OUTPUT directory where graphics on exons will be stored
###         (data/GFF_to_GENE_files/GRAPH_GFF)
###
###   OUTPUT:
###      1/ OUTPUT directory containing 1 file/species containing gene annotation
###      2/ OUTPUT directory containing 1 directory/species containing stats graphics on exon
###
###   Name: from_GFF_to_GENE.py       Author: Yoann Anselmetti
###   Creation date: 2015/07/20       Last modification: 2020/11/05
###


from sys import argv, exit
from re import search
from os import close, path, makedirs, listdir, mkdir
from datetime import datetime
from shutil import rmtree
import errno

import itertools
import matplotlib
matplotlib.use('Agg')
 
import numpy as np 
from matplotlib import pyplot as plt



def hist_Exon_size(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=35,range=(0,3500))
   plt.title("Distribution of Exon size")
   plt.xlabel("Exon size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0, xmax=3500)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Exon_size.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_Intron_size(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=35,range=(0,35000))
   plt.title("Distribution of Intron size")
   plt.xlabel("Intron size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0, xmax=35000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Intron_size.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_Intron_size_zoom(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=50,range=(0,1000))
   plt.title("Distribution of Intron size")
   plt.xlabel("Intron size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0, xmax=1000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Intron_size_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

#def hist_Intron_size_zoom_zoom(list_to_plot,OUTPUT,species):
#   plt.hist(list_to_plot,bins=10,range=(0,100))
#   plt.title("Distribution of Intron size")
#   plt.xlabel("Intron size")
#   plt.ylabel("Frequency")
#   plt.xlim(xmin=0, xmax=100)
#   plt.ylim(ymin=0)
#   fig_name=OUTPUT+"/"+species+"_distrib_Intron_size_zoom_zoom.pdf"
#   plt.tight_layout()
#   plt.savefig(fig_name,format='pdf')
#   plt.cla()

def hist_nbExon_per_gene(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=35,range=(0,35))
   plt.title("Distribution of Exon number per gene")
   plt.xlabel("Exon number per gene")
   plt.ylabel("Frequency")
   plt.xlim(xmin=1, xmax=35)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_nbExon_per_gene.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def graph_Exon_size(dict_to_plot,OUTPUT,species):
   N=len(list(dict_to_plot.keys()))
   index = np.arange(N)
   width = .75
   list_mean=list()
   list_SD=list()

   list_keys=list(dict_to_plot.keys())
   list_keys.sort()
   for spe in list_keys:
#      print(str(dict_to_plot[species][0])+" "+str(dict_to_plot[species][1]))
      list_mean.append(dict_to_plot[spe][0])
      list_SD.append(dict_to_plot[spe][1])
   plt.bar(index, list_mean, width, color="green", yerr=list_SD)
   plt.title("Mean Exon size per species")
   plt.ylabel("Mean Exon size", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=1000)
   fig_name=OUTPUT+"/"+species+"_mean_Exon_size_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_Intron_size(dict_to_plot,OUTPUT,species):
   N=len(list(dict_to_plot.keys()))
   index = np.arange(N)
   width = .75
   list_mean=list()
   list_SD=list()

   list_keys=list(dict_to_plot.keys())
   list_keys.sort()
   for spe in list_keys:
      list_mean.append(dict_to_plot[spe][0])
      list_SD.append(dict_to_plot[spe][1])
   plt.bar(index, list_mean, width, color="green", yerr=list_SD)
   plt.title("Mean Intron size per species")
   plt.ylabel("Mean Intron size", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=7500)
   fig_name=OUTPUT+"/"+species+"_mean_Intron_size_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_Intron_size_zoom(dict_to_plot,OUTPUT,species):
   N=len(list(dict_to_plot.keys()))
   index = np.arange(N)
   width = .75
   list_mean=list()
   list_SD=list()

   list_keys=list(dict_to_plot.keys())
   list_keys.sort()
   for spe in list_keys:
      list_mean.append(dict_to_plot[spe][0])
      list_SD.append(dict_to_plot[spe][1])
   plt.bar(index, list_mean, width, color="green", yerr=list_SD)
   plt.title("Mean Intron size per species")
   plt.ylabel("Mean Intron size", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=1500)
   fig_name=OUTPUT+"/"+species+"_mean_Intron_size_per_species_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_nbExon_per_gene(dict_to_plot,OUTPUT,species):
   N=len(list(dict_to_plot.keys()))
   index = np.arange(N)
   width = .75
   list_mean=list()
   list_SD=list()

   list_keys=list(dict_to_plot.keys())
   list_keys.sort()
   for spe in list_keys:
      list_mean.append(dict_to_plot[spe][0])
      list_SD.append(dict_to_plot[spe][1])
   plt.bar(index, list_mean, width, color="green", yerr=list_SD)
   plt.title("Mean Exon number per gene per species")
   plt.ylabel("Mean Exon number per gene", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin =0)
   fig_name=OUTPUT+"/"+species+"_mean_nbExon_per_gene_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def mean(table):
    return sum(table, 0.0) / len(table)

def variance(table):
    m=mean(table)
    return mean([(x-m)**2 for x in table])

def SD(table):
   return variance(table)**0.5


def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise



if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   INPUT_dir=argv[1]
   OUTPUT_dir_GENE=argv[2]
   OUTPUT_dir_GRAPH=argv[3]


   # To be sure than directory have no "/" to the end of the path
   r=search('^(.*)/$', INPUT_dir)
   if r:
      INPUT_dir=r.group(1)
   # To be sure than directory have no "/" to the end of the path
   r=search('^(.*)/$', OUTPUT_dir_GENE)
   if r:
      OUTPUT_dir_GENE=r.group(1)
   # To be sure than directory have no "/" to the end of the path
   r=search('^(.*)/$', OUTPUT_dir_GRAPH)
   if r:
      OUTPUT_dir_GRAPH=r.group(1)

   # If OUTPUT_dir_GENE exists: Remove directory and recreate it
   if path.exists(OUTPUT_dir_GENE):
      rmtree(OUTPUT_dir_GENE)
   mkdir_p(OUTPUT_dir_GENE)
   # If OUTPUT_dir_GRAPH exists: Remove directory and recreate it
   if path.exists(OUTPUT_dir_GRAPH):
      rmtree(OUTPUT_dir_GRAPH)
   mkdir_p(OUTPUT_dir_GRAPH)

   # Get list of files contained in input directory
   list_files=listdir(INPUT_dir)

   nb_gene_tot=0
   nb_exon_tot=0
   dict_spe_exon_size={}
   dict_spe_intron_size={}
   dict_spe_nbExon_per_gene={}
   list_exon_size=list()
   list_intron_size=list()
   list_nbExon_per_gene=list()

   gID=""
   # Browse GFF files contained in INPUT directory
   for gff_file in sorted(list_files):
      r_spe=search('^(.*)_sorted\.gff3$', gff_file)
      if r_spe:
         species=r_spe.group(1)
         print("\nFor species "+species+":")
         contig_name=""
         start_pos=0
         end_pos=0
         gene_ID=""
         nb_exon=0
         nb_gene_tot_spe=0
         nb_exon_tot_spe=0
         list_exon_size_spe=list()
         list_intron_size_spe=list()
         list_nbExon_per_gene_spe=list()

         # Get Exon position of a Gene
         list_Exon_pos=list()

         # Path of output files
         OUTPUT_Gene_position=OUTPUT_dir_GENE+"/"+species+".txt"
         # Open output files to write in it
         output_Gpos=open(OUTPUT_Gene_position,'w')
         output_Gpos.write("#species\tctg/scaff\tgene\torientation_gene\tstart_gene\tend_gene\t#exons\texons_position\n")

         # Browse GFF file
         INPUT_GFF=INPUT_dir+"/"+gff_file
         input_file=open(INPUT_GFF,'r')

         #######
         ### USE BIOPYTHON
         #######
         for line in input_file:
            r=search('^([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\n\t]*)\n$', line)
            if r:
               seqname=r.group(1)
               source=r.group(2)
               feature=r.group(3)
               start=r.group(4)
               end=r.group(5)
               score=r.group(6)
               strand=r.group(7)
               frame=r.group(8)
               gID=r.group(9)

               # Get Exon size and add it at list of Intron size
               exon_size=int(end)-int(start)+1
               list_exon_size_spe.append(exon_size)

               # If it's not the first gene of the species
               if gene_ID!="":
                  # If we change of gene_ID
                  if gene_ID!=gID:
                     # Print END position of the last gID
                     output_Gpos.write(str(end_pos)+"\t"+str(nb_exon))
                     first_exon=True
                     for exon in list_Exon_pos:
                        if first_exon:
                           output_Gpos.write("\t"+exon)
                           first_exon=False
                        else:
                           output_Gpos.write(":"+exon)
                     output_Gpos.write("\n")
                     del list_Exon_pos[:]
                     # Get nb_exon in gID and store it in list of Exon number per gID
                     list_nbExon_per_gene_spe.append(nb_exon)
                     # Print START position of the current gID in OUTPUT GENE file
                     output_Gpos.write(species+"\t"+seqname+"\t"+gID+"\t"+strand+"\t"+start+"\t")
                     nb_exon=0
                     nb_gene_tot_spe+=1
                  # If we don't change of gene_ID
                  else:
                     # Get Intron size and add it at list of Intron size
                     intron_size=int(start)-end_pos-1
                     list_intron_size_spe.append(intron_size)
#                     output_intron.write(str(intron_size)+"+")
               # If it's the first gID of the species
               else:
                  # Print START position of the current gID
                  output_Gpos.write(species+"\t"+seqname+"\t"+gID+"\t"+strand+"\t"+start+"\t")
                  nb_gene_tot_spe=1
               nb_exon+=1
               nb_exon_tot_spe+=1
               list_Exon_pos.append(start+"-"+end)


               # Update variables with data from the line
               contig_name=seqname
               start_pos=int(start)
               end_pos=int(end)
               gene_ID=gID

            else:
               exit("ERROR in file "+species+" wrong line written "+line)

         # Write END position of the last gene
         output_Gpos.write(str(end_pos)+"\t"+str(nb_exon))
         first_exon=True
         for exon in list_Exon_pos:
            if first_exon:
               output_Gpos.write("\t"+exon)
               first_exon=False
            else:
               output_Gpos.write(":"+exon)
         output_Gpos.write("\n")
         del list_Exon_pos[:]
         output_Gpos.close()

         # Get Exon number for the last gene
         list_nbExon_per_gene_spe.append(nb_exon)
         # Fill dictionnary for the different statistics
         dict_spe_exon_size[species]=[mean(list_exon_size_spe),SD(list_exon_size_spe)]
         dict_spe_intron_size[species]=[mean(list_intron_size_spe),SD(list_intron_size_spe)]
         dict_spe_nbExon_per_gene[species]=[mean(list_nbExon_per_gene_spe),SD(list_nbExon_per_gene_spe)]

         #Store species list() and variable in all_species list() and variable
         list_exon_size+=list_exon_size_spe
         list_intron_size+=list_intron_size_spe
         list_nbExon_per_gene+=list_nbExon_per_gene_spe
         nb_exon_tot+=nb_exon_tot_spe
         nb_gene_tot+=nb_gene_tot_spe

         # Draw histogram for different statistics on Exon and Intron
         OUTPUT_dir_GRAPH_Spe=OUTPUT_dir_GRAPH+"/"+species
         mkdir_p(OUTPUT_dir_GRAPH_Spe)
         hist_Exon_size(list_exon_size_spe,OUTPUT_dir_GRAPH_Spe,species)
         hist_Intron_size(list_intron_size_spe,OUTPUT_dir_GRAPH_Spe,species)
         hist_Intron_size_zoom(list_intron_size_spe,OUTPUT_dir_GRAPH_Spe,species)
#         hist_Intron_size_zoom_zoom(list_intron_size_spe,OUTPUT_dir_GRAPH_Spe,species)
         hist_nbExon_per_gene(list_nbExon_per_gene_spe,OUTPUT_dir_GRAPH_Spe,species)

         # Print summary available on current species
         print("\t- Gene number:"+str(nb_gene_tot_spe))
         print("\t- Exon number:"+str(nb_exon_tot_spe))
         print("\t- Average Exon size:"+str(mean(list_exon_size_spe)))
         print("\t- SD Exon size:"+str(SD(list_exon_size_spe)))
         print("\t- Average interExon size:"+str(mean(list_intron_size_spe)))
         print("\t- SD interExon size:"+str(SD(list_intron_size_spe)))
         print("\t- Average Exon number per Gene:"+str(mean(list_nbExon_per_gene_spe)))
         print("\t- SD Exon number per Gene:"+str(SD(list_nbExon_per_gene_spe)))

      else:
         exit('Name of file '+gff_file+' is bad written!!!')

   # Draw histogram for different statistics on Exon and Intron
   OUTPUT_dir_GRAPH_Spe=OUTPUT_dir_GRAPH+"/ALL_species"
   mkdir_p(OUTPUT_dir_GRAPH_Spe)
   hist_Exon_size(list_exon_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   hist_Intron_size(list_intron_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   hist_Intron_size_zoom(list_intron_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
#   hist_Intron_size_zoom_zoom(list_intron_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   hist_nbExon_per_gene(list_nbExon_per_gene,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_Exon_size(dict_spe_exon_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_Intron_size(dict_spe_intron_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_Intron_size_zoom(dict_spe_intron_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_nbExon_per_gene(dict_spe_nbExon_per_gene,OUTPUT_dir_GRAPH_Spe,"ALL_species")

   # Print summary available on ALL species
   print("\nFor ALL species:")
   print("\t- Gene number:"+str(nb_gene_tot))
   print("\t- Exon number:"+str(nb_exon_tot))
   print("\t- Average Exon size:"+str(mean(list_exon_size)))
   print("\t- SD Exon size:"+str(SD(list_exon_size)))
   print("\t- Average interExon size:"+str(mean(list_intron_size)))
   print("\t- SD interExon size:"+str(SD(list_intron_size)))
   print("\t- Average Exon number per Gene:"+str(mean(list_nbExon_per_gene)))
   print("\t- SD Exon number per Gene:"+str(SD(list_nbExon_per_gene)))

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time-start_time))
