#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###
###   Goal:
###      - Detect included genes and remove them
###      - Produce an unique GENE file for all species considered
###      - Produce statistics graphs on genes
###
###   INPUT:
###      1- INPUT directory containing sorted GENE files
###         (data/GFF_to_GENE_files/sorted_GENE)
###         (data/GFF_to_GENE_files/filtered_GENE)
###      2- OUTPUT directory path where results will be stored
###         (data/GFF_to_GENE_files/without_filter)
###         (data/GFF_to_GENE_files/with_filter)
###   OUTPUT:
###      - OUTPUT directory containing 1 directory for each species composed of statistics graphs on genes
###      - OUTPUT files:
###         + OUTPUT annotation gene file => Without included genes and modification of exon position, localized on overlapping zones
###            (ALL_species_GENE_file)   
###         + OUTPUT inclusion gene file
###            (ALL_species_Inclusion_file)
###         + ALL_species_Overlap_file
###            (OUTPUT overlap gene file)
###
###   Name: detect_includedGenes.py     Author: Yoann Anselmetti
###   Creation date: 2015/07/20         Last modification: 2020/11/05
###

from sys import argv, exit
from re import search, match
from os import close, path, makedirs, listdir, mkdir
from datetime import datetime
import subprocess
import errno

import itertools
import matplotlib
matplotlib.use('Agg')
 
import numpy as np 
from matplotlib import pyplot as plt

def hist_nbGene_contig(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=80,range=(0,4000))
   plt.title("Distribution of Gene number per contig.")
   plt.xlabel("Gene number")
   plt.ylabel("Contig number")
   plt.xlim(xmin=0,xmax=4000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_GeneNb_per_contig.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_nbGene_contig_zoom(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=50,range=(0,250))
   plt.title("Distribution of Gene number per contig.")
   plt.xlabel("Gene number")
   plt.ylabel("Contig number")
   plt.xlim(xmin=0,xmax=250)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_GeneNb_per_contig_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_nbGene_contig_zoom_zoom(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=50,range=(0,50))
   plt.title("Distribution of Gene number per contig.")
   plt.xlabel("Gene number")
   plt.ylabel("Contig number")
   plt.xlim(xmin=0,xmax=50)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_GeneNb_per_contig_zoom_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def hist_Gene_size(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=150,range=(0,150000))
   plt.title("Distribution of Gene size.")
   plt.xlabel("Gene size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=150000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Gene_size.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_Gene_size_zoom(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=200,range=(0,20000))
   plt.title("Distribution of Gene size.")
   plt.xlabel("Gene size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=20000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Gene_size_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def hist_Gene_size_Exon(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=200,range=(0,20000))
   plt.title("Distribution of Gene size.")
   plt.xlabel("Gene size (only exon)")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=20000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Gene_size_Exon.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def hist_interGene_size(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=50,range=(0,500000))
   plt.title("Distribution of interGene size.")
   plt.xlabel("interGene size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=500000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_interGene_size.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_interGene_size_zoom(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=100,range=(0,10000))
   plt.title("Distribution of interGene size.")
   plt.xlabel("interGene size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=10000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_interGene_size_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def hist_overlap_size(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=400,range=(0,800000))
   plt.title("Distribution of Overlap size.")
   plt.xlabel("Overlap size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=80000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Overlap_size.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_overlap_size_zoom(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=40,range=(0,2000))
   plt.title("Distribution of Overlap size.")
   plt.xlabel("Overlap size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=2000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Overlap_size_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def hist_Exon_overlap_ratio_min(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=20,range=(0,100))
   plt.title("Distribution of Min Overlap ratio on Exon.")
   plt.xlabel("Min Overlap ratio")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=100)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_min_Overlap_ratio_on_Exon.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_Exon_overlap_ratio_max(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=20,range=(0,100))
   plt.title("Distribution of Max Overlap ratio on Exon.")
   plt.xlabel("Max Overlap ratio")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=100)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_max_Overlap_ratio_on_Exon.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_Exon_overlap_ratio(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=20,range=(0,100))
   plt.title("Distribution of Overlap ratio on Exon.")
   plt.xlabel("Overlap ratio")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=100)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Overlap_ratio_on_Exon.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def hist_inclusion_size(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=400,range=(0,800000))
   plt.title("Distribution of Inclusion size.")
   plt.xlabel("Inclusion size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=80000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Inclusion_size.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_inclusion_size_zoom(list_to_plot,OUTPUT,species):
   plt.hist(list_to_plot,bins=40,range=(0,2000))
   plt.title("Distribution of Inclusion size.")
   plt.xlabel("Inclusion size")
   plt.ylabel("Frequency")
   plt.xlim(xmin=0,xmax=2000)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_Inclusion_size_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def hist_nbOverlap_gene(list_to_plot,OUTPUT,species):
   index = np.arange(5)
   width = 0.5
   plt.hist(list_to_plot,bins=4,range=(1,5))
   plt.title("Distribution of Gene number per Overlap number.")
   plt.xlabel("Overlap number")
   plt.ylabel("Gene number")
   plt.xlim(xmin=1,xmax=5)
   plt.xticks(index+width,[i for i in index])
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_OverlapNb_per_gene.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def hist_nbInclusion_gene(list_to_plot,OUTPUT,species):
   index = np.arange(11)
   width = 0.5
   plt.hist(list_to_plot,bins=10,range=(1,11))
   plt.title("Distribution of Gene number per Inclusion number.")
   plt.xlabel("Inclusion number")
   plt.ylabel("Gene number")
   plt.xlim(xmin=1,xmax=11)
   plt.xticks(index+width,[i for i in index])
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_distrib_InclusionNb_per_gene.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()



def graph_ctgNb(dict_to_plot,OUTPUT,species):
   N=len(list(dict_to_plot.keys()))
   index = np.arange(N)
   width = .75
   list_ctgNb=list()

   list_keys=list(dict_to_plot.keys())
   list_keys.sort()
   for spe in list_keys:
      list_ctgNb.append(dict_to_plot[spe])

   plt.bar(index, list_ctgNb, width, color="green")
   plt.title("Contig number (with Gene) per species")
   plt.ylabel("Contig number (with Gene)", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_contig_number_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def graph_nbGene_contig(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Gene number per contig per species")
   plt.ylabel("Mean Gene number per contig", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_mean_Gene_number_per_contig_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_nbGene_contig_zoom(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Gene number per contig per species")
   plt.ylabel("Mean Gene number per contig", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=250)
   fig_name=OUTPUT+"/"+species+"_mean_Gene_number_per_contig_per_species_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_nbGene_contig_zoom_zoom(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Gene number per contig per species")
   plt.ylabel("Mean Gene number per contig", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=50)
   fig_name=OUTPUT+"/"+species+"_mean_Gene_number_per_contig_per_species_zoom_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_Gene_size(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Gene size per species")
   plt.ylabel("Mean Gene size", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=20000)
   fig_name=OUTPUT+"/"+species+"_mean_Gene_size_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_Gene_size_zoom(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Gene size per species")
   plt.ylabel("Mean Gene size", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=7500)
   fig_name=OUTPUT+"/"+species+"_mean_Gene_size_per_species_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_Gene_size_zoom_zoom(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Gene size per species")
   plt.ylabel("Mean Gene size", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=5000)
   fig_name=OUTPUT+"/"+species+"_mean_Gene_size_per_species_zoom_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_Gene_size_Exon(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Gene size (Exon only) per species")
   plt.ylabel("Mean Gene size (Exon only)", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_mean_Gene_size_Exon_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_interGene_size(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean interGene size per species")
   plt.ylabel("Mean interGene size", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=60000)
   fig_name=OUTPUT+"/"+species+"_mean_interGene_size_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_interGene_size_zoom(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean interGene size per species")
   plt.ylabel("Mean interGene size", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=20000)
   fig_name=OUTPUT+"/"+species+"_mean_interGene_size_per_species_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_overlap_size(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean overlap size per species")
   plt.ylabel("Mean overlap size", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=5000)
   fig_name=OUTPUT+"/"+species+"_mean_overlap_size_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_Exon_minRatio_Overlap(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Min overlap ratio Exon per species")
   plt.ylabel("Mean Min overlap ratio", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=100)
   fig_name=OUTPUT+"/"+species+"_mean_min_overlap_ratio_Exon_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_Exon_maxRatio_Overlap(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Max overlap ratio Exon per species")
   plt.ylabel("Mean Max overlap ratio", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=100)
   fig_name=OUTPUT+"/"+species+"_mean_max_overlap_ratio_Exon_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_Exon_Ratio_Overlap(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean overlap ratio Exon per species")
   plt.ylabel("Mean overlap ratio", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=100)
   fig_name=OUTPUT+"/"+species+"_mean_overlap_ratio_Exon_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_inclusion_size(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean inclusion size per species")
   plt.ylabel("Mean inclusion size", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0,ymax=5000)
   fig_name=OUTPUT+"/"+species+"_mean_inclusion_size_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_nbOverlap_gene(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Overlap number per Gene per species")
   plt.ylabel("Mean Overlap number per Gene", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_mean_nbOverlap_per_gene_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def graph_nbInclusion_gene(dict_to_plot,OUTPUT,species):
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
   plt.title("Mean Inclusion number per Gene per species")
   plt.ylabel("Mean Inclusion number per Gene", fontsize=12.5)
   plt.xlabel("Species", fontsize=12.5)
   plt.xticks(index+width/2., [spe for spe in list_keys], rotation=90, fontsize=10)
#   plt.yticks(np.arange(0, 101, 10), fontsize=16)
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/"+species+"_mean_nbInclusion_per_gene_per_species.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def uniq(seq):
   # Not order preserving
   keys = {}
   for e in seq:
       keys[e] = 1
   return list(keys.keys())

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise

def mean(table):
    return sum(table, 0.0) / len(table)

def variance(table):
    m=mean(table)
    return mean([(x-m)**2 for x in table])

def SD(table):
   return variance(table)**0.5


if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   INPUT_dir=argv[1]
   OUTPUT_dir=argv[2]

   # To be sure than directory have no "/" to the end of the path
   r=search('^(.*)/$', OUTPUT_dir)
   if r:
      OUTPUT_dir=r.group(1)
   # To be sure than directory have no "/" to the end of the path
   r=search('^(.*)/$', INPUT_dir)
   if r:
      INPUT_dir=r.group(1)

   # Create OUTPUT_dir if not existing
   if not path.exists(OUTPUT_dir):
      mkdir_p(OUTPUT_dir)

   # Get list of files contained in input directory
   list_files=listdir(INPUT_dir)

   # Variable 
   nb_contig_tot=0
   nb_gene_tot=0
   nb_exon_tot=0
   nb_overlap_tot=0
   nb_inclusion_tot=0
   nb_gene_overlap_tot=0
   nb_gene_inclusion_tot=0
   nb_gene_overlap_AND_OR_inclusion=0
   nb_gene_overlap_AND_inclusion=0

   # List
   list_nbContig_per_spe=list()
   list_nbG_per_contig=list()
   list_gene_size=list()
   list_gene_size_Exon=list()
   list_interGene_size=list()
   list_overlap_size=list()
   list_inclusion_size=list()
   list_nbOverlap_per_gene=list()
   list_nbInclusion_per_gene=list()
   list_Exon_minRatio_Overlap=list()
   list_Exon_maxRatio_Overlap=list()
   list_Exon_Ratio_Overlap=list()

   # dictionary
   dict_spe_ctgNb={}
   dict_spe_nbG_per_contig={}
   dict_spe_Gene_size={}
   dict_spe_Gene_size_Exon={}
   dict_spe_interGene_size={}
   dict_spe_overlap_size={}
   dict_spe_inclusion_size={}
   dict_spe_nbOverlap_gene={}
   dict_spe_nbInclusion_gene={}
   dict_spe_Exon_minRatio_Overlap={}
   dict_spe_Exon_maxRatio_Overlap={}
   dict_spe_Exon_Ratio_Overlap={}

   # Open GENE output files in Write
   OUTPUT_GENE=OUTPUT_dir+"/ALL_species_GENE_file"
   output_G=open(OUTPUT_GENE,'w')
   # Open OVERLAP output files in Write
   OUTPUT_Overlap=OUTPUT_dir+"/ALL_species_Overlap_file"
   output_Ov=open(OUTPUT_Overlap,'w')
   # Open INCLUSION output files in Write
   OUTPUT_Inclusion=OUTPUT_dir+"/ALL_species_Inclusion_file"
   output_Inc=open(OUTPUT_Inclusion,'w')


   OUTPUT_dir_GRAPH=OUTPUT_dir+"/GRAPH"

   # Browse GENE files contained in INPUT directory
   for gene_file in sorted(list_files):
      r_spe1=search('^(.*)_sorted.txt$', gene_file)
      r_spe2=search('^(.*)_filtered.txt$', gene_file)
      if not (r_spe1 or r_spe2):
         exit("GENE file "+gene_file+" is bad written should be $(name_species)_sorted.txt OR $(name_species)_filtered.txt")
      else:
         name_spe=""
         if r_spe1:
            name_spe=r_spe1.group(1)
         else:
            name_spe=r_spe2.group(1)

         # Create OUTPUT_dir_GRAPH if not existing
         if not path.exists(OUTPUT_dir_GRAPH+"/"+name_spe):
            mkdir_p(OUTPUT_dir_GRAPH+"/"+name_spe)

         # Browse GENE file
         INPUT_GENE_file=INPUT_dir+"/"+gene_file
         input_file=open(INPUT_GENE_file,'r')

         # Variable
         contig_name=""
         nb_contig=0
         nb_gene=0
         start_pos=0
         end_pos=0
         nb_gene_tot_spe=0
         nb_exon_tot_spe=0
         nb_overlap=0
         nb_inclusion=0
         nb_gene_overlap_spe=0
         nb_gene_inclusion_spe=0
         nb_gene_overlap_AND_OR_inclusion_spe=0
         nb_gene_overlap_AND_inclusion_spe=0

         # list() / dict{}
         list_nbG_per_contig_spe=list()
         list_gene_size_spe=list()
         list_gene_size_Exon_spe=list()
         list_interGene_size_spe=list()
         list_overlap_size_spe=list()
         list_inclusion_size_spe=list()
         list_nbOverlap_per_gene_spe=list()
         list_nbInclusion_per_gene_spe=list()
         list_Exon_minRatio_Overlap_spe=list()
         list_Exon_maxRatio_Overlap_spe=list()
         list_Exon_Ratio_Overlap_spe=list()
         # Dictionaries to store gene and determine overlap or/and inclusion in the gene(s) stored in "dic_gene_ENDpos" + Get the number of overlap and inclusion where genes are involved
         dic_gene_ENDpos={}
         dic_gene_pos={}
         dic_gene_before_clean={}
         dic_gene_after_clean={}
         dic_gene_nbInclusion_spe={}
         dic_gene_nbOverlap_spe={}

         
         # Browse sorted GENE file of "name_spe"
         for line in input_file:
            r=search('^([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)\n$', line)
            if r:
#               print(line)
               species=r.group(1)
               contig=r.group(2)
               gene=r.group(3)
               orientation=r.group(4)
               start=r.group(5)
               stop=r.group(6)
               nb_exon=r.group(7)
               exon_pos=r.group(8)

#               print(exon_pos)

               if species!="#species":
#################################
### Get infos on CURRENT GENE ###
#################################
                  bool_overlap=False
                  bool_inclusion=False
                  dic_gene_pos[gene]=start+"\t"+stop

                  # Store gene before cleaning
                  dic_gene_before_clean[gene]=species+"\t"+contig+"\t"+gene+"\t"+orientation+"\t"+start+"\t"+stop+"\t"+nb_exon+"\t"+exon_pos

                  # Store gene size (start 1st exon until stop of last exon)
                  list_gene_size_spe.append(int(stop)-int(start)+1)

                  # Get position of exons
                  list_Exon_pos_CUR_gene=exon_pos.split(":")
                  # print(list_Exon_pos_CUR_gene)

                  # Get cumulative size of exons present in current gene
                  CUR_gene_size_exon=0
                  for exon in list_Exon_pos_CUR_gene:
   #                  print(exon)
                     r_exon=search('^([0-9]*)-([0-9]*)$',exon)
                     if r_exon:
                        start_exon=r_exon.group(1)
                        stop_exon=r_exon.group(2)
                        CUR_gene_size_exon+=int(stop_exon)-int(start_exon)+1
   #                     print(CUR_gene_size_exon)
                     else:
                        exit("Error exon position "+exon+" is bad written!!!")
                  # Store gene size (Sum of all exons size)
                  list_gene_size_Exon_spe.append(CUR_gene_size_exon)


                  # IF NOT THE 1st CONTIG
                  if contig_name!="":
##############################
### If we change of contig ###
##############################
                     if contig_name!=contig:
                        list_nbG_per_contig_spe.append(nb_gene)
                        nb_gene=0
                        nb_contig+=1

                        # Gene stored in dic_gene_ENDpos can't overlap gene on an other contig => Then, clear "dic_gene_ENDpos"
                        dic_gene_ENDpos.clear()

                        # Print new position of gene in dic_gene_after_clean in GENE file output
                        for gene_ID in list(dic_gene_after_clean.keys()):
                           # Case or all exon from a gene are overlap by gene overlaping this gene
                           if dic_gene_after_clean[gene_ID]=="":
                              output_G.write(dic_gene_before_clean[gene_ID]+"\t"+dic_gene_pos[gene_ID]+"\t0\t\n")
                              del dic_gene_pos[gene_ID]
                              del dic_gene_before_clean[gene_ID]
                              del dic_gene_after_clean[gene_ID]
                           else:
                              r_gene=search('^([^-]*).*-([0-9]*)$',dic_gene_after_clean[gene_ID])
                              if r_gene:
                                 exon_nb=len(dic_gene_after_clean[gene_ID].split(":"))
                                 output_G.write(dic_gene_before_clean[gene_ID]+"\t"+dic_gene_pos[gene_ID]+"\t"+str(exon_nb)+"\t"+dic_gene_after_clean[gene_ID]+"\n")
                                 del dic_gene_pos[gene_ID]
                                 del dic_gene_before_clean[gene_ID]
                                 del dic_gene_after_clean[gene_ID]
                              else:
                                 exit("X Exon positions of gene "+gene_ID+" are bad written: "+dic_gene_after_clean[gene_ID])

####################################
### If we DON'T change of contig ###
####################################
                     else:
   #                     if len(dic_gene_ENDpos)>2:
   #                        print("Number of gene stored in dic_gene_ENDpos = "+str(len(dic_gene_ENDpos)))
                        # Browse list of genes to see if we have reach the end of the gene (if not there is an ovelap or an inclusion)
                        for gene_ID in list(dic_gene_ENDpos.keys()):
                           # Get STOP position of Gene stored in dic_gene_ENDpos
                           stop_gene=""
                           r_stop=search('^.*-([0-9]*)$',dic_gene_ENDpos[gene_ID])
                           if r_stop:
                              stop_gene=r_stop.group(1)
                              if int(stop_gene)<int(start):
                                 del dic_gene_ENDpos[gene_ID]
                                 if gene_ID in list(dic_gene_after_clean.keys()):
                                    # Print new position of gene_ID in GENE file output
                                    # Case or all exon from a gene are overlap by gene overlaping this gene
                                    if dic_gene_after_clean[gene_ID]=="":
                                       output_G.write(dic_gene_before_clean[gene_ID]+"\t"+dic_gene_pos[gene_ID]+"\t0\t\n")
                                       del dic_gene_pos[gene_ID]
                                       del dic_gene_before_clean[gene_ID]
                                       del dic_gene_after_clean[gene_ID]
                                    else:
                                       if match('^([^-]*).*-([0-9]*)$',dic_gene_after_clean[gene_ID]):
                                          exon_nb=len(dic_gene_after_clean[gene_ID].split(":"))
                                          output_G.write(dic_gene_before_clean[gene_ID]+"\t"+dic_gene_pos[gene_ID]+"\t"+str(exon_nb)+"\t"+dic_gene_after_clean[gene_ID]+"\n")
                                          del dic_gene_pos[gene_ID]
                                          del dic_gene_before_clean[gene_ID]
                                          del dic_gene_after_clean[gene_ID]
                                       else:
                                          exit("XX Exon positions of gene "+gene_ID+" are bad written: "+dic_gene_after_clean[gene_ID])

##################################################
### In this case, gene_ID overlap current gene ###
##################################################
                              else:
                                 bool_overlap=True
#######################################################
### In this case, current gene included in gene_ID
#######################################################
                                 if int(stop)<=int(stop_gene):
                                    bool_inclusion=True
##########################################################################################################################################
### TO DO => To reduce number of genes included. Several possibilities:                                                              
###   1- Introduce threshold to determine INCLUDED gene => Gene is included if included gene has a size >= X% of the Including gene   
###   2- Analyze if exon overlap between genes INCLUDED gene => If overlap then gene is included else, it's NOT
###   Use definition of article:   Overlapping genes in vertebrate genomes, Computational Biology and Chemistry, 2005           
##########################################################################################################################################
                                    dic_gene_before_clean.pop(gene, None)
                                    # Size inclusion = size of gene included
                                    inclusion_size=int(stop)-int(start)+1
                                    list_inclusion_size_spe.append(inclusion_size)
   #                                 print("Gene "+gene+" is included in gene "+gene_ID)
                                    # Write gene in OUTPUT file containg gene involved in Inclusion (Format: Inclusor_Gene_ID Included_Gene_ID) 
                                    output_Inc.write(gene_ID+"\t"+gene+"\n")
                                    nb_inclusion+=1
                                    # Store gene_ID in dictionary of gene involved in Inclusion
                                    if gene_ID in dic_gene_nbInclusion_spe:
                                       dic_gene_nbInclusion_spe[gene_ID]+=1
                                    else:
                                       dic_gene_nbInclusion_spe[gene_ID]=1
                                    # Store gene in dictionary of gene involved in Inclusion
                                    if gene in dic_gene_nbInclusion_spe:
                                       dic_gene_nbInclusion_spe[gene]+=1
                                    else:
                                       dic_gene_nbInclusion_spe[gene]=1
                                    
                                 else:
                                    # Overlap size on all Gene (Exon+Intron)
                                    overlap_size=int(stop_gene)-int(start)+1
                                    list_overlap_size_spe.append(overlap_size)

                                    # Overlap size of Exon of a gene on complete other gene (Exon only)
                                    l_ov_gene1=0   # Gene in 5' (stored gene)
                                    l_ov_gene2=0   # Gene in 3' (current gene)

                                    # Get position of exons in Gene stored in dic_gene_ENDpos
                                    list_Exon_pos_STORE_gene=dic_gene_ENDpos[gene_ID].split(":")

                                    # Get the size of Exon overlap of gene1 (stored) on gene2 (current)
                                    STORE_gene_size_exon=0
                                    for exon in list_Exon_pos_STORE_gene:
                                       r_exon=search('^([0-9]*)-([0-9]*)$',exon)
                                       if r_exon:
                                          start_exon=r_exon.group(1)
                                          stop_exon=r_exon.group(2)
                                          STORE_gene_size_exon+=int(stop_exon)-int(start_exon)+1
                                          # If start exon of stored gene > start on current gene: all the exon overlap
                                          if int(start_exon)>int(start):
                                             l_ov_gene1+=int(stop_exon)-int(start_exon)+1
                                          else:
                                             # If not and stop exon of stored gene > start of current gene: part of the exon overlap
                                             if int(stop_exon)>int(start):
                                                l_ov_gene1+=int(stop_exon)-int(start)+1
                                       else:
                                          exit("Error exon position "+exon+" is bad written!!!")

                                    # If not then the overlapping gene was also an included gene so he has to be removed
                                    if gene_ID in list(dic_gene_after_clean.keys()):
                                       # Get position of exons in Gene stored in dic_gene_ENDpos
                                       list_Exon_pos_STORE_gene=dic_gene_after_clean[gene_ID].split(":")
   #                                    print list_Exon_pos_STORE_gene

                                       # To change position of exons after overlap cleaning
                                       new_exon_pos_store=""
                                       for exon in list_Exon_pos_STORE_gene:
   #                                       print exon
                                          r_exon=search('^([0-9]*)-([0-9]*)$',exon)
   #                                       print "gene_ID "+gene_ID+" exon "+exon+" start "+start+" stop "+stop
                                          if r_exon:
                                             start_exon=r_exon.group(1)
                                             stop_exon=r_exon.group(2)
                                             if int(start_exon)<int(start):
                                                # If stop exon of stored gene > start of current gene: part of the exon overlap
                                                if int(stop_exon)>int(start):
                                                   end_exon=int(start)-1
                                                   if new_exon_pos_store=="":
                                                      new_exon_pos_store+=start_exon+"-"+str(end_exon)
                                                   else:
                                                      new_exon_pos_store+=" "+start_exon+"-"+str(end_exon)
                                                # Else no overlap of the exon of stored gene on current gene
                                                else:
                                                   if new_exon_pos_store=="":
                                                      new_exon_pos_store+=exon
                                                   else:
                                                      new_exon_pos_store+=":"+exon
                                          else:
                                             exit("Error exon position "+exon+" is bad written!!!")

                                       # Store new positions of gene_ID (stored gene) in dic_gene_after_clean
   #                                    print "NEW "+new_exon_pos_store
                                       r_gene=search('^([0-9]*)\t[0-9]*$',dic_gene_pos[gene_ID])
                                       if r_gene:
                                          start_G=r_gene.group(1)
                                          dic_gene_pos[gene_ID]=start_G+"\t"+str(int(start)-1)
   #                                       print gene_ID+" "+start_G+"\t"+str(int(start)-1)
                                       else:
                                          exit("dic_gene_pos[gene_ID] is bad written")
                                       dic_gene_after_clean[gene_ID]=new_exon_pos_store

                                    # Get the size of Exon overlap of gene2 (current gene) on gene1 (stored)
                                    new_exon_pos_cur=""
                                    for exon in list_Exon_pos_CUR_gene:
                                       r_exon=search('^([0-9]*)-([0-9]*)$',exon)
                                       if r_exon:
                                          start_exon=r_exon.group(1)
                                          stop_exon=r_exon.group(2)
                                          # If start exon of current gene < stop position of stored gene: Exon overlap stored gene
                                          if int(start_exon)<int(stop_gene):
                                             # If stop exon < stop position of stored gene: All exon overlap
                                             if int(stop_exon)<int(stop_gene):
                                                l_ov_gene2+=int(stop_exon)-int(start_exon)+1
                                             # If stop exon > stop position of stored gene: Part of exon overlap
                                             else:
                                                l_ov_gene2+=int(stop_gene)-int(start_exon)+1
                                                start_pos_exon_cur=int(stop_gene)+1
                                                if new_exon_pos_cur=="":
                                                   new_exon_pos_cur+=str(start_pos_exon_cur)+"-"+stop_exon
                                                else:
                                                   new_exon_pos_cur+=" "+str(start_pos_exon_cur)+"-"+stop_exon
                                          # Else exon of current gene don't overlap store gene'
                                          else:
                                                if new_exon_pos_cur=="":
                                                   new_exon_pos_cur+=exon
                                                else:
                                                   new_exon_pos_cur+=":"+exon
                                       else:
                                          exit("Error exon position "+exon+" is bad written!!!")
                                    # Store new positions of gene (current gene) in dic_gene_after_clean
                                    r_gene=search('^[0-9]*\t([0-9]*)$',dic_gene_pos[gene])
                                    if r_gene:
                                       stop_G=r_gene.group(1)
                                       dic_gene_pos[gene]=str(int(stop_gene)+1)+"\t"+stop_G
   #                                    print gene+" "+str(int(stop_gene)+1)+"\t"+stop_G
                                    else:
                                       exit("dic_gene_pos[gene_ID] os bad written")
                                    if not bool_inclusion:
                                       dic_gene_after_clean[gene]=new_exon_pos_cur

   #                                 print gene_ID+" "+str(l_ov_gene1)+" | "+gene+" "+str(l_ov_gene2)

                                    Overlap_ratio_gene1=float(l_ov_gene1)*100/float(STORE_gene_size_exon)
                                    Overlap_ratio_gene2=float(l_ov_gene2)*100/float(CUR_gene_size_exon)

   #                                 print "Exon_Overlap_ratio "+gene_ID+" "+str(Overlap_ratio_gene1)+" | Exon_Overlap_ratio "+gene+" "+str(Overlap_ratio_gene2)
   #                                 print "max: "+str(max(Overlap_ratio_gene1,Overlap_ratio_gene2))+" | min: "+str(min(Overlap_ratio_gene1,Overlap_ratio_gene2))

                                    list_Exon_maxRatio_Overlap_spe.append(max(Overlap_ratio_gene1,Overlap_ratio_gene2))
                                    list_Exon_minRatio_Overlap_spe.append(min(Overlap_ratio_gene1,Overlap_ratio_gene2))

   #                                 print "Gene "+gene+" overlap on the end of "+gene_ID+" on "+str(overlap_size)+" nucleotides"
                                    output_Ov.write(gene_ID+"\t"+gene+"\t"+str(overlap_size)+"\n")
                                    nb_overlap+=1
                                    # Store gene_ID in dictionary of gene involved in Overlap
                                    if gene_ID in dic_gene_nbOverlap_spe:
                                       dic_gene_nbOverlap_spe[gene_ID]+=1
                                    else:
                                       dic_gene_nbOverlap_spe[gene_ID]=1
                                    # Store gene in dictionary of gene involved in Overlap
                                    if gene in dic_gene_nbOverlap_spe:
                                       dic_gene_nbOverlap_spe[gene]+=1
                                    else:
                                       dic_gene_nbOverlap_spe[gene]=1
                           else:
                              exit("Exon positions of gene "+gene_ID+" are bad written: "+dic_gene_ENDpos[gene_ID])
                        if not bool_overlap:
                           # If CUR gene isn't overlapping with stored gene(s) exon position are the same as in original file for 5'end
                           dic_gene_after_clean[gene]=exon_pos
                           interGene_size=int(start)-end_pos
                           list_interGene_size_spe.append(interGene_size)

########################################
### If 1ST CONTIG OF CURRENT SPECIES ###
########################################
                  else:
                     nb_contig=1

####################################
### ADD STATS ON CURRENT SPECIES ###
####################################
                  if not bool_inclusion: 
                     dic_gene_ENDpos[gene]=exon_pos
                  nb_exon_tot_spe+=int(nb_exon)
                  nb_gene_tot_spe+=1
                  nb_gene+=1


                  if not bool_inclusion:
                     if not gene in list(dic_gene_after_clean.keys()):
                        dic_gene_after_clean[gene]=exon_pos
                     if not gene in list(dic_gene_pos.keys()):
                        dic_gene_pos[gene]=start+"\t"+stop

                  # Update variables with data from the line
                  contig_name=contig
                  start_pos=int(start)
                  end_pos=int(stop)
            else:
               exit("ERROR, line:\n"+line+" in file "+INPUT_GENE_file+" is incorrectly written\n")

         dic_gene_ENDpos.clear()

         # Print new position of gene in dic_gene_after_clean in GENE file output
         for gene_ID in list(dic_gene_after_clean.keys()):
            # Case or all exon from a gene are overlap by gene overlaping this gene
            if dic_gene_after_clean[gene_ID]=="":
               output_G.write(dic_gene_before_clean[gene_ID]+"\t"+dic_gene_pos[gene_ID]+"\t0\t\n")
               del dic_gene_pos[gene_ID]
               del dic_gene_before_clean[gene_ID]
               del dic_gene_after_clean[gene_ID]
            else:
               gene_pos=""
               if match('^([^-]*).*-([0-9]*)$',dic_gene_after_clean[gene_ID]):
                  nb_exon=len(dic_gene_after_clean[gene_ID].split(":"))
                  output_G.write(dic_gene_before_clean[gene_ID]+"\t"+dic_gene_pos[gene_ID]+"\t"+str(nb_exon)+"\t"+dic_gene_after_clean[gene_ID]+"\n")
                  del dic_gene_pos[gene_ID]
                  del dic_gene_before_clean[gene_ID]
                  del dic_gene_after_clean[gene_ID]
               else:
                  exit("XXX Exon positions of gene "+gene_ID+" are bad written: "+dic_gene_after_clean[gene_ID])


#######################################
#### STATISTICS on CURRENT species ####
#######################################

      # Fill dictionary for the CURRENT species
      for overlap in list(dic_gene_nbOverlap_spe.keys()):
         list_nbOverlap_per_gene_spe.append(dic_gene_nbOverlap_spe[overlap])
      for inclusion in list(dic_gene_nbInclusion_spe.keys()):
         list_nbInclusion_per_gene_spe.append(dic_gene_nbInclusion_spe[inclusion])

      nb_gene_overlap_spe=len(list_nbOverlap_per_gene_spe)
      nb_gene_inclusion_spe=len(list_nbInclusion_per_gene_spe)
      nb_gene_overlap_AND_OR_inclusion_spe=len(uniq(list(dic_gene_nbOverlap_spe.keys())+list(dic_gene_nbInclusion_spe.keys())))
      list_gene_overlap_AND_inclusion=list(set(dic_gene_nbOverlap_spe.keys()).intersection(list(dic_gene_nbInclusion_spe.keys())))
      nb_gene_overlap_AND_inclusion_spe=len(list_gene_overlap_AND_inclusion)

      dict_spe_ctgNb[name_spe]=nb_contig


      # Statistics SUMMARY on the CURRENT species
      print("")
      print("##############################")
      print("### "+name_spe)
      print("##############################")
      print("\t- Contig number:"+str(nb_contig))
      print("\t- Gene number:"+str(nb_gene_tot_spe))
      print("\t- Exon number:"+str(nb_exon_tot_spe))
      print("\t- Average Gene number per contig:"+str(mean(list_nbG_per_contig_spe)))
      print("\t- SD Gene number per contig:"+str(SD(list_nbG_per_contig_spe)))
      print("\t- Average Gene size:"+str(mean(list_gene_size_spe)))
      print("\t- SD Gene size:"+str(SD(list_gene_size_spe)))
      print("\t- Average Gene size (Exon only):"+str(mean(list_gene_size_Exon_spe)))
      print("\t- SD Gene size (Exon only):"+str(SD(list_gene_size_Exon_spe)))
      print("\t- Average interGene size:"+str(mean(list_interGene_size_spe)))
      print("\t- SD interGene size:"+str(SD(list_interGene_size_spe)))
      if len(list_overlap_size_spe)>0:
         print("\t- Overlap number:"+str(nb_overlap))
         print("\t- Gene number involved in Overlap:"+str(nb_gene_overlap_spe))
         print("\t- Average Overlap size:"+str(mean(list_overlap_size_spe)))
         print("\t- SD Overlap size:"+str(SD(list_overlap_size_spe)))
      else:
         print("\t- NO Gene overlap for this species")
      if len(list_inclusion_size_spe)>0:
         print("\t- Inclusion number:"+str(nb_inclusion))
         print("\t- Gene number involved in Inclusion:"+str(nb_gene_inclusion_spe))
         print("\t- Average Inclusion size:"+str(mean(list_inclusion_size_spe)))
         print("\t- SD Inclusion size:"+str(SD(list_inclusion_size_spe)))
      else:
         print("\t- NO Gene inclusion for this species")
      print("\t- Gene number involved in Overlap AND/OR Inclusion:"+str(nb_gene_overlap_AND_OR_inclusion_spe))
      print("\t- Gene number involved in Overlap AND Inclusion:"+str(nb_gene_overlap_AND_inclusion_spe))


      # Print DATA on OVERLAP and INCLUSION
#      print "\nGene(s) involved in OVERLAP and INCLUSION:"
#      for g in list_gene_overlap_AND_inclusion:
#         print "\t"+g

      print("\nOVERLAP:")
      if bool(dic_gene_nbOverlap_spe):
#         for ov in dic_gene_nbOverlap_spe.keys():
#            if dic_gene_nbOverlap_spe[ov]>2:
#               print ov+" : "+str(dic_gene_nbOverlap_spe[ov])
         dic_gene_nbOverlap_spe.clear()
         # Print distribution Gene number per Overlap number
         dist_nbOv_gene_spe={}
         for elem in list_nbOverlap_per_gene_spe:
            if elem in dist_nbOv_gene_spe:
               dist_nbOv_gene_spe[elem]+=1
            else:
               dist_nbOv_gene_spe[elem]=1
         for elem in dist_nbOv_gene_spe:
            print("\t"+str(dist_nbOv_gene_spe[elem])+" gene(s) have "+str(elem)+" overlap(s)")
         dist_nbOv_gene_spe.clear()
      else:
         print("\tNO gene overlap in this species")

      print("\nINCLUSION:")
      if bool(dic_gene_nbInclusion_spe):
#         for inc in dic_gene_nbInclusion_spe.keys():
#            if dic_gene_nbInclusion_spe[inc]>3:
#               print inc+" : "+str(dic_gene_nbInclusion_spe[inc])
         dic_gene_nbInclusion_spe.clear()
         # Print distribution Gene number per Overlap number
         dist_nbInc_gene_spe={}
         for elem in list_nbInclusion_per_gene_spe:
            if elem in dist_nbInc_gene_spe:
               dist_nbInc_gene_spe[elem]+=1
            else:
               dist_nbInc_gene_spe[elem]=1
         for elem in dist_nbInc_gene_spe:
            print("\t"+str(dist_nbInc_gene_spe[elem])+" gene(s) have "+str(elem)+" inclusion(s)")
         dist_nbInc_gene_spe.clear()
      else:
         print("\tNO gene inclusion in this species")


      # Draw histogram for different statistics on Genome Structure
      OUTPUT_dir_GRAPH_Spe=OUTPUT_dir_GRAPH+"/"+name_spe
      mkdir_p(OUTPUT_dir_GRAPH_Spe)

      # Contig number / Gene number distribution
      hist_nbGene_contig(list_nbG_per_contig_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      hist_nbGene_contig_zoom(list_nbG_per_contig_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      hist_nbGene_contig_zoom_zoom(list_nbG_per_contig_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      # Gene size distributon
      hist_Gene_size(list_gene_size_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      hist_Gene_size_zoom(list_gene_size_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      # Gene size distributon (Exon only)
      hist_Gene_size_Exon(list_gene_size_Exon_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
#      hist_Gene_size_Exon_zoom(list_gene_size_Exon_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      # interGene size distributon
      hist_interGene_size(list_interGene_size_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      hist_interGene_size_zoom(list_interGene_size_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      # Overlap size distributon
      if len(list_overlap_size_spe)>0:
         hist_overlap_size(list_overlap_size_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
         hist_overlap_size_zoom(list_overlap_size_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      # Gene ratio in Overlap (Min, Max & ALL)
         hist_Exon_overlap_ratio_min(list_Exon_minRatio_Overlap_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
         hist_Exon_overlap_ratio_max(list_Exon_maxRatio_Overlap_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
         list_Exon_Ratio_Overlap_spe=list_Exon_minRatio_Overlap_spe+list_Exon_maxRatio_Overlap_spe
         hist_Exon_overlap_ratio(list_Exon_Ratio_Overlap_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      # Inclusion size distributon
      if len(list_inclusion_size_spe)>0:
         hist_inclusion_size(list_inclusion_size_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
         hist_inclusion_size_zoom(list_inclusion_size_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      # Overlap number / Gene number distribution
      if len(list_nbOverlap_per_gene_spe)>0:
         hist_nbOverlap_gene(list_nbOverlap_per_gene_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
#         hist_inclusion_size_zoom(list_nbOverlap_per_gene_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
      # Inclusion number / Gene number distribution
      if len(list_nbInclusion_per_gene_spe)>0:
         hist_nbInclusion_gene(list_nbInclusion_per_gene_spe,OUTPUT_dir_GRAPH_Spe,name_spe)
#         hist_inclusion_size_zoom(list_nbInclusion_per_gene_spe,OUTPUT_dir_GRAPH_Spe,name_spe)




###################################
###   STATISTICS on ALL species ###
###################################

      #Store species list() and variable in ALL species list() and variable

      # Fill variable
      nb_contig_tot+=nb_contig
      nb_gene_tot+=nb_gene_tot_spe
      nb_exon_tot+=nb_exon_tot_spe
      nb_overlap_tot+=nb_overlap
      nb_inclusion_tot+=nb_inclusion
      nb_gene_overlap_tot+=nb_gene_overlap_spe
      nb_gene_inclusion_tot+=nb_gene_inclusion_spe
      nb_gene_overlap_AND_OR_inclusion+=nb_gene_overlap_AND_OR_inclusion_spe
      nb_gene_overlap_AND_inclusion+=nb_gene_overlap_AND_inclusion_spe

      # Fill list()
      list_nbContig_per_spe.append(nb_contig)
      list_nbG_per_contig+=list_nbG_per_contig_spe
      list_gene_size+=list_gene_size_spe
      list_gene_size_Exon+=list_gene_size_Exon_spe
      list_interGene_size+=list_interGene_size_spe
      list_overlap_size+=list_overlap_size_spe
      list_Exon_minRatio_Overlap+=list_Exon_minRatio_Overlap_spe
      list_Exon_maxRatio_Overlap+=list_Exon_maxRatio_Overlap_spe
      list_Exon_Ratio_Overlap+=list_Exon_Ratio_Overlap_spe
      list_inclusion_size+=list_inclusion_size_spe
      list_nbOverlap_per_gene+=list_nbOverlap_per_gene_spe
      list_nbInclusion_per_gene+=list_nbInclusion_per_gene_spe

      # Fill dictionary for statistics on ALL SPECIES
      dict_spe_nbG_per_contig[species]=[mean(list_nbG_per_contig_spe),SD(list_nbG_per_contig_spe)]
      dict_spe_Gene_size[species]=[mean(list_gene_size_spe),SD(list_gene_size_spe)]
      dict_spe_Gene_size_Exon[species]=[mean(list_gene_size_Exon_spe),SD(list_gene_size_Exon_spe)]
      dict_spe_interGene_size[species]=[mean(list_interGene_size_spe),SD(list_interGene_size_spe)]
      if len(list_overlap_size_spe)>0:
         dict_spe_overlap_size[species]=[mean(list_overlap_size_spe),SD(list_overlap_size_spe)]
         dict_spe_Exon_minRatio_Overlap[species]=[mean(list_Exon_minRatio_Overlap_spe),SD(list_Exon_minRatio_Overlap_spe)]
         dict_spe_Exon_maxRatio_Overlap[species]=[mean(list_Exon_maxRatio_Overlap_spe),SD(list_Exon_maxRatio_Overlap_spe)]
         dict_spe_Exon_Ratio_Overlap[species]=[mean(list_Exon_Ratio_Overlap_spe),SD(list_Exon_Ratio_Overlap_spe)]
      if len(list_inclusion_size_spe)>0:
         dict_spe_inclusion_size[species]=[mean(list_inclusion_size_spe),SD(list_inclusion_size_spe)]
      if len(list_nbOverlap_per_gene_spe)>0:
         dict_spe_nbOverlap_gene[species]=[mean(list_nbOverlap_per_gene_spe),SD(list_nbOverlap_per_gene_spe)]
      if len(list_nbInclusion_per_gene_spe)>0:
         dict_spe_nbInclusion_gene[species]=[mean(list_nbInclusion_per_gene_spe),SD(list_nbInclusion_per_gene_spe)]


   output_G.close()
   output_Ov.close()
   output_Inc.close()   


   # Print summary available on ALL species
   print("\n###########\nALL_SPECIES\n###########")
   print("\t- Contig number:"+str(nb_contig_tot))
   print("\t- Gene number:"+str(nb_gene_tot))
   print("\t- Exon number:"+str(nb_exon_tot))
   print("\t- Average Contig number per species:"+str(mean(list_nbContig_per_spe)))
   print("\t- SD Contig number per species:"+str(SD(list_nbContig_per_spe)))
   print("\t- Average Gene number per contig:"+str(mean(list_nbG_per_contig)))
   print("\t- SD Gene number per contig:"+str(SD(list_nbG_per_contig)))
   print("\t- Average Gene size:"+str(mean(list_gene_size)))
   print("\t- SD Gene size:"+str(SD(list_gene_size_spe)))
   print("\t- Average Gene size (Exon only):"+str(mean(list_gene_size_Exon)))
   print("\t- SD Gene size (Exon only):"+str(SD(list_gene_size_Exon)))
   print("\t- Average interGene size:"+str(mean(list_interGene_size)))
   print("\t- SD interGene size:"+str(SD(list_interGene_size)))
   if len(list_overlap_size)>0:
      print("\t- Overlap number:"+str(nb_overlap_tot))
      print("\t- Gene number involved in Overlap:"+str(nb_gene_overlap_tot))
      print("\t- Average Overlap size:"+str(mean(list_overlap_size)))
      print("\t- SD Overlap size:"+str(SD(list_overlap_size)))
   else:
      print("\t- NO Gene overlap for this species")
   if len(list_inclusion_size)>0:
      print("\t- Inclusion number:"+str(nb_inclusion_tot))
      print("\t- Gene number involved in Inclusion:"+str(nb_gene_inclusion_tot))
      print("\t- Average Inclusion size:"+str(mean(list_inclusion_size)))
      print("\t- SD Inclusion size:"+str(SD(list_inclusion_size)))
   else:
      print("\t- NO Gene inclusion for this species")
   print("\t- Gene number involved in Overlap AND/OR Inclusion:"+str(nb_gene_overlap_AND_OR_inclusion))
   print("\t- Gene number involved in Overlap AND Inclusion:"+str(nb_gene_overlap_AND_inclusion))


   # Draw histogram for different statistics on Genome Structure
   OUTPUT_dir_GRAPH_Spe=OUTPUT_dir_GRAPH+"/ALL_species"
   mkdir_p(OUTPUT_dir_GRAPH_Spe)

   # Contig number / Gene number distributon
   hist_nbGene_contig(list_nbG_per_contig,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   hist_nbGene_contig_zoom(list_nbG_per_contig,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   hist_nbGene_contig_zoom_zoom(list_nbG_per_contig,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   # Gene size distributon
   hist_Gene_size(list_gene_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   hist_Gene_size_zoom(list_gene_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   # Gene size distributon (Exon only)
   hist_Gene_size_Exon(list_gene_size_Exon,OUTPUT_dir_GRAPH_Spe,"ALL_species")
#   hist_Gene_size_Exon_zoom(list_gene_size_Exon,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   # interGene size distributon
   hist_interGene_size(list_interGene_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   hist_interGene_size_zoom(list_interGene_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   # Overlap size distributon
   if len(list_overlap_size)>0:
      hist_overlap_size(list_overlap_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
      hist_overlap_size_zoom(list_overlap_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   # Gene ratio in Overlap (Min, Max & ALL)
      hist_Exon_overlap_ratio_min(list_Exon_minRatio_Overlap,OUTPUT_dir_GRAPH_Spe,"ALL_species")
      hist_Exon_overlap_ratio_max(list_Exon_maxRatio_Overlap,OUTPUT_dir_GRAPH_Spe,"ALL_species")
      list_Exon_Ratio_Overlap=list_Exon_minRatio_Overlap+list_Exon_maxRatio_Overlap
      hist_Exon_overlap_ratio(list_Exon_Ratio_Overlap,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   # Inclusion size distributon
   if len(list_inclusion_size)>0:
      hist_inclusion_size(list_inclusion_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
      hist_inclusion_size_zoom(list_inclusion_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   # Overlap number / Gene number distribution
   if len(list_nbOverlap_per_gene)>0:
      hist_nbOverlap_gene(list_nbOverlap_per_gene,OUTPUT_dir_GRAPH_Spe,"ALL_species")
#      hist_nbOverlap_gene_zoom(list_nbOverlap_per_gene,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   # Inclusion number / Gene number distribution
   if len(list_nbInclusion_per_gene)>0:
      hist_nbInclusion_gene(list_nbInclusion_per_gene,OUTPUT_dir_GRAPH_Spe,"ALL_species")
#      hist_nbInclusion_gene_zoom(list_nbInclusion_per_gene,OUTPUT_dir_GRAPH_Spe,"ALL_species")

   # Comparison between species
   graph_ctgNb(dict_spe_ctgNb,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_nbGene_contig(dict_spe_nbG_per_contig,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_nbGene_contig_zoom(dict_spe_nbG_per_contig,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_nbGene_contig_zoom_zoom(dict_spe_nbG_per_contig,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_Gene_size(dict_spe_Gene_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_Gene_size_zoom(dict_spe_Gene_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_Gene_size_zoom_zoom(dict_spe_Gene_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_Gene_size_Exon(dict_spe_Gene_size_Exon,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_interGene_size(dict_spe_interGene_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_interGene_size_zoom(dict_spe_interGene_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_overlap_size(dict_spe_overlap_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_Exon_minRatio_Overlap(dict_spe_Exon_minRatio_Overlap,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_Exon_maxRatio_Overlap(dict_spe_Exon_maxRatio_Overlap,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_Exon_Ratio_Overlap(dict_spe_Exon_Ratio_Overlap,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_inclusion_size(dict_spe_inclusion_size,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_nbOverlap_gene(dict_spe_nbOverlap_gene,OUTPUT_dir_GRAPH_Spe,"ALL_species")
   graph_nbInclusion_gene(dict_spe_nbInclusion_gene,OUTPUT_dir_GRAPH_Spe,"ALL_species")


   # Sort OUTPUT annotation gene file by species, then contigs and genes.
   out=open("buffer_file","w")
   out.write("#species\tctg\tgene\torientation_gene\tstart_gene_old\tend_gene_old\t#exons_old\texons_position_old\t#exons_new\texons_position_new\n")
   out.close()
   command_line="sort -k1d,1d -k2d,2d -k6n,6n "+OUTPUT_GENE+" >> buffer_file; mv buffer_file "+OUTPUT_GENE
   subprocess.call(command_line,shell=True)

   end_time = datetime.now()
   print('\nDuration: {}'.format(end_time - start_time))
