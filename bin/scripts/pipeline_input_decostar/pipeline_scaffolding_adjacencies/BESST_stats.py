#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      From BESST file score generate stats on scaffolded contigs links
###
###   INPUT:
###      1- BESST file for contigs links analysis
###      2- OUTPUT directory containing all output file
###
###   OUTPUT files:
###      - Stats file on link number par contigs pair
###      - Stats file on gap size
###
###   Name: BESST_stats.py                      Author: Yoann Anselmetti
###   Creation date: 2015/10/21                 Last modification: 2017/07/19
###


from sys import argv
from re import search
from os import close, path, makedirs, remove
import fileinput
import subprocess
import itertools
import matplotlib
matplotlib.use('Agg')
 
import numpy as np 
from matplotlib import pyplot as plt

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

def hist_gap_size(list_to_plot,OUTPUT):
#   index = np.arange(11)
#   width = 0.5
   plt.hist(list_to_plot,bins=100,range=(-10000,10000))
   plt.title("Distribution of gap size between contigs pair.")
   plt.xlabel("Gap size")
   plt.ylabel("Contigs pair number")
   plt.xlim(xmin=-10000,xmax=10000)
#   plt.xticks(index*(-1000)+width,[i*(-1000) for i in index])
#   plt.xticks(index*1000+width,[i*1000 for i in index])
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/distrib_gap_size.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_gap_size_zoom(list_to_plot,OUTPUT):
#   index = np.arange(11)
#   width = 0.5
   plt.hist(list_to_plot,bins=100,range=(-5000,5000))
   plt.title("Distribution of gap size between contigs pair.")
   plt.xlabel("Gap size")
   plt.ylabel("Contigs pair number")
   plt.xlim(xmin=-5000,xmax=5000)
#   plt.xticks(index*500+width,[i*500 for i in index])
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/distrib_gap_size_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()


def hist_linkNb_contigsPair(list_to_plot,OUTPUT):
   index = np.arange(11)
   width = 0.5
   plt.hist(list_to_plot,bins=100,range=(1,101))
   plt.title("Distribution of Link number per Contigs pair number.")
   plt.xlabel("Link number")
   plt.ylabel("Contigs pair number")
   plt.xlim(xmin=1,xmax=101)
   plt.xticks(index*10+width,[i*10 for i in index])
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/distrib_link_nb_per_contigs_pair.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_linkNb_contigsPair_zoom(list_to_plot,OUTPUT):
   index = np.arange(6)
   width = 0.5

   plt.hist(list_to_plot,bins=25,range=(1,26))
   plt.title("Distribution of Link number per Contigs pair number.")
   plt.xlabel("Link number")
   plt.ylabel("Contigs pair number")
   plt.xlim(xmin=1,xmax=26)
   plt.xticks(index*5+width,[i*5 for i in index])
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/distrib_link_nb_per_contigs_pair_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_linkNb_contigsPair_zoom_zoom(list_to_plot,OUTPUT):
   index = np.arange(11)
   width = 0.5

   plt.hist(list_to_plot,bins=10,range=(1,11))
   plt.title("Distribution of Link number per Contigs pair number.")
   plt.xlabel("Link number")
   plt.ylabel("Contigs pair number")
   plt.xlim(xmin=1,xmax=11)
   plt.xticks(index+width,[i for i in index])
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/distrib_link_nb_per_contigs_pair_zoom_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

# Recovery of input parameters
BESST_file=argv[1]
OUTPUT_DIR=argv[2]

file2=OUTPUT_DIR+"/stats_gap_size"
file3=OUTPUT_DIR+"/stats_contigsPair_link"


# To be sure than directory have no "/" to the end of the path
r=search('^(.*)/$', OUTPUT_DIR)
if r:
   OUTPUT_DIR=r.group(1)

# Create OUTPUT_DIR if not existing
if not path.exists(OUTPUT_DIR):
   mkdir_p(OUTPUT_DIR)


print "\nBrowse BESST score file analyze statistics on link number per contigs pairs AND on gap size"
# Creation of dictionary "dict_link" to link pair of link contigs with list of PE/MP reads linking them
list_link_nb=list()
list_gap_size=list()
gap_size_min=10000
gap_size_max=-10000
link_nb_min=10000
link_nb_max=0
tot_link_nb=0
IN=open(BESST_file,'r')
for line in sorted(IN):
   r=search("([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\t]*)[\t]([^\t]*)",line)
   if r:
      ctg1=r.group(1)
      ori1=r.group(2)
      ctg2=r.group(3)
      ori2=r.group(4)
      gap=r.group(5)
      score1=r.group(6)
      score2=r.group(7)
      link_nb=r.group(8)

      if ctg1!="scf1/ctg1":
         # print line,
         list_gap_size.append(float(gap))
         if gap_size_max<float(gap):
            gap_size_max=float(gap)
         if gap_size_min>float(gap):
            gap_size_min=float(gap)
         list_link_nb.append(int(link_nb))
         if link_nb_max<int(link_nb):
            link_nb_max=int(link_nb)
         if link_nb_min>int(link_nb):
            link_nb_min=int(link_nb)
         tot_link_nb+=int(link_nb)
   else:
      exit("Error line in BESST file score: "+BESST_file+" are incorrectly written!!!")
IN.close()
#remove(file1)


hist_gap_size(list_gap_size,OUTPUT_DIR)
hist_gap_size_zoom(list_gap_size,OUTPUT_DIR)
#hist_gap_size_zoom_zoom(list_gap_size,OUTPUT_DIR)

OUT2=open(file2,'w')
OUT2.write("PAIR_NB\t"+str(len(list_gap_size))+"\n")
OUT2.write("GAP_SIZE_MIN\t"+str(gap_size_min)+"\n")
OUT2.write("GAP_SIZE_MAX\t"+str(gap_size_max)+"\n")
OUT2.write("GAP_SIZE_MEAN\t"+str(mean(list_gap_size))+"\n")
OUT2.write("GAP_SIZE_SD\t"+str(SD(list_gap_size))+"\n")
OUT2.write("GAP_SIZE_DISTRIB\t")
bool_first=True
for gap_size in list_gap_size:
   if bool_first:
      OUT2.write(str(gap_size))
      bool_first=False
   else:
      OUT2.write(" "+str(gap_size))
OUT2.close()


hist_linkNb_contigsPair(list_link_nb,OUTPUT_DIR)
hist_linkNb_contigsPair_zoom(list_link_nb,OUTPUT_DIR)
hist_linkNb_contigsPair_zoom_zoom(list_link_nb,OUTPUT_DIR)

OUT3=open(file3,'w')
OUT3.write("PAIR_NB\t"+str(len(list_link_nb))+"\n")
OUT3.write("LINK_NB\t"+str(tot_link_nb)+"\n")
OUT3.write("LINK_NB_MIN\t"+str(link_nb_min)+"\n")
OUT3.write("LINK_NB_MAX\t"+str(link_nb_max)+"\n")
OUT3.write("LINK_NB_MEAN\t"+str(mean(list_link_nb))+"\n")
OUT3.write("LINK_NB_SD\t"+str(SD(list_link_nb))+"\n")
OUT3.write("LINK_NB_DISTRIB\t")
bool_first=True
for nb_link in list_link_nb:
   if bool_first:
      OUT3.write(str(nb_link))
      bool_first=False
   else:
      OUT3.write(" "+str(nb_link))
OUT3.close()

print "\nStatistics on gap size:"
print "   - There are "+str(len(list_gap_size))+" pairs of contigs linked:"
print "      + Min gap size: "+str(gap_size_min)
print "      + Max gap size: "+str(gap_size_max)
print "      + Mean gap size: "+str(mean(list_gap_size))
print "      + SD gap size: "+str(SD(list_gap_size))

print "\nStatistics on linked contigs:"
print "   - There are "+str(len(list_link_nb))+" pairs of contigs linked by "+str(tot_link_nb)+" paired reads:"
print "      + Min number of links per contigs pair: "+str(link_nb_min)
print "      + Max number of links per contigs pair: "+str(link_nb_max)
print "      + Mean number of links per contigs pair: "+str(mean(list_link_nb))
print "      + SD number of links per contigs pair: "+str(SD(list_link_nb))
