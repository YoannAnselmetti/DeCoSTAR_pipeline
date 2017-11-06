#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Produce statistics on mapping from BAM files
###
###   INPUT:
###      1- BAM file to analyze for statistics on mapping
###      2- OUTPUT directory containing all output file
###      3- Orientation library
###
###   OUTPUT files:
###      - PE linking ctg | Format: PE_ID ctg1 ctg2
###      - ctg linked with list PE | Format: ctg1 ctg2 list_PE links_nb
###      - Stats file on insert size of PE/MP reads from BAM file
###      - Insert size distribution histogram of PE/MP reads from BAM file
###      - Statisics file on contigs pairs links
###
###   Name: BAM_stats.py                     Author: Yoann Anselmetti
###   Creation date: 2015/10/13              Last modification: 2017/03/01
###


from sys import argv, stdout
from re import search
from os import close, path, makedirs, remove
import fileinput
import subprocess
import itertools
import matplotlib
matplotlib.use('Agg')
 
import numpy as np 
from matplotlib import pyplot as plt

stdout.flush()

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


# def plot_insert_size_distrib(list_insert_size,list_FR_occ,list_RF_occ,list_Tandem_occ,bool_ori,OUTPUT):
#    if bool_ori=="no_tandem":
#       p1 = plt.fill(list_insert_size, list_FR_occ, 'r-')
#       p2 = plt.fill(list_insert_size, list_RF_occ, 'g-')
#       plt.legend((p1[0], p2[0]),('FR', 'RF'), loc='upper right', fontsize=15)
#    elif bool_ori=="tandem":
#       p1 = plt.fill(list_insert_size, list_FR_occ, 'r-')
#       p2 = plt.fill(list_insert_size, list_RF_occ, 'g-')
#       p3 = plt.fill(list_insert_size, list_Tandem_occ, 'b-')
#       plt.legend((p1[0], p2[0], p3[0]),('FR', 'RF', 'Tandem'), loc='upper right', fontsize=15)
#    else:
#       exit("Error: bool_ori should have a value!!!")
#    plt.title("Distribution of insert size per library orientation.")
#    plt.xlabel("Insert size")
#    plt.ylabel("Paired reads number")
#    fig_name=OUTPUT+"/distrib_insert_size_per_library_orientation.pdf"
#    plt.tight_layout()
#    plt.savefig(fig_name,format='pdf')
#    plt.cla()


def hist_linkNb_contigsPair(list_to_plot,OUTPUT):
   index = np.arange(11)
   width = 0.5
   plt.hist(list_to_plot,bins=100,range=(1,1001))
   plt.title("Distribution of Link number per Contigs pair number.")
   plt.xlabel("Link number")
   plt.ylabel("Contigs pair number")
   plt.xlim(xmin=1,xmax=1001)
   plt.xticks(index*100+(width),[i*100 for i in index])
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/distrib_link_nb_per_contigs_pair.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_linkNb_contigsPair_zoom(list_to_plot,OUTPUT):
   index = np.arange(11)
   width = 0.5
   plt.hist(list_to_plot,bins=100,range=(1,101))
   plt.title("Distribution of Link number per Contigs pair number.")
   plt.xlabel("Link number")
   plt.ylabel("Contigs pair number")
   plt.xlim(xmin=1,xmax=101)
   plt.xticks(index*10+(width),[i*10 for i in index])
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/distrib_link_nb_per_contigs_pair_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

def hist_linkNb_contigsPair_zoom_zoom(list_to_plot,OUTPUT):
   index = np.arange(6)
   width = 0.5
   plt.hist(list_to_plot,bins=25,range=(1,26))
   plt.title("Distribution of Link number per Contigs pair number.")
   plt.xlabel("Link number")
   plt.ylabel("Contigs pair number")
   plt.xlim(xmin=1,xmax=26)
   plt.xticks(index*5+(width),[i*5 for i in index])
   plt.ylim(ymin=0)
   fig_name=OUTPUT+"/distrib_link_nb_per_contigs_pair_zoom_zoom.pdf"
   plt.tight_layout()
   plt.savefig(fig_name,format='pdf')
   plt.cla()

# Recovery of input parameters
BAM_file=argv[1]
OUTPUT_DIR=argv[2]
orientation=argv[3]

file1=OUTPUT_DIR+"/PE_contigsLinked_file"
file2=OUTPUT_DIR+"/contigsLinked_listPE_file"
file3=OUTPUT_DIR+"/stats_contigsPair_link"

CollectInsertSizeMetrics="/share/apps/bin/picard-tools-1.61/CollectInsertSizeMetrics.jar"
# samtools="/share/apps/bin/samtools-1.3/bin/samtools"
   
# To be sure than directory have no "/" to the end of the path
r=search('^(.*)/$', OUTPUT_DIR)
if r:
   OUTPUT_DIR=r.group(1)

# Create OUTPUT_DIR if not existing
if not path.exists(OUTPUT_DIR):
   mkdir_p(OUTPUT_DIR)

print "\nCreate statictics file and histogram distribution of insert size of PE/MP reads contained in BAM file "+BAM_file+" with picard-tools"
# Creation statictics file and histogram distribution of insert size of PE/MP reads contained in BAM file using picard-tools
if orientation=="fr":
   max_insert_size=500
   bash3="module load java1.8; module load R-3.3.1; java -jar "+CollectInsertSizeMetrics+" I="+BAM_file+" O="+OUTPUT_DIR+"/stats_insert_size_"+str(max_insert_size)+" H="+OUTPUT_DIR+"/hist_insert_size_"+str(max_insert_size)+" W="+str(max_insert_size)
   print "=> "+bash3
   subprocess.check_output(bash3,shell=True)

   max_insert_size=100000
   bash3="module load java1.8; module load R-3.3.1; java -jar "+CollectInsertSizeMetrics+" I="+BAM_file+" O="+OUTPUT_DIR+"/stats_insert_size_"+str(max_insert_size)+" H="+OUTPUT_DIR+"/hist_insert_size_"+str(max_insert_size)+" W="+str(max_insert_size)
   print "=> "+bash3
   subprocess.check_output(bash3,shell=True)


elif orientation=="rf":
   max_insert_size=10000
   bash3="module load java1.8; module load R-3.3.1; java -jar "+CollectInsertSizeMetrics+" I="+BAM_file+" O="+OUTPUT_DIR+"/stats_insert_size_"+str(max_insert_size)+" H="+OUTPUT_DIR+"/hist_insert_size_"+str(max_insert_size)+" W="+str(max_insert_size)
   print "=> "+bash3
   subprocess.check_output(bash3,shell=True)

   # Generation of insert size distribution for the different orientations in RF library
   stats_file=OUTPUT_DIR+"/stats_insert_size_"+str(max_insert_size)




   #############################################################################
   ### NO MORE NEED SINCE PICARD TOOLS DEAL LIBRARY WITH SEVERAL ORIENTATION ###
   #############################################################################
   # print "   => Generation of insert size distribution for the different orientations in RF library from "+stats_file+" file"
   # insert_size_file=open(stats_file,'r')
   # bool_ori=""
   # list_insert_size=list()
   # list_FR_occ=list()
   # list_RF_occ=list()
   # list_Tandem_occ=list()
   # for line in insert_size_file:
   #    if bool_ori == "tandem" :
   #       r_tandem=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\n]*)\n$",line)
   #       if r_tandem:
   #          insert_size=int(r_tandem.group(1))
   #          FR_occ=int(r_tandem.group(2))
   #          RF_occ=int(r_tandem.group(3))
   #          Tandem_occ=int(r_tandem.group(4))

   #          list_insert_size.append(insert_size)
   #          list_FR_occ.append(FR_occ)
   #          list_RF_occ.append(RF_occ)
   #          list_Tandem_occ.append(Tandem_occ)
   #    elif bool_ori=="no_tandem":
   #       r_notandem=search("^([^\t]*)\t([^\t]*)\t([^\n]*)\n$",line)
   #       if r_notandem:
   #          insert_size=int(r_notandem.group(1))
   #          FR_occ=int(r_notandem.group(2))
   #          RF_occ=int(r_notandem.group(3))

   #          list_insert_size.append(insert_size)
   #          list_FR_occ.append(FR_occ)
   #          list_RF_occ.append(RF_occ)
   #    if "insert_size\tAll_Reads.fr_count\tAll_Reads.rf_count\tAll_Reads.tandem_count" in line:
   #       bool_ori="tandem"
   #    elif "insert_size\tAll_Reads.fr_count\tAll_Reads.rf_count" in line:
   #       bool_ori="no_tandem"
   # plot_insert_size_distrib(list_insert_size,list_FR_occ,list_RF_occ,list_Tandem_occ,bool_ori,OUTPUT_DIR)


else:
   exit("Error the 3rd parameter should be \"fr\" for Paired-End or fosill library OR \"rf\" for Mate-Pairs library!!!")


print "\nCreation of buffer file where each line correspond to a read linking 2 contigs | Format: ID_PE_reads ctg1 ctg2"
# Creation of buffer file where each line correspond to a read linking 2 contigs | Format: ID_PE_reads ctg1 ctg2
bash4="samtools view "+BAM_file+" | awk '$7!=\"=\" {print $0}' | grep -v -P \"\\t\\*\\t\" |  awk \'{print $1,$3,$7}\' | sort > "+file1
print "=> "+bash4
subprocess.check_output(bash4,shell=True)


print "\nBrowse buffer file to create dictionary of contigs pair with list of PE/MP reads linking the pair"
# Creation of dictionary "dict_link" to link pair of link contigs with list of PE/MP reads linking them
dict_link={}
IN=open(file1,'r')
for line in IN:
   no_PE=True
   r=search("^([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)\n$",line)
   if r:
      PE=r.group(1)
      ctg1=r.group(2)
      ctg2=r.group(3)
      link1=ctg1+"\t"+ctg2
      link2=ctg2+"\t"+ctg1
      if link1 in dict_link:
         for elem in dict_link[link1]:
            # Paired reads have the same ID, here if one the the paired reads have been added in the dictionary the other will be filtered
            if elem==PE:
              no_PE=False
         if no_PE:
            dict_link[link1].append(PE)
      elif link2 in dict_link:
         for elem in dict_link[link2]:
            if elem==PE:
              no_PE=False
         if no_PE:
            dict_link[link2].append(PE)
      else:
         dict_link[link1]=list()
         dict_link[link1].append(PE)
IN.close()
remove(file1)

list_link_nb=list()
tot_link_nb=0
link_nb_min=10000
link_nb_max=0
print "\nBrowse dictionary previously created to generate contigs link file | Format: ctg1 ctg2 list_PE/MP_reads nb_links"
# Browse dictionary "dict_link" previously created to generate contigs link file | Format: ctg1 ctg2 list_PE_reads nb_links
OUT=open(file2,'w')
for key in dict_link.keys():
#   print key+"\t",
   OUT.write(key+"\t")
   first_pe=True
   for pe_id in dict_link[key]:
      if first_pe:
#         print pe_id,
         OUT.write(pe_id)
         first_pe=False
      else:
#         print " "+pe_id,
         OUT.write(" "+pe_id)
#   print "\t"+str(len(dict_link[key]))
   OUT.write("\t"+str(len(dict_link[key]))+"\n")
   list_link_nb.append(len(dict_link[key]))
   if link_nb_max<len(dict_link[key]):
      link_nb_max=len(dict_link[key])
   if link_nb_min>len(dict_link[key]):
      link_nb_min=len(dict_link[key])
   tot_link_nb+=len(dict_link[key])
OUT.close()
remove(file2)

print "\nCreation distribution graph for links number per contigs pair"
hist_linkNb_contigsPair(list_link_nb,OUTPUT_DIR)
hist_linkNb_contigsPair_zoom(list_link_nb,OUTPUT_DIR)
hist_linkNb_contigsPair_zoom_zoom(list_link_nb,OUTPUT_DIR)

print "\nCreation statitics file on links number per contigs pair"
OUT2=open(file3,'w')
OUT2.write("PAIR_NB\t"+str(len(list_link_nb))+"\n")
OUT2.write("LINK_NB\t"+str(tot_link_nb)+"\n")
OUT2.write("LINK_NB_MIN\t"+str(link_nb_min)+"\n")
OUT2.write("LINK_NB_MAX\t"+str(link_nb_max)+"\n")
OUT2.write("LINK_NB_MEAN\t"+str(mean(list_link_nb))+"\n")
OUT2.write("LINK_NB_SD\t"+str(SD(list_link_nb))+"\n")
OUT2.write("LINK_NB_DISTRIB\t")
bool_first=True
for nb_link in list_link_nb:
   if bool_first:
      OUT2.write(str(nb_link))
      bool_first=False
   else:
      OUT2.write(" "+str(nb_link))
OUT2.close()
remove(file3)

print "\nStatistics on linked contigs:"
print "   - There are "+str(len(list_link_nb))+" pairs of contigs linked by "+str(tot_link_nb)+" paired reads:"
print "      + Min number of links per contigs pair: "+str(link_nb_min)
print "      + Max number of links per contigs pair: "+str(link_nb_max)
print "      + Mean number of links per contigs pair: "+str(mean(list_link_nb))
print "      + SD number of links per contigs pair: "+str(SD(list_link_nb))

print "\nStatistics reads mapping:"
# Get number of reads total and unmapped reads
bash1="samtools view "+BAM_file+" | wc -l"
print "=> Compute total number of reads: "+bash1
tot_nb_read=subprocess.check_output(bash1,shell=True)
bash2="samtools view "+BAM_file+" | grep -P \"\\t\\*\\t\" | wc -l"
print "=> Compute number of unmapped reads: "+bash2
nb_read_unmap=subprocess.check_output(bash2,shell=True)
# Compute % of unmapped reads on total number of reads
pourc_unmap_read=float(nb_read_unmap)/float(tot_nb_read)*100.0
print "   - There are "+str(int(tot_nb_read))+" reads in file : "+BAM_file
print "   - "+str(pourc_unmap_read)+" % of reads have been unmapped!!!"
