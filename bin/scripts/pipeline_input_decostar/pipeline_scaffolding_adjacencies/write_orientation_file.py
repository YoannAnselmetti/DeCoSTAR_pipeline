#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
###   Goal:
###      Write in sequencing library orientation file for genome scaffolding with BESST
###
###   INPUT:
###      1- Directory containing BAM stats files obtained with Picard tools
###         (data/DATA_SEQ/STATS/MAPPING/Bowtie2_k100/TRIM/)
###      2- Species name
###         (Gallus_gallus)
###      3- Library orientation file
###         (data/DATA_SEQ/orientation_libraries)
###
###   OUTPUT:
###      - Sequencing library orientation file with informatins on species given as input
###
###   Name: write_orientation_file.py             Author: Yoann Anselmetti
###   Creation date: 2017/11/09               Last modification: 2020/11/05
###

from sys import argv
from re import search
from os import close, path, makedirs, listdir, stat
from datetime import datetime
from collections import namedtuple   #New in version 2.6

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise

def read_BAM_stats_file(BAMstats_file):
    dict_SRX_BAMstats=dict()
    stats_file=open(BAMstats_file,"r")
    for line in stats_file:
        r=search("^([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t ]*)[\t ]+([^\t\n ]*)\n$",line)
        if r:
            print(line)
            median_size=r.group(1)
            median_abs_dev=r.group(2)
            mean=r.group(5)
            sd=r.group(6)
            reads_nb=r.group(7)
            orientation=r.group(8)

            stats=MEAN_SD(mean.replace(",","."),sd.replace(",","."))
            dict_SRX_BAMstats[orientation]=stats

    stats_file.close()
    return dict_SRX_BAMstats



################
###   MAIN   ###
################
if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   STATS_DIR=argv[1]
   species=argv[2]
   orientation_file=argv[3]

   MEAN_SD=namedtuple("MEAN_SD",["mean","sd"])
   STATS=namedtuple("STATS",["ori","mean","sd"])
   

   dict_SRX=dict()
   for SRX in listdir(STATS_DIR+"/"+species):
      stats_file=STATS_DIR+"/"+species+"/"+SRX+"/stats_insert_size_50000"
      dict_BAMstats_file=read_BAM_stats_file(stats_file)
      if "RF" in dict_BAMstats_file:
         mean=dict_BAMstats_file["RF"].mean
         sd=dict_BAMstats_file["RF"].sd
         stats=STATS("rf",mean,sd)
         dict_SRX[SRX]=stats
      elif "FR" in dict_BAMstats_file:
         mean=dict_BAMstats_file["FR"].mean
         sd=dict_BAMstats_file["FR"].sd
         stats=STATS("fr",mean,sd)
         dict_SRX[SRX]=stats
      else:
         exit("ERROR, the orientation of the insert size library "+SRX+" should be \"FR\" or \"RF\"!!!")


   output_file = open(orientation_file,"a")
   if stat(orientation_file).st_size == 0:
      output_file.write("#species\tSRX_ID\torientation\tmean_insert_size\tSD_insert_size\n")
   for SRX in sorted(dict_SRX):
      ori=dict_SRX[SRX].ori
      mean=dict_SRX[SRX].mean
      sd=dict_SRX[SRX].sd
      output_file.write(species+"\t"+SRX+"\t"+ori+"\t"+mean+"\t"+sd+"\n")
   output_file.close()


   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
