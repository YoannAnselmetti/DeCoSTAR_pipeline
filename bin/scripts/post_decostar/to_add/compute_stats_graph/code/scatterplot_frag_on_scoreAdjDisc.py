#! /usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Create scatter plot of genome fragmentation on cumulative sum of a posteriori scores of adjacencies discrarded during genome linearization
###   INPUT:
###      1- File of adjacencies discarded during linearization
###         (results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_0.1_M1_disc)
###         (results/decostar/Xtopo_RAW/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1_0.1_M1_disc)
###         (results/decostar/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1_0.1_M1_disc)
###      2- SCJ stats file (contains #scaffolds / species)
###         (results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_0.1_M1_scj_stats)
###         (results/decostar/Xtopo_RAW/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1_0.1_M1_scj_stats)
###         (results/decostar/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1_0.1_M1_scj_stats)
###      3- CTG file (contains #scaffolds / species)
###         (data/data_DeCoSTAR/CTG_file)
###      4- Chromosome file (contains #scaffolds / species)
###         (data/INPUT_DATA/18Anopheles_species)
###      5- Species file to get correspondance between EXTANT species number ID and name
###         (results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1.species.txt)
###      6- OUTPUT directory where graph will be stored
###         (figures/scatterplot/Xtopo/scatterplot_Xtopo)
###         (figures/scatterplot/Xtopo_RAW/scatterplot_Xtopo_RAW)
###         (figures/scatterplot/WGtopo/scatterplot_WGtopo)
###
###   OUTPUT:   (RUN in ~??min)
###      - OUTPUT file with new adj with CTG info proposed by ARt-DeCo_seq (Format: species ctg1 ctg2 ori_ctg1 ori_ctg2 score_BESST)
###
###   Name: scatterplot_frag_on_scoreAdjDisc.py      Author: Yoann Anselmetti
###   Creation date: 2017/03/07                      Last modification: 2017/11/03
###

from sys import argv, stdout
from re import search
from os import close, listdir, path, makedirs
from datetime import datetime
from collections import namedtuple   #New in version 2.6
import errno
import subprocess

import itertools
import matplotlib
matplotlib.use('Agg')
import numpy as np 
from matplotlib import pyplot as plt
from matplotlib.pyplot import cm

stdout.flush()



def mkdir_p(dir_path):
    output_dir=path.dirname(path.abspath(dir_path))
    try:
        makedirs(output_dir)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and path.isdir(output_dir):
            pass
        else:
            raise



def parse_disc_adj_file(disc_adj_file):
    dict_spe_SUMscore=dict()
    file=open(disc_adj_file,"r")
    for line in file:
        r=search("^([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)\n$",line)
        if r:
            spe=r.group(1)
            gene1=r.group(2)
            gene2=r.group(3)
            ori1=r.group(4)
            ori2=r.group(5)
            prior_score=r.group(6)
            post_score=r.group(7)

            if spe!="#species":
                spe_ID=int(spe)
                # Increment gene degree of gene1 and gene2 
                if not spe_ID in dict_spe_SUMscore:
                    dict_spe_SUMscore[spe_ID]=0
                dict_spe_SUMscore[spe_ID]+=float(post_score)

    file.close()

    return dict_spe_SUMscore


def read_species_file(species_file):
    dict_ID_name_spe=dict()
    file=open(species_file,"r")
    for line in file:
        r=search("^([^\t ]*)[\t ]*([^\t ]*)\n$",line)
        if r:
            spe_ID=r.group(1)
            spe_name=r.group(2)
            dict_ID_name_spe[int(spe_ID)]=spe_name
    file.close()
    return dict_ID_name_spe





def parse_res_file(stats_file):
    dict_spe_scaff=dict()
    file=open(stats_file,"r")
    for line in file:
        r=search("^([AE][0-9][0-9]?)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)\n$",line)
        if r:
            spe=r.group(1)
            gene=r.group(2)
            geneX=r.group(3)
            geneA=r.group(4)
            gene_alone=r.group(5)
            scaffolds=r.group(6)
            scaffX=r.group(7)
            scaffA=r.group(8)
            scj=r.group(9)
            scjX=r.group(10)
            scjA=r.group(11)
            dup=r.group(12)

            spe_ID=int(spe[1:])
            dict_spe_scaff[spe_ID]=int(scaffolds)

    file.close()

    return dict_spe_scaff



def parse_CTG_CHR_files(CTG_file,CHR_file):
    dict_spe_chr=dict()
    chr_file=open(CHR_file,"r")
    for line in chr_file:
        spe=line.split()[0]
        chrom= line.split()[1]
        dict_spe_chr[spe]=int(chrom)
    chr_file.close()
    
    dict_spe_nbCTG=dict()
    ctg_file=open(CTG_file,"r")
    for line in ctg_file:
        spe=line.split()[0]
        ctg=line.split()[1]
        if spe!="#species":
            if not spe in dict_spe_nbCTG:
                dict_spe_nbCTG[spe]=0
            dict_spe_nbCTG[spe]+=1
    ctg_file.close()

    dict_spe_scaff_chr=dict()
    for species in dict_spe_chr:
        chrom=dict_spe_chr[species]
        scaff=dict_spe_nbCTG[species]

        init=INIT(int(scaff),int(chrom))
        dict_spe_scaff_chr[species]=init

    return dict_spe_scaff_chr



def scatterplot_statsScaff_on_scoreSumDiscAdj(dict_to_plot,OUTPUT,zoom):
    n=len(dict_to_plot)
    color=cm.rainbow(np.linspace(0,1,n))
    i=0
    factor=2.0
    for spe in sorted(dict_to_plot):
        label=spe.replace('Anopheles_','A.')
        X=dict_to_plot[spe].X
        Y1=dict_to_plot[spe].Y1
        Y2=dict_to_plot[spe].Y2
        size=dict_to_plot[spe].size
        # plt.scatter(X, Y1, s=size*factor, c=color[i], alpha=0.75)
        plt.plot([X,X], [Y1,Y2], '-', c=color[i], lw=1)
        plt.scatter(X, Y2, s=size*factor, c=color[i], alpha=0.75, label=label)
        i+=1

    plt.grid(True)
    plt.title("Scatter plot")
    plt.xlabel("SUM of \"a posteriori\" scores of discarded ADJ\n(after genome linearization)\n")
    plt.ylabel("#Scaffolds")


    plt.legend(fontsize=8,loc='best')

    # scatter1 = plt.scatter((0,0), (0,0), s=1*factor, c="slategray", alpha=0.75, label="1%")
    # scatter2 = plt.scatter((0,0), (0,0), s=2*factor, c="slategray", alpha=0.75, label="5%")
    # scatter3 = plt.scatter((0,0), (0,0), s=5*factor, c="slategray", alpha=0.75, label="10%")
    # scatter4 = plt.scatter((0,0), (0,0), s=25*factor, c="slategray", alpha=0.75, label="25%")
    # scatter5 = plt.scatter((0,0), (0,0), s=50*factor, c="slategray", alpha=0.75, label="50%")
    # scatter6 = plt.scatter((0,0), (0,0), s=75*factor, c="slategray", alpha=0.75, label="75%")
    # scatter7 = plt.scatter((0,0), (0,0), s=100*factor, c="slategray", alpha=0.75, label="100%")
    # legend1 = plt.legend((scatter1,scatter2,scatter3,scatter4,scatter5,scatter6,scatter7),('1%','5%','10%','25%','50%','100%'),numpoints=1, fontsize=8, loc='best')
    # plt.gca().add_artist(legend1)

    fig_name=""
    if zoom==1:
        # Default dispay of the graph
        fig_name=OUTPUT+".pdf"
    elif zoom==2:
        plt.xlim(xmin=-5,xmax=405)
        plt.ylim(ymin=0,ymax=5000)
        fig_name=OUTPUT+"_zoom.pdf"
    elif zoom==3:
        plt.xlim(xmin=-1,xmax=25)
        plt.ylim(ymin=0,ymax=2000)
        fig_name=OUTPUT+"_zoom_zoom.pdf"

    plt.tight_layout()
    plt.savefig(fig_name,format='pdf')
    plt.cla()




################
###   MAIN   ###
################
if __name__ == '__main__':

    start_time = datetime.now()

    disc_adj_file=argv[1]
    res_file=argv[2]
    CTG_file=argv[3]
    CHR_file=argv[4]
    species_file=argv[5]
    output_graph=argv[6]

    mkdir_p(output_graph)

    INIT=namedtuple("INIT",["scaff","chr"])
    PLOT=namedtuple("PLOT",["X","Y1","Y2","size"])

    dict_spe_SUMscore=parse_disc_adj_file(disc_adj_file)
    dict_spe_init=parse_CTG_CHR_files(CTG_file,CHR_file)
    dict_spe_scaff=parse_res_file(res_file)
    dict_ID_name_spe=read_species_file(species_file)


    for elem in dict_ID_name_spe:
        print str(elem)+" -> "+dict_ID_name_spe[elem]


    dict_spe_scaffANDsumDisc=dict()
    for spe in sorted(dict_spe_scaff):
        if not spe in dict_spe_SUMscore:
            dict_spe_SUMscore[spe]=0.0

        res_scaff=dict_spe_scaff[spe]
        SUMscore=dict_spe_SUMscore[spe]
        if spe in dict_ID_name_spe:
            spe_name=dict_ID_name_spe[spe]

            init_stats=dict_spe_init[spe_name]
            init_scaff=float(init_stats.scaff)
            chrom=float(init_stats.chr)
            scaff_improv=(init_scaff-float(res_scaff))/(init_scaff-chrom)*100.0

            if scaff_improv>0.0:
                stats=PLOT(SUMscore,init_scaff,res_scaff,scaff_improv)
                dict_spe_scaffANDsumDisc[spe_name]=stats

            print "For EXTANT species "+spe_name+":\n\t"+str(res_scaff)+" scaffolds ("+str(scaff_improv)+"% of scaffolding improvment) AND a cumulative SUM of discarded ADJ during genome linearization of "+str(SUMscore)
        else:
            print "For ANCESTRAL species "+str(spe)+":\n\t"+str(res_scaff)+" scaffolds AND a cumulative SUM of discarded ADJ during genome linearization of "+str(SUMscore)


    scatterplot_statsScaff_on_scoreSumDiscAdj(dict_spe_scaffANDsumDisc,output_graph,1)
    scatterplot_statsScaff_on_scoreSumDiscAdj(dict_spe_scaffANDsumDisc,output_graph,2)
    scatterplot_statsScaff_on_scoreSumDiscAdj(dict_spe_scaffANDsumDisc,output_graph,3)


    end_time = datetime.now()
    print('\nDuration: {}'.format(end_time - start_time))
