#! /usr/bin/env python
# -*- coding: utf-8 -*-
###                                                                            
###   Goal:                                                                    
###      Draw species tree with ETE toolkit with annotation to nodes to add evolution or assembly statistics (#gene, #rearrangements, #scaffolds, ...)
###                                                                            
###   INPUT:                                                                   
###      1- Species tree file
###         (data/INPUT_DATA/Anopheles_species_tree_X_topology.nwk)
###         (data/INPUT_DATA/Anopheles_species_tree_WG_topology.nwk)
###      2- SCJ stats file
###         (results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_0.1_M1_scj_stats)
###         (results/decostar/Xtopo_RAW/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1_0.1_M1_scj_stats)
###         (results/decostar/Xtopo-scaff/DeCoSTAR_Anopheles_Xtopo-scaff_Boltz_0.1_0.1_M1_scj_stats)
###         (results/decostar/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1_0.1_M1_scj_stats)
###      3- New adjacencies file
###         (results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_0.1_M1_kept)
###         (results/decostar/Xtopo_RAW/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1_0.1_M1_kept)
###         (results/decostar/Xtopo-scaff/DeCoSTAR_Anopheles_Xtopo-scaff_Boltz_0.1_0.1_M1_kept)
###         (results/decostar/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1_0.1_M1_kept)
###      4- OUTPUT graph file
###         (figures/ETE_species_trees/Xtopo+scaff/ETE_species_tree_Xtopo+scaff.pdf)
###         (figures/ETE_species_trees/Xtopo_RAW/ETE_species_tree_Xtopo_RAW.pdf)
###         (figures/ETE_species_trees/Xtopo-scaff/ETE_species_tree_Xtopo-scaff.pdf)
###         (figures/ETE_species_trees/WGtopo+scaff/ETE_species_tree_WGtopo+scaff.pdf)
###                                                                            
###   OUTPUT:                                                                  
###      - Create annotation gene file for ARt-DeCo_seq instance creation      
###                                                                            
###   Name: draw_tree_ETE.py            Author: Yoann Anselmetti     
###   Creation date: 2015/11/11         Last modification: 2017/11/03
###                                                                            

from sys import argv, exit
from re import search, match
from os import close, path, makedirs
from collections import namedtuple   #New in version 2.6
import subprocess
import errno
from datetime import datetime
import random
from ete3 import Tree, TreeStyle, NodeStyle, faces, TextFace, AttrFace, CircleFace, PieChartFace



def mkdir_p(dir_path):
    output_dir=path.dirname(path.abspath(dir_path))
    try:
        makedirs(output_dir)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and path.isdir(output_dir):
            pass
        else:
            raise


def parse_new_adj_file(new_adj_file):
    dict_spe_degree=dict()
    dict_adj=dict()
    file=open(new_adj_file,"r")
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
                # Increment adjacencies number of current species
                if not spe_ID in dict_adj:
                    dict_adj[spe_ID]=0
                dict_adj[spe_ID]+=1

                # Increment gene degree of gene1 and gene2 
                if not spe_ID in dict_spe_degree:
                    dict_spe_degree[spe_ID]=dict()
                if not gene1 in dict_spe_degree[spe_ID]:
                    dict_spe_degree[spe_ID][gene1]=0
                if not gene2 in dict_spe_degree[spe_ID]:
                    dict_spe_degree[spe_ID][gene2]=0
                dict_spe_degree[spe_ID][gene1]+=1
                dict_spe_degree[spe_ID][gene2]+=1

    file.close()


    dict_degree=dict()
    for spe in dict_spe_degree:
        if not spe in dict_degree:
            dict_degree[spe]=dict()
        for gene in dict_spe_degree[spe]:
            adj_degree=dict_spe_degree[spe][gene]
            # To sum all genes with adj degree>=3
            if adj_degree>=3:
                adj_degree=3
            if not adj_degree in dict_degree[spe]:
                dict_degree[spe][adj_degree]=0

            dict_degree[spe][adj_degree]+=1

    # for spe in dict_degree:
    #     print "For species "+spe
        # for degree in dict_degree[spe]:
        #     gene_nb=dict_degree[spe][degree]
        #     print "\t- "+str(gene_nb)+" genes have an adjacency degree equivalent to "+str(degree)

    return dict_degree,dict_adj







def parse_stats_file(stats_file,factor,dict_degree,dict_adj):
    dict_stats=dict()
    file=open(stats_file,"r")
    for line in file:
        r=search("^([AE][0-9]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)[\t ]*([^\t ]*)\n$",line)
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
            sclX=r.group(10)
            scjA=r.group(11)
            dup=r.group(12)

            if spe!="#species":
                spe_ID=int(spe[1:])
                if spe_ID in dict_degree:
                    dict_degree[spe_ID][0]=int(gene_alone)
                else:
                    exit("ERROR, species "+str(spe_ID)+" should be present in dict_degree[]!!!")
                adj_nb=dict_adj[spe_ID]
                scj_per_adj=float(scj)/float(adj_nb)
                stats=STATS(float(float(gene)/factor),float(scaffolds),float(scj_per_adj))
                dict_stats[spe_ID]=stats
                # print str(spe_ID)+":\n\t",
                # print stats
    file.close()
    return dict_stats,dict_degree



def modif_species_tree(tree,dict_stats,dict_degree):
    i=0
    for node in tree.traverse(strategy='postorder'):
        # print node.name+":"
        stats=dict_stats[i]
        node.dist=stats.scj*1000.0

        adj0=float(dict_degree[i][0])
        adj1=float(dict_degree[i][1])
        adj2=float(dict_degree[i][2])
        if 3 in dict_degree[i]:
            adj3=float(dict_degree[i][3])
        else:
            adj3=0.0

        gene_tot=float(adj0+adj1+adj2+adj3)
        adj0=adj0/gene_tot*100.0
        adj1=adj1/gene_tot*100.0
        adj2=adj2/gene_tot*100.0
        adj3=adj3/gene_tot*100.0

        node.add_features(id=i, size=stats.gene, scaff=stats.scaff, adj0=adj0, adj1=adj1, adj2=adj2, adj3=adj3)
        # node.add_features(id=i, size=stats.gene, scaff=stats.scaff)

        i+=1
    return tree




def draw_chartpie(node,position,opacity):
    # Creates a PieChart face
    # - diameter corresponds to gene number (node.size)
    # - section corresponds to adjacency degree of genes (adj0,adj1,adj2 & adj3)
    P = PieChartFace([node.adj0,node.adj1,node.adj2,node.adj3], node.size, node.size, colors=["DarkOrange","Cyan","Chartreuse","Red"], line_color=None)

    # Set opacity of the PieChart
    P.opacity = opacity
    # Position PieChart to the node
    faces.add_face_to_node(P, node, 0, position=position)


def draw_circle(node,position,opacity,color):
    # Creates a Circle face
    # - diameter corresponds to gene number (node.size)
    # - color corresponds to scaffolding completness (color)
    P = CircleFace(node.size, color=color)

    # # Let's make the sphere transparent
    P.opacity = opacity
    # # And place as a float face over the tree
    faces.add_face_to_node(P, node, 0, position=position)




def layout(node):
    node.img_style["size"]=0
    opacity=0.75
    scaff=node.scaff
    color=""

    if scaff<10:
        color="Chartreuse"
    elif scaff<50:
        color="LimeGreen"
    elif scaff<100:
        color="ForestGreen"
    elif scaff<250:
        color="Gold"
    elif scaff<500:
        color="DarkOrange"
    elif scaff<1000:
        color="Tomato"
    elif scaff<2500:
        color="OrangeRed"
    elif scaff<5000:
        color="Red"
    else:
        color="Black"



    if node.is_leaf():
        name=node.name.replace("_"," ")

        # Add node name to leaf nodes
        N = faces.TextFace(name, fsize=8, fgcolor="black")
        faces.add_face_to_node(N, node, 0, position="branch-right")

        # position values available: “branch-right”, “branch-top”, “branch-bottom”, “float”, “float-behind” and “aligned”
        position="branch-right"
        if bool_chartpie:
            draw_chartpie(node,position,opacity)
        else:
            draw_circle(node,position,opacity)

    else:
        # position values available: “branch-right”, “branch-top”, “branch-bottom”, “float”, “float-behind” and “aligned”
        position="float"
        if bool_chartpie:
            draw_chartpie(node,position,opacity)
        else:
            draw_circle(node,position,opacity,color)



    node.img_style["size"]=0





def draw_tree(species_tree,dict_stats,factor,dict_degree,bool_chartpie):
    t_init = Tree(species_tree)
    t_modif = modif_species_tree(t_init,dict_stats,dict_degree)

    # Create an empty TreeStyle
    ts = TreeStyle()

    # Set our custom layout function
    ts.layout_fn = layout

    # Draw a tree
    ts.mode = "r"

    # We will add node names manually
    ts.show_leaf_name = False

    # Show branch data
    ts.show_branch_length = False
    # ts.show_branch_support = True

    ts.scale=2

    return t_modif, ts




################
###   MAIN   ###
################
if __name__ == '__main__':

    start_time = datetime.now()

    species_tree=argv[1]
    stats_file=argv[2]
    newAdj_file=argv[3]
    output_graph=argv[4]

    mkdir_p(output_graph)


    DEGREE=namedtuple("DEGREE",["adj0","adj1","adj2","adj3"])
    STATS=namedtuple("STATS",["gene","scaff","scj"])


    factor=500.0
    dict_degree,dict_adj=parse_new_adj_file(newAdj_file)
    dict_stats,dict_degree=parse_stats_file(stats_file,factor,dict_degree,dict_adj)


    bool_chartpie=True
    t, ts = draw_tree(species_tree,dict_stats,factor,dict_degree,bool_chartpie)
    t.render(output_graph, w=2500, dpi=600, tree_style=ts)
    # t.show(tree_style=ts)


    # Get and print execution time of the script
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
