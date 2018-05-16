__author__="Cedric Chauve"
__email__="cedric.chauve@sfu.ca"
__date__="October 2, 2016"

import sys
import os
import math
from operator import itemgetter
import networkx as nx
from random import shuffle

USAGE="python linearize_genomes.py <adjacencies_file> [<species> <threshold> <output_prefix> <algorithm> <observed> <epsilon>]\n"
USAGE+="  adjacencies_file: list of adjacencies in the DeCo* format\n"
USAGE+="  species:          list of species whose genome will be linearized\n"
USAGE+="                    either a species id, or ALL for all species, or a file with one species per line (first field is th species id, comments can follow)\n"
USAGE+="  threshold:        floating number in [0,1], all adjacencies with a posterior score < threshold are discarded\n"
USAGE+="  output_prefix:    path+prefix of output files\n"
USAGE+="  algorithm:        M1/M2/F1/F2 (see below for an explanation)\n"
USAGE+="  observed:         floating number in [0,1] representing the score of observed adjacencies\n"
USAGE+="  epsilon:          floating number in [0,1] used to decrease the score and avoid an infinite  value for the log(1-score) when score=1\n"
USAGE+="\n"
USAGE+="This script takes a list of adjacencies and a set of species and linearize the genome of these species by discarding adjacencies to keep a subset that is compatible"
USAGE+="with a set of linear genome segments\n"
USAGE+="DeCo* adjacency format: species_id gene1 gene2 orientation_gene1(+/-) orientation_gene2(+/-) prior_score([0,1]) posterior_score([0,1])\n"
USAGE+="An adjacency is observed if it is seen in the initial contigs. The parameter <observed> indicates the score of such adjacencies (expected to be 1) and the"
USAGE+="algorithm assumes that all adjacencies with such a score are observed\n"
USAGE+="Algorithms F1/F2:\n"
USAGE+="  1. filtering out all adjencies of score < threshold\n"
USAGE+="  2. filtering out all adjacencies that link two gene extremities adjacent to more than 1 other extremity (overloaded gene extremities)\n"
USAGE+="  3. filtering out all adjacencies that are incident to one overloaded gene extremity\n" 
USAGE+="  4. for each resulting circular scaffold, discarding an adjacency of minimum weight to make it linear\n"
USAGE+="  Variant F1: any adjacency can be discarded.\n"
USAGE+="  Variant F2: observed adjacencies can not be discarded.\n"
USAGE+="Algorithms M1/M2:\n"
USAGE+="  1. filtering out all adjencies of score < threshold\n"
USAGE+="  2. A maximum weight matching (MWM) is applied to extract a subset of adjacencies\n"
USAGE+="  4. for each resulting circular scaffold, discarding an adjacency of minimum weight to make it linear\n"
USAGE+="  Variant M1: any adjacency can be discarded.\n"
USAGE+="  Variant M2: observed adjacencies can not be discarded.\n"
USAGE+="\n"
USAGE+="The output is composed of \n"
USAGE+="  a conflict-free adjacencies file\n"
USAGE+="  a scaffold file, where genes are ordered along scaffolds post-linearization"
USAGE+="  a discarded adjacencies file where each discarded adjacency is labeled by the stage where it has been filtered out \n"
USAGE+="  (1 to 4 for algorithms F1/F2 and 1, MWM or 4 for algorithms M1/M2)\n"


# --------------------------------------------------------------------------------------------
MAX_WEIGHT=0     # Maximum weight of an adjacency
GENE_LIST=[]     # List of all considered genes
NB_GENES=0       # Number of genes
NB_ADJS=0        # Number of adjacencies
OBSERVED=1.0     # Default score for observed adjacencies (true value to be read as a parameter)
EPSILON=0.000001 # Default value of the score correction to avoid infinite log-score (true value to be read as a parameter)

# DATA STRUCTURE FOR ADJACENCY: (species,gene_name,extremity('h','t'),gene_name,extremity,prior,weight,posterior)
# where weight is -log(1-posterior)
def adj_species(adj):
    return(adj[0])
def adj_gene1(adj):
    return(adj[1])
def adj_gene2(adj):
    return(adj[3])
def adj_ext1(adj):
    return(adj[2])
def adj_ext2(adj):
    return(adj[4])
def adj_weight(adj):
    return(adj[5])
def adj_prior(adj):
    return(adj[6])
def adj_posterior(adj):
    return(adj[7])
def adj_observed(adj): # True/False
    return(adj[6]==OBSERVED)
def adj_genes(adj):
    g1=adj_gene1(adj)
    g2=adj_gene2(adj)
    return((g1,g2))
def adj_genes_exts(adj):
    (g1,ext1)=(adj_gene1(adj),adj_ext1(adj))
    (g2,ext2)=(adj_gene2(adj),adj_ext2(adj))
    return(((g1,ext1),(g2,ext2)))

# Translating the orientation of the contigs into the corresponding gene extremity involved in an adjacency
GENE1_ORIENT2EXT={'-':'t','+':'h','?':'?'}
GENE2_ORIENT2EXT={'-':'h','+':'t','?':'?'}
GENE1_EXT2ORIENT={'t':'-','h':'+','?':'?'}
GENE2_EXT2ORIENT={'h':'-','t':'+','?':'?'}

def write_adj(adj):
    return(adj_species(adj)+" "+adj_gene1(adj)+" "+adj_gene2(adj)+" "+GENE1_EXT2ORIENT[adj_ext1(adj)]+" "+GENE2_EXT2ORIENT[adj_ext2(adj)]+" "+str(adj_prior(adj))+" "+str(adj_posterior(adj)))

def complement_sign(sign):
    if sign=='+':
        return('-')
    elif sign=='-':
        return('+')
    elif sign=='?':
        return('?')

def complement_ext(ext):
    if ext=='t':
        return('h')
    elif ext=='h':
        return('t')
    elif ext=='?':
        return('?')

# --------------------------------------------------------------------------------------------                                                                                
# Reading the adjacencies in the DeCo* format 
# file_stream = opened file
# output = adj_list, list of adjacencies
# updates the GENES_LIST and MAX_WEIGHT global variables
def read_adjacencies(file_stream,species):
    global MAX_WEIGHT
    global GENES_LIST
    global NB_GENES
    global NB_ADJS
    genes_list_aux=[] # List of genes with duplicates
    res_adj_list=[] # List of adjacencies, returned
    for adj in file_stream:
        if adj[0]!="#": # "#" is used for comments
            adj1=adj.rstrip().split()
            if adj1[0]==species:
                assert(float(adj1[5])<=1.0 and float(adj1[6])<=1.0), "ERROR: score greater than 1" 
                if float(adj1[6]) > 1.0:
                    w=-math.log1p(-(1.0-EPSILON))
                else:
                    w=-math.log1p(-(float(adj1[6])-EPSILON))
                genes_list_aux.append(adj1[1])
                genes_list_aux.append(adj1[2])
                res_adj_list.append((adj1[0],adj1[1],GENE1_ORIENT2EXT[adj1[3]],adj1[2],GENE2_ORIENT2EXT[adj1[4]],w,float(adj1[5]),float(adj1[6])))
                NB_ADJS+=1
                if (w>MAX_WEIGHT):
                    MAX_WEIGHT=w
    GENES_LIST=list(set(genes_list_aux))
    GENES_LIST.sort()
    NB_GENES=len(GENES_LIST)
    return(res_adj_list)

# Computes the degree of each gene of a given list given a set of adjacencies
def compute_gene_degrees(adj_list):
    gene_degrees={}   # number of neighbours of each gene
    for g in GENES_LIST:
        gene_degrees[g]=0
    for adj in adj_list:
        (g1,g2)=adj_genes(adj)
        gene_degrees[g1]+=1
        gene_degrees[g2]+=1
    return(gene_degrees)

# Computes the degree of each gene extremity of a given list given a set of adjacencies
def compute_ext_degrees(adj_list):
    ext_degrees={}   # number of neighbours of each gene
    for g in GENES_LIST:
        ext_degrees[(g,'h')]=0
        ext_degrees[(g,'t')]=0
    for adj in adj_list:
        (g1ext1,g2ext2)=adj_genes_exts(adj)
        ext_degrees[g1ext1]+=1
        ext_degrees[g2ext2]+=1
    return(ext_degrees)

# Compute the list of neighbours of a gene extremities
# A neighbours is a pair (gene_name,extremity)
def compute_ext_ngbs(adj_list):
    neighbours={} # List of all genextremities linked to a given gene extremity
    for g in GENES_LIST:
        neighbours[(g,'h')]=[]
        neighbours[(g,'t')]=[]
    for adj in adj_list:
        (g1ext1,g2ext2)=adj_genes_exts(adj)
        neighbours[g1ext1].append(g2ext2)
        neighbours[g2ext2].append(g1ext1)
    return(neighbours)

# --------------------------------------------------------------------------------------------
# Filtering functions: each returns a pair of lists of adjacencies (kept_adjacencies, discarded_adjacencies)

# Filtering out all unobserved adjacencies
def filter_unobserved(adj_list):
    res_kept_adjs=[]
    res_disc_adjs=[]
    for adj in adj_list:
        if adj_observed(adj):
            res_kept_adjs.append(adj)
        else:
            res_disc_adjs.append(adj)
    res_kept_adjs.sort()
    res_disc_adjs.sort()
    return((res_kept_adjs,res_disc_adjs))

# Filtering out all adjacencies of score  < threshold
# Assumption: observed extant adjacencies have a score of 1 and will never be filtered
# param 1 = list of adjacencies
# param 2 = filtering threshold
def filter_1(adj_list,threshold):
    res_kept_adjs=[]
    res_disc_adjs=[]
    for adj in adj_list:
        if adj_weight(adj)<threshold:
            res_disc_adjs.append(adj)
        else:
            res_kept_adjs.append(adj)
    res_kept_adjs.sort()
    res_disc_adjs.sort()
    return((res_kept_adjs,res_disc_adjs))

# Filtering stegaes 2 and 3 implemented into a single fucntion
# param 1 = list of adjacencies
# param 2 = keep observed adjacencies (True/False)
# param 3 = '2': Filtering out all adjacencies involved in two conflicts
# param 3 = '3': Filtering out all adjacencies involved in one conflict
def filter_23(adj_list,keep_obs,p3):    
    res_kept_adjs=[]    # Kept adjacencies
    res_disc_adjs=[]    # Discarded adjacencies
    degrees={}          # Degree of each gene extremity
    # 1. Computing the degree of each gene extremity
    ext_degrees=compute_ext_degrees(adj_list)
    # 2. Filtering out
    for adj in adj_list:
        (g1ext1,g2ext2)=adj_genes_exts(adj)
        if p3==2:
            test_disc=ext_degrees[g1ext1]>1 and ext_degrees[g2ext2]>1 and (keep_obs==False or adj_observed(adj)==False)
        elif p3==3:
            test_disc=(ext_degrees[g1ext1]>1 or ext_degrees[g2ext2]>1) and (keep_obs==False or adj_observed(adj)==False)
        if test_disc:
            res_disc_adjs.append(adj)
        else:
            res_kept_adjs.append(adj)
    res_kept_adjs.sort()
    res_disc_adjs.sort()
    return((res_kept_adjs,res_disc_adjs))


# --------------------------------------------------------------------------------------------
# Computing the scaffolds, as connected components of the set of conflict-free adjacencies
# seen as a graph
# Recursive graph exploration
# param 1 = starting gene (assumed to have not been assigned to a scaffold so far)
# param 2 = ID of the current scaffold/connected compoenent
# param 3 = array of neighbours of genes (indexed by gene names)
# param 4 = array indicating the scaffold/connected component associated to each gene
def find_cc_rec(gene,cc_id,ngbs,cc):
    cc[gene]=cc_id
    for (g,ext) in ngbs[(gene,'h')]:
        if cc[g]==0:
            find_cc_rec(g,cc_id,ngbs,cc)
    for (g,ext) in ngbs[(gene,'t')]:
        if cc[g]==0:
            find_cc_rec(g,cc_id,ngbs,cc)

# Initialisation and main function
# Returns the number of scaffolds and an array indexed by genes that indicate in which scaffold each belongs to
# param 1 = adjacencies list
# Returns (nb_cc,cc), nb_cc=number of connected components, cc=array indexed by genes, with ID of the connected comp.
def compute_scaffolds_circ(adj_list):
    sys.setrecursionlimit(max(1000,10+len(adj_list)))
    connected_components={} # Each gene receives the unique ID of its connected component
    # 1. Computing a list of neighbour for each gene extremity
    neighbours=compute_ext_ngbs(adj_list)
    for g in GENES_LIST:
        connected_components[g]=0
    # 2. Recursive identification of connected components
    cc_id=1
    for g in GENES_LIST:
        if connected_components[g]==0:
            find_cc_rec(g,cc_id,neighbours,connected_components)
            cc_id+=1
    return((cc_id-1,connected_components))

# Linearizing circular scaffolds by removing a min-weight adjacency per circular scaffold
# param 1 = adjacencies list
# param 2 = array indicating the scaffold/connected component associated to each gene
# param 3 = number of scaffolds (including gthe ones reduced to a single gene)
# param 4 = genes copy numbers 
def filter_4(adj_list,cc,nb_ccs):
    res_kept_adjs=[] # Kept adjacencies
    res_disc_adjs=[] # Discarded adjacencies
    # 1. Recording the degree of each gene
    gene_degrees=compute_gene_degrees(adj_list)
    # 2. Recording the scaffolds with two genes of degree 1
    nb_ext_ccs={} # Number of extremity of scaffolds (0 or 2)
    for c in range(1,nb_ccs+1):
        nb_ext_ccs[c]=0
    for adj in adj_list:
        (g1,g2)=adj_genes(adj)
        if gene_degrees[g1]==1:
            nb_ext_ccs[cc[g1]]+=1
        if gene_degrees[g2]==1:
            nb_ext_ccs[cc[g2]]+=1
    # 3. Checking that each scaffold has indeed value 0 or 4 
    for c in range(1,nb_ccs+1):
        if nb_ext_ccs[c]!=0 and nb_ext_ccs[c]!=2:
            print("ERROR TYPE 1\n")
    # 4. Computing the min weight of each scaffolds
    min_weights={} # Min weight of adjacencies of each scaffold, indexed by scaffolds ID
    for c in range(1,nb_ccs+1):
        min_weights[c]=MAX_WEIGHT+1
    for adj in adj_list:
        g1=adj_gene1(adj)
        if adj_weight(adj)<min_weights[cc[g1]]:
            min_weights[cc[g1]]=adj_weight(adj)    
    # 4. Filtering
    linearized={} # Status of each scaffold, indexed by scaffolds ID
    for c in range(1,nb_ccs+1):
        if nb_ext_ccs[c]==0: # Circular scaffold
            linearized[c]=False
        else:
            linearized[c]=True
    for adj in adj_list:
        g1=adj_gene1(adj)
        c1=cc[g1]
        if linearized[c1]==False and adj_weight(adj)==min_weights[c1]:
            res_disc_adjs.append(adj)
            linearized[c1]=True
        else:
            res_kept_adjs.append(adj)
    res_kept_adjs.sort()
    res_disc_adjs.sort()
    return((res_kept_adjs,res_disc_adjs))

# --------------------------------------------------------------------------------------------
# Ordering genes along scaffolds
# Assume that the genes have all copy number 1 and all scaffolds are linear
# param 1 = adjacencies list
# param 2 = array indicating the scaffold/connected component associated to each gene
# param 3 = number of scaffolds (including gthe ones reduced to a single gene)                 
def order_genes(adj_list,cc,nb_ccs):
    gene_degrees=compute_gene_degrees(adj_list)
    # Finding a starting gene for each scaffold
    cc_start={} # Record one starting point for each scaffold
    cc_list=sorted(list(cc.keys()))
    for g in cc_list:
        cc_start[cc[g]]=(g,'t') # Used to handle scaffolds with one gene
    for adj in adj_list:
        ((g1,ext1),(g2,ext2))=adj_genes_exts(adj)
        # If an adjacency has a gene of degree 1, this gene starts the scaffold.
        if gene_degrees[g1]==1:
            cc_start[cc[g1]]=(g1,complement_ext(ext1))
        elif gene_degrees[g2]==1:
            cc_start[cc[g2]]=(g2,complement_ext(ext2))
    # Recording the neighbours of gene extremities
    neighbours=compute_ext_ngbs(adj_list)
    # Creating the scaffolds
    scaffolds={} # Array of scaffolds, indexed by scaffold ID
    for c in range(1,nb_ccs+1):
        current_gene=cc_start[c] # First gene of the scaffold
        scaffolds[c]=[] # A scaffold is an ordered list of oriented genes
        scaffolds[c].append(current_gene) # We add the first gene of the scaffold
        current_gene=(current_gene[0],complement_ext(current_gene[1])) # Current gene: other extremity
        while len(neighbours[current_gene[0],current_gene[1]])==1: # We stop when a gene has no neighbour
            current_gene=neighbours[current_gene[0],current_gene[1]][0]
            scaffolds[c].append(current_gene)
            current_gene=(current_gene[0],complement_ext(current_gene[1]))
    return(scaffolds)

def write_scaffolds(scaffolds,species):
    res_str=""
    scaffolds_list=sorted(list(scaffolds.keys()))
    for c in scaffolds_list:
        scf=scaffolds[c]
        for (g,ext) in scf:
            res_str+=species+" "+str(c)+" "+g+" "+GENE2_EXT2ORIENT[ext]+"\n"
    return(res_str)

# --------------------------------------------------------------------------------------------
# Maximum Weight Matching: Creating a weighted graph from an adjacencies list
# param1 = list of adjacencies
# param2 = keep observed adjacencies (True/False)

def MWM_create_graph(adj_list,keep_observed):
    G = nx.Graph()
    for g in GENES_LIST:
        G.add_node((g,'t'))
        G.add_node((g,'h'))
    edges_list=[]
    for adj in adj_list:
        (g1ext1,g2ext2)=adj_genes_exts(adj)
        if keep_observed==False or adj_observed(adj)==False:
            edges_list.append((g1ext1,g2ext2,adj_weight(adj)))
        else:
            edges_list.append((g1ext1,g2ext2,float(NB_ADJS)*MAX_WEIGHT+1.0))
    #shuffle(edges_list)
    G.add_weighted_edges_from(edges_list)
    return(G)

# Returning the kept and discarded edges using a MWM algorithm
def MWM_filter(G,adj_list):
    res_kept_adjs=[] # Returned list of kept adjacencies
    res_disc_adjs=[] # Returned list of discarded adjacencies
    ADJ_KEPT_MWM={}  # Indexed by adjacencies: records if kept or not by MWM
    ADJ_INFO={}      # Adj associuated to a pair of gene extremities
    # Initializing the data structures
    for adj in adj_list:
        (g1ext1,g2ext2)=adj_genes_exts(adj)
        ADJ_KEPT_MWM[adj]=False
        ADJ_KEPT_MWM[adj]=False
        ADJ_INFO[g1ext1,g2ext2]=adj
        ADJ_INFO[g2ext2,g1ext1]=adj
    # Maximum Weight Matching algorithm
    M=nx.max_weight_matching(G)
    # Translating the matching as returned by the MWM fucntion
    M_LIST=sorted(list(M.keys()))
    for m1 in M_LIST:
        m2=M[m1]
        ADJ_KEPT_MWM[ADJ_INFO[(m1,m2)]]=True
    # Creating the output
    for adj in adj_list:
        if ADJ_KEPT_MWM[adj]==True:
            res_kept_adjs.append(adj)
        else:
            res_disc_adjs.append(adj)
    res_kept_adjs.sort()
    res_disc_adjs.sort()
    return((res_kept_adjs,res_disc_adjs))

# --------------------------------------------------------------------------------------------
# Function processing a single genome

def process_genome(SPECIES,ADJ_FILE_STREAM,THRESHOLD,ALGORITHM,OUTPUT_PREFIX):
    print("Species "+SPECIES+" --------------------------")
    # Reading adjacencies
    print("  Reading adjacencies")
    ADJ_LIST=read_adjacencies(ADJ_FILE_STREAM,SPECIES)
    print("  Number of adjacencies: "+str(len(ADJ_LIST)))

    # Observed adjacencies
    print("  Observed adjacencies")
    (ADJ_OBSERVED,ADJ_UNOBSERVED)=filter_unobserved(ADJ_LIST)
    (NB_SCFS_OBS,SCFS_OBS)=compute_scaffolds_circ(ADJ_OBSERVED)
    OBS_SCAFFOLDS=order_genes(ADJ_OBSERVED,SCFS_OBS,NB_SCFS_OBS)

    # Filtering adjacencies
    print("  Filter 1")
    (ADJ_KEPT_1,ADJ_DISC_1)=filter_1(ADJ_LIST,THRESHOLD)
    print("    Number of adjacencies: "+str(len(ADJ_KEPT_1)))
    if ALGORITHM[0]=='F':
        print("  Filter 2")
        (ADJ_KEPT_2,ADJ_DISC_2)=filter_23(ADJ_KEPT_1,ALGORITHM[1]=='2',2)
        print("    Number of adjacencies: "+str(len(ADJ_KEPT_2)))
        print("  Filter 3")
        (ADJ_KEPT_3,ADJ_DISC_3)=filter_23(ADJ_KEPT_2,ALGORITHM[1]=='2',3)
        print("    Number of adjacencies: "+str(len(ADJ_KEPT_3)))
    elif ALGORITHM[0]=='M':
        print("  Creating graph")
        G=MWM_create_graph(ADJ_KEPT_1,ALGORITHM[1]=='2')
        print("  Maximum Weight Matching")
        (ADJ_KEPT_3,ADJ_DISC_3)=MWM_filter(G,ADJ_KEPT_1)
        print("    Number of adjacencies: "+str(len(ADJ_KEPT_3)))
    print("  Creating scaffolds (including circular ones)")
    (NB_SCFS,SCFS)=compute_scaffolds_circ(ADJ_KEPT_3)
    print("  Filter 4")
    (ADJ_KEPT_4,ADJ_DISC_4)=filter_4(ADJ_KEPT_3,SCFS,NB_SCFS)
    print("    Number of adjacencies: "+str(len(ADJ_KEPT_4)))
    print("  Ordering genes along scaffolds")
    SCAFFOLDS=order_genes(ADJ_KEPT_4,SCFS,NB_SCFS)

    # Writing output
    print("  Writing output")
    print("    Observed scaffolds")
    output_obs_scfs.write(write_scaffolds(OBS_SCAFFOLDS,SPECIES))
    print("    Improved scaffolds")
    output_scfs.write(write_scaffolds(SCAFFOLDS,SPECIES))
    print("    Kept adjacencies")
    for adj in ADJ_KEPT_4:
        output_kept.write(write_adj(adj)+"\t"+str(SCFS[adj_gene1(adj)])+"\n")
    print("    Discarded adjacencies")
    for adj in ADJ_DISC_1:
        output_disc.write(write_adj(adj)+"\t1\n")
    if ALGORITHM[0]=='F':
        for adj in ADJ_DISC_2:
            output_disc.write(write_adj(adj)+"\t2\n")
        for adj in ADJ_DISC_3:
            output_disc.write(write_adj(adj)+"\t3\n")
    elif ALGORITHM[0]=='M':
        for adj in ADJ_DISC_3:
            output_disc.write(write_adj(adj)+"\tMWM\n")
    for adj in ADJ_DISC_4:
        output_disc.write(write_adj(adj)+"\t4\n")

# --------------------------------------------------------------------------------------------
# Checking the output result is consistent:
# - every kept adjacency is consistent with the scaffolds
# - every adjacency observed in a scaffold is in the kept adjacencies
# - no discarded adjacency can be added to the scaffolds

def check_results(kept_adjacencies_file, discarded_adjacencies_file, scaffolds_file):
    print("Checking consistency of the results---------------------")    
    # Reading the scaffolds and recording the positions of all genes
    scaffolds_stream=open(scaffolds_file,"r").readlines()
    genes_info={}  # Indexed by genes, fields = (species, scaffold_id, sign, last_in_its_scaffold, kept_adj)
    prev_sp=""
    for scf in scaffolds_stream:
        if scf[0]!="#":
            scf1=scf.rstrip().split()
            (sp,scfid,gene,sign)=(scf1[0],scf1[1],scf1[2],scf1[3])
            if sp!=prev_sp or scfid!=prev_scfid:
                pos=1
                if prev_sp!="": # Updxating the last field of the previous gene, which is last in its scaffold
                    (sp1,scfid1,pos1,sign1,last1,seen1)=genes_info[prev_gene]
                    genes_info[prev_gene]=(sp1,scfid1,pos1,sign1,True,True)
            else:
                pos+=1
            genes_info[gene]=(species,scfid,pos,sign,False,False) # Last field: True if last gene of its scaffold
            (prev_sp,prev_scfid,prev_gene)=(sp,scfid,gene)
    (sp1,scfid1,pos1,sign1,last1,seen1)=genes_info[prev_gene]
    genes_info[prev_gene]=(sp1,scfid1,pos1,sign1,True,True)

    # Reading the kept adjacencies to check they all appear in the scaffolds, updating the genes_info for each seen adjacency
    print("  Checking that kept adjacencies are observed in the scaffolds")
    kept_adjacencies_stream=open(kept_adjacencies_file,"r").readlines()
    for adj in kept_adjacencies_stream:
        if adj[0]!="#":
            adj1=adj.rstrip().split()
            (sp,gene1,gene2,sign1,sign2)=(adj1[0],adj1[1],adj1[2],adj1[3],adj1[4])
            test1=(genes_info[gene1][0]==genes_info[gene2][0]   and genes_info[gene1][1]==genes_info[gene2][1]) # Same species and scaffold
            test2=(genes_info[gene1][2]==genes_info[gene2][2]-1 and genes_info[gene1][3]==sign1 and genes_info[gene2][3]==sign2) # Same adjacency
            test3=(genes_info[gene1][2]==genes_info[gene2][2]+1 and genes_info[gene1][3]==complement_sign(sign1) and genes_info[gene2][3]==complement_sign(sign2)) # Reversed adjacency
            if not (test1 and (test2 or test3)):
                print("ERROR: unobserved adjacency in the scaffolds\t"+adj)
            else:
                if genes_info[gene1][2]==genes_info[gene2][2]-1:
                    (sp1,scfid1,pos1,sign1,last1,seen1)=genes_info[gene1]
                    genes_info[gene1]=(sp1,scfid1,pos1,sign1,last1,True)
                else:
                    (sp1,scfid1,pos1,sign1,last1,seen1)=genes_info[gene2]
                    genes_info[gene2]=(sp1,scfid1,pos1,sign1,last1,True)

    # Checking that all scaffold adjacencies are in the kept adjacencies
    print("  Checking that scaffolds adjacencies are in the kept adjacencies")
    genes=sorted(genes_info.keys())
    for g in genes:
        g_info=genes_info[g]
        if g_info[4]==False and g_info[5]==False:
            print("ERROR: scaffold adjacency unobserved in the kept adjacencies\t"+g_info[0]+"."+g_info[1]+"."+g)

    # Reading the discarded adjacency to check that none could join two scaffolds
    print("  Checking that discarded adjacencies can not be added back")
    discarded_adjacencies_stream=open(discarded_adjacencies_file,"r").readlines()
    for adj in discarded_adjacencies_stream:
        if adj[0]!="#":
            adj1=adj.rstrip().split()
            (sp,gene1,gene2,sign1,sign2)=(adj1[0],adj1[1],adj1[2],adj1[3],adj1[4])
            test1=(genes_info[gene1][0]==genes_info[gene2][0]   and genes_info[gene1][1]==genes_info[gene2][1]) # Same species and scaffold 
            test2=(genes_info[gene1][4]==True and genes_info[gene2][3]==1 and genes_info[gene1][3]==sign1 and genes_info[gene2][3]==sign2)
            test3=(genes_info[gene1][4]==True and genes_info[gene2][4]==True and genes_info[gene1][3]==sign1 and genes_info[gene2][3]==complement_sign(sign2))
            test4=(genes_info[gene1][3]==1 and genes_info[gene2][4]==True and genes_info[gene1][3]==complement_sign(sign1) and genes_info[gene2][3]==complement_sign(sign2))
            test5=(genes_info[gene1][3]==1 and genes_info[gene2][3]==1 and genes_info[gene1][3]==complement_sign(sign1) and genes_info[gene2][3]==sign2)
            if test1 and (test2 or test3 or test4 or test5):
                print("ERROR: observable discarded adjacency linking scaffolds "+genes_info[gene1][0]+"."+genes_info[gene1][1]+"-"+genes_info[gene2][1]+"\t"+adj)

# --------------------------------------------------------------------------------------------
# Main function

# Reading parameters
assert(len(sys.argv)>=2), USAGE
P_ADJ_FILE_NAME=sys.argv[1]
assert(os.path.isfile(P_ADJ_FILE_NAME)), "Unknown adjacencies file"
if len(sys.argv) >= 3:
    P_SPECIES=sys.argv[2]
else:
    P_SPECIES="ALL"
if len(sys.argv) >= 4:
    P_THRESHOLD=-math.log1p(-float(sys.argv[3]))
else:
    P_THRESHOLD=-math.log1p(-EPSILON)
if len(sys.argv) >= 5:
    P_OUTPUT_PREFIX=sys.argv[4]
else:
    P_OUTPUT_PREFIX="linear_scaffolds"
if len(sys.argv) >= 6:
    P_ALGORITHM=sys.argv[5]  
else:
    P_ALGORITHM="M1"
assert (P_ALGORITHM in ['F1','F2','M1','M2']), "Unknown algorithm"
if len(sys.argv) >= 7:
    OBSERVED=float(sys.argv[6])
if len(sys.argv) >= 8:
    EPSILON=float(sys.argv[7])


# Reading all species from the adjacencies file name
ADJ_FILE_STREAM=open(P_ADJ_FILE_NAME,"r").readlines()
ALL_SPECIES_AUX={}
for adj in ADJ_FILE_STREAM:
    if adj[0]!="#": # "#" is used for comments
        l1 = adj.rstrip().split()
        species=l1[0]
        ALL_SPECIES_AUX[species]=1
ALL_SPECIES=sorted(list(ALL_SPECIES_AUX.keys()))

# Reading the species list
SPECIES_LIST=[]
if os.path.isfile(P_SPECIES): # If the species parameter is a file, we read the species from it
    SPECIES_STREAM=open(P_SPECIES,"r").readlines()
    for l in SPECIES_STREAM:
        if l[0]!="#":
            SPECIES_LIST.append(l.rstrip().split()[0])
elif P_SPECIES!="ALL": # One species is passed in parameters
    assert(P_SPECIES in ALL_SPECIES), "Unknown species"
    SPECIES_LIST.append(P_SPECIES)
else:
    SPECIES_LIST=list(ALL_SPECIES)
SPECIES_LIST.sort()

# Processing the species
output_obs_scfs=open(P_OUTPUT_PREFIX+"_"+P_ALGORITHM+"_obs_scaffolds","w")
output_obs_scfs.write("#species scaffold gene orientation\n")
output_scfs=open(P_OUTPUT_PREFIX+"_"+P_ALGORITHM+"_scaffolds","w")
output_scfs.write("#species scaffold gene orientation\n")
output_kept=open(P_OUTPUT_PREFIX+"_"+P_ALGORITHM+"_kept","w")
output_kept.write("#species gene1 gene2 orientation1 orientation2 prior_score posterior_score scaffold_id\n")
output_disc=open(P_OUTPUT_PREFIX+"_"+P_ALGORITHM+"_disc","w")
output_disc.write("#species gene1 gene2 orientation1 orientation2 prior_score posterior_score filtering_stage\n")
for species in SPECIES_LIST:
    process_genome(species,ADJ_FILE_STREAM,P_THRESHOLD,P_ALGORITHM,P_OUTPUT_PREFIX)
output_obs_scfs.close()
output_scfs.close()
output_kept.close()
output_disc.close()

check_results(P_OUTPUT_PREFIX+"_"+P_ALGORITHM+"_kept",P_OUTPUT_PREFIX+"_"+P_ALGORITHM+"_disc",P_OUTPUT_PREFIX+"_"+P_ALGORITHM+"_scaffolds")
