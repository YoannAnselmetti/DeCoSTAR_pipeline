__author__="Cedric Chauve"
__email__="cedric.chauve@sfu.ca"
__date__="November 3, 2016"

import sys
import os

# This script computes an SCJ distance between a pair of potentially fragmented genomes
# Input:
# - an orthologous families file, in the format
#   ancestor_species_id descendant_species_id ancestral_gene(gene_name or NIL) descendant_gene(gene_name or NIL)
# - the path to the scaffold file
# - a list of branches (pairs of species ancestor descendant), with "ALL" for all branches seen in the families file
# - a filter option 0/1/2/3
#   0=no filter, 1=filter genes with no orthologs, 2=filter genes with 0 or >1 orthologs
#   3=filter 2 + filter genes that have no conserved adjacency (single-gene synteny block)
#   where a synteny block is a set of genes with similar order and orientations
#   this option also returns the list of the length of the synteny blocks (standard output)
#
# Scaffolding files are in the format
#   species scaffold gene orientation<TAB><extra:optional>
#   with possibly other lines such as # ... (comment) or > ... (scaffold header)
#
# Output:
# - a file of non-conserved adjacencies in the format
#   species_id gene_1 gene_2 loss_reason
#                            NO for no ortholog
#                            DSCJ if at least one gene has no ortholog
#                            SCJ if orthologs exist for both genes but do not form an adjacency or potential adjacency in the other species
#                            PSCJ if orthologs exist for both genes but do not form an adjacency but do form a potential adjacency in the other species
# A potential adjacency exist if two orthologs are at extremities of scaffolds (in the proper orientation)
#
# Assumption: all genes present in the orthogroup file are also present in the scaffolds
# If they are not they are added by the script

# ----------------------------------------------------------------------------------------
# Reading ortholog groups

# Takes as input the path to the orthogroups file
# Returns an array, indexed by the pairs (species,gene), of a list of the orthologs
def read_orthogroups(F_OG_stream,anc,desc,F_SCF_stream):
    OG={}
    # Initialisation
    for l in F_OG_stream:
        l1=l.rstrip().split()
        (ancestor,descendant,anc_gene,desc_gene)=(l1[0],l1[1],l1[2],l1[3])
        if anc==ancestor and desc==descendant:
            OG[(ancestor,anc_gene)]=[]
            OG[(descendant,desc_gene)]=[]
    for l in F_SCF_stream:
        if l[0]!="#":
            l1=l.rstrip().split()
            (species,gene)=(l1[0],l1[2])
            if species==anc or species==desc:
                OG[(species,gene)]=[]
    # Filling the tables
    for l in F_OG_stream:
        l1=l.rstrip().split()
        (ancestor,descendant,anc_gene,desc_gene)=(l1[0],l1[1],l1[2],l1[3])        
        if anc==ancestor and desc==descendant:
            if desc_gene!="NIL":
                OG[(ancestor,anc_gene)].append((descendant,desc_gene))
            if anc_gene!="NIL":
                OG[(descendant,desc_gene)].append((ancestor,anc_gene))
    return(OG)

# ----------------------------------------------------------------------------------------
# Reading a scaffold file to record all observed adjacencies

# Takes as input a scaffold file path
# Return a list of quadruples(species,scaffold,gene,sign)
# Add single-gene scaffolds for genes in orthogroups not seen in scaffolds
# Two options: 1. filter genes with no orthologs in the other species (Filter=True)
#              2. keep such genes (filter=False)

SIGN2EXT1={"+":"tail", "-":"head"}
SIGN2EXT2={"+":"head", "-":"tail"}
SWITCHSIGN={"+":"-","-":"+"}

EXTREMITY=("XX","NIL","+") # dummy neighbour of a fragment extremity

def check_filter(sp,gene,f,OG):
    filter1=(f==1 and len(OG[(sp,gene)])==0)
    filter2=(f in [2,3] and (len(OG[(sp,gene)])!=1 or (len(OG[(sp,gene)])==1 and len(OG[OG[(sp,gene)][0]])!=1)))
    return((filter1,filter2))
        
def read_scaffold(F_stream,OG,f,species):
    SCFS={} # Array encoding the scaffolds
    # Recording genes incorporated in the scaffolds array
    present_genes={}
    OG_LIST=sorted(list(OG.keys()))
    for (sp,gene) in OG_LIST: # Initialization of the present genes array
        (filter1,filter2)=check_filter(sp,gene,f,OG)
        if sp==species and not (filter1 or filter2):
            present_genes[(sp,gene)]=False
    # Encoding the scaffolds
    index=0
    for l in F_stream:
        if l[0]!="#" and l[0]!=">": # gene line
            l1=l.rstrip().split()
            if l1[0]==species:
                (sp,scf,gene,sign)=(l1[0],l1[1],l1[2],l1[3])
                (filter1,filter2)=check_filter(sp,gene,f,OG)
                if not (filter1 or filter2):
                    SCFS[index]=(sp,scf,gene,sign)
                    index+=1
                    present_genes[(sp,gene)]=True
    # Adding single-gene scaffolds
    scf_id=1 # For ID of added scaffolds
    for (sp,gene) in OG_LIST: # Loop on all genes from orthogroups
        (filter1,filter2)=check_filter(sp,gene,f,OG)
        if sp==species and (not (filter1 or filter2)) and present_genes[(sp,gene)]==False:
            SCFS[index]=(sp,"SGSCF"+str(scf_id),gene,"+")
            index+=1
            scf_id+=1
    return(SCFS)

# Apply filter 3 to scaffolds, to remove single-gene synteny blocks and generate the distribution of synteny block length
# Assume filter 2 has already been applied

## Auxiliary function that associates a unique id to each gene such that genes in 1-to-1 orthogroup have the same id
def og_assign_id(og):
    og_id=0
    og_id_array={}
    og_id_inv={}
    og_list=sorted(list(og.keys()))
    for (sp1,g1) in og_list:
        if len(og[(sp1,g1)])==1:
            [(sp2,g2)]=og[(sp1,g1)]
            if sp1<sp2:
                og_id_array[g1]=og_id
                og_id_array[g2]=og_id
                og_id_inv[og_id]={}
                og_id_inv[og_id][sp1]=g1
                og_id_inv[og_id][sp2]=g2
                og_id+=1
    return((og_id_array,og_id_inv))
                
## Auxiliary function that returns a list of adjacencies and their positions in a scaffold array scf of length lg
def scaffolds_2_list(scf,lg,og_id):
    adj_list=[]   # Records the adjacencies
    adj_pos={}    # Position of an adjacency
    adj_status={} # Status of an adjacency at a given position(set at 0 here)
    # Scanning genome
    for i in range(1,lg):
        (sp,scf1,g1,sign1)=scf[i-1]
        ext1=SIGN2EXT2[sign1]
        (sp,scf2,g2,sign2)=scf[i]        
        ext2=SIGN2EXT1[sign2]
        if (scf1==scf2):
            if g1<g2:
                adj_list.append((og_id[g1],ext1,og_id[g2],ext2))
            else:
                adj_list.append((og_id[g2],ext2,og_id[g1],ext1))
        adj_pos[(g1,ext1,g2,ext2)]=i
        adj_pos[(g2,ext2,g1,ext1)]=i
        adj_status[i]=False
    return((adj_list,adj_pos,adj_status))

## Actually filtering a scaffold to keep only genes with at least one conserved adjacency
def filter3(in_scf,lg,status):
    out_scf={}
    newindex=0
    for i in range(0,lg):
        (sp,scf,g,sign)=in_scf[i] 
        if (i<lg-1 and status[i+1]==True) or (i>0 and status[i]==True):
            out_scf[newindex]=(sp,scf,g,sign)
            newindex+=1
    return(out_scf)
            
def filter_scaffolds_3(og,scf_anc,scf_desc,anc,desc):
    lg_anc=len(scf_anc)   # number of genes in the ancestor scaffolds
    lg_desc=len(scf_desc) # number of genes in the descendant scaffolds
    assert(lg_anc==lg_desc), "ERROR: after filter 2 ancestor and descendants do not have the same gene content"
    # Assigning id to ogs
    (og_id,og_id_inv)=og_assign_id(og)
    
    # Scanning ancestor genome
    (adj_anc_list,adj_anc_pos,adj_anc_status)=scaffolds_2_list(scf_anc,lg_anc,og_id)
    # Scanning descendant genome
    (adj_desc_list,adj_desc_pos,adj_desc_status)=scaffolds_2_list(scf_desc,lg_desc,og_id)

    # Intersection of both lists:
    adj_common_list = list(set(adj_anc_list) & set(adj_desc_list))
    
    # Marking conserved adjacencies in the status arrays
    for (id1,ext1,id2,ext2) in adj_common_list:
        (g1,g2)=(og_id_inv[id1][anc],og_id_inv[id2][anc])
        anc_pos=adj_anc_pos[(g1,ext1,g2,ext2)]
        adj_anc_status[anc_pos]=True
        (g1,g2)=(og_id_inv[id1][desc],og_id_inv[id2][desc])        
        desc_pos=adj_desc_pos[(g1,ext1,g2,ext2)]
        adj_desc_status[desc_pos]=True

    # Creating the new scaffolds
    new_scf_anc=filter3(scf_anc,lg_anc,adj_anc_status)
    new_scf_desc=filter3(scf_desc,lg_desc,adj_desc_status)
    str1="\tpre-filtering for single-gene blocks #genes "+str(len(scf_anc))+" #adj anc "+str(len(adj_anc_list))+" #adj desc "+str(len(adj_desc_list))
    str2="\tpost-filtering for single-gene blocks #genes "+str(len(new_scf_anc))+" #adj "+str(len(adj_common_list))
    print(str1+str2)

    # Generating the distribution of the length of the synteny blocks
    lg_synt_blocks_dist_aux={}
    for lg in range(0,lg_anc):
        lg_synt_blocks_dist_aux[lg]=0
    current_block_lg=1
    for i in range(1,lg_anc):
        if adj_anc_status[i]==True:
            current_block_lg+=1
        else:
           lg_synt_blocks_dist_aux[current_block_lg]+=1
           current_block_lg=1
    lg_synt_blocks_dist=[]
    lg_synt_blocks_dist_aux_list=sorted(list(lg_synt_blocks_dist_aux.keys()))
    for lg in lg_synt_blocks_dist_aux_list:
        if (lg>1 and lg_synt_blocks_dist_aux[lg]>0):
            lg_synt_blocks_dist.append((lg,lg_synt_blocks_dist_aux[lg]))
    lg_synt_blocks_dist.sort()

    return((new_scf_anc,new_scf_desc,lg_synt_blocks_dist))
           
# SCF = array encoding the (possibly filtered for genes with orthologs) scaffolds
def adjacencies_from_scaffold(SCF):
    ADJS={} # Array of the adjacencies, indexed by gene extremities
    current_scf="0" # Assumption: scaffold ID can not be the string "0"
    prev_gene="NIL"
    prev_sign="+"
    SCF_LIST=sorted(list(SCF.keys()))
    for i in SCF_LIST:
        (species,scf,gene,sign)=SCF[i]
        if scf!=current_scf: # starting new scaffold
            ADJS[(species,prev_gene,SIGN2EXT2[prev_sign])]=EXTREMITY # Last gene of previous scaffold
            ADJS[(species,gene,SIGN2EXT1[sign])]=EXTREMITY           # First gene of current scaffold
        else: # continuing the current scaffold
            ADJS[(species,prev_gene,SIGN2EXT2[prev_sign])]=(species,gene,SIGN2EXT1[sign]) # Previous gene
            ADJS[(species,gene,SIGN2EXT1[sign])]=(species,prev_gene,SIGN2EXT2[prev_sign]) # Current gene
        current_scf=scf
        prev_gene=gene
        prev_sign=sign
    if len(SCF_LIST)>0:
        ADJS[(species,prev_gene,SIGN2EXT2[prev_sign])]=EXTREMITY # Last gene of last scaffold
    return(ADJS)            

# ----------------------------------------------------------------------------------------
# Processing a scaffold file to record all broken adjacencies

# Takes as input a scaffold, an orthogroups array, and adjacencies array for the other species
# Returns a list of lost adjacencies, each with a status in {NO,SCJ,PSCJ}

def process_scaffold(SCF,OG,ADJS):
    SCJ_ADJS=[]     # List of non-conserved adjacencies
    current_scf="0" # Assumption: scaffold ID can not be the string "0"
    prev_gene="NIL" 
    prev_sign="+"
    SCF_LIST=sorted(list(SCF.keys()))
    for i in SCF_LIST: # Loop on the array encoding all scaffolds    
        (species,scf,gene,sign)=SCF[i]
        if scf==current_scf: # continuing the current scaffold, adjacency to check
            og_gene=OG[(species,gene)]
            og_prev_gene=OG[(species,prev_gene)]
            scj_status="SCJ" # By default, the adjacency is lost
            if len(og_gene)==0 or len(og_prev_gene)==0:
                scj_status="DSCJ" # The adjacency is lost due to gene content difference
            else:
                for (s1,g1) in og_gene:
                    G1=(s1,g1,SIGN2EXT1[sign])
                    for (s2,g2) in og_prev_gene:
                        G2=(s2,g2,SIGN2EXT2[prev_sign])
                        #print(str(G1)+" "+str(G2)+" "+str(len(og_gene))+" "+str(len(og_prev_gene)))
                        if ADJS[G1]==G2:
                            scj_status="NO" # The adjacency is conserved
                        elif scj_status=="SCJ" and ADJS[G1]==EXTREMITY and ADJS[G2]==EXTREMITY:
                            scj_status="PSCJ"
            if scj_status!="NO":
                SCJ_ADJS.append((species,prev_gene,prev_sign,gene,sign,scj_status))
        current_scf=scf
        prev_gene=gene
        prev_sign=sign
    return(SCJ_ADJS)

# -------------------------------------------------------------------------------------
# Main body

P_OG_FILE_NAME=sys.argv[1]          # Orthogroups file
assert(os.path.isfile(P_OG_FILE_NAME)), "ERROR: Unknown orthologous relations file"
P_BRANCHES_FILE_NAME=sys.argv[2]    # List of branches where to compute SCJ operations
P_SCF_FILE_NAME=sys.argv[3]         # Name of scaffold file
assert(os.path.isfile(P_SCF_FILE_NAME)), "ERROR: Unknown scaffolds file"
P_FILTER=int(sys.argv[4])           # 0 (no filter), 1 (filter genes with no ortholog), 2 (filter also genes with more than 1 ortholog)
assert(P_FILTER in [0,1,2,3]), "Error: the filter parameter should be 0, 1, 2 or 3"

P_OG_FILE_STREAM=open(P_OG_FILE_NAME,"r").readlines()
P_SCF_FILE_STREAM=open(P_SCF_FILE_NAME,"r").readlines()
if P_BRANCHES_FILE_NAME=="ALL": # Reading branches in the orthogroups file
    BRANCHES_AUX={}
    for l in P_OG_FILE_STREAM:
        if l[0]!="#":
            l1=l.rstrip().split()
            BRANCHES_AUX[(l1[0],l1[1])]=1
    BRANCHES=sorted(list(BRANCHES_AUX.keys()))
elif os.path.isfile(P_BRANCHES_FILE_NAME): # Reading the branches in a file
    P_BRANCHES_FILE_STREAM=open(P_BRANCHES_FILE_NAME,"r").readlines()
    BRANCHES={}
    for l in P_BRANCHES_FILE_STREAM:
        if l[0]!="#":
            l1=l.rstrip().split()
            BRANCHES.append((l1[0],l1[1]))
else:
     print("ERROR: Unknown branches file")       

OUTPUT_FILE_NAME=sys.argv[5]
OUTPUT_FILE_STREAM=open(OUTPUT_FILE_NAME,"w")
OUTPUT_FILE_STREAM.write("#ancestor descendant species\tgene1 gene2 sign1 sign2\tscj_status\n")

for (ANCESTOR,DESCENDANT) in BRANCHES:
    print("#Branch "+ANCESTOR+"-->"+DESCENDANT+"-----------------")
    OG=read_orthogroups(P_OG_FILE_STREAM,ANCESTOR,DESCENDANT,P_SCF_FILE_STREAM)
    SCFS_ANC_AUX=read_scaffold(P_SCF_FILE_STREAM,OG,P_FILTER,ANCESTOR)
    SCFS_DESC_AUX=read_scaffold(P_SCF_FILE_STREAM,OG,P_FILTER,DESCENDANT)
    if (P_FILTER==3):
        (SCFS_ANC,SCFS_DESC,SB_LG_DIST)=filter_scaffolds_3(OG,SCFS_ANC_AUX,SCFS_DESC_AUX,ANCESTOR,DESCENDANT)
        str_sb_lg_dist='\tSynt.Blocks.distribution\t'
        for (lg,nb) in SB_LG_DIST:
            str_sb_lg_dist+="["+str(lg)+","+str(nb)+"] "
        print(str_sb_lg_dist)
    else:
        (SCFS_ANC,SCFS_DESC)=(SCFS_ANC_AUX,SCFS_DESC_AUX)
    ADJS_ANC=adjacencies_from_scaffold(SCFS_ANC)
    ADJS_DESC=adjacencies_from_scaffold(SCFS_DESC)
    LOST_ADJS_ANC=process_scaffold(SCFS_ANC,OG,ADJS_DESC)
    LOST_ADJS_DESC=process_scaffold(SCFS_DESC,OG,ADJS_ANC)
    
    OUTPUT_FILE_STREAM.write("#"+ANCESTOR+" "+DESCENDANT+" "+str(len(SCFS_ANC.keys()))+" "+str(len(SCFS_DESC.keys()))+"\n")
    for (species,gene1,sign1,gene2,sign2,scj_status) in LOST_ADJS_ANC:
        OUTPUT_FILE_STREAM.write(ANCESTOR+" "+DESCENDANT+" "+species+"\t"+gene1+" "+gene2+" "+sign1+" "+sign2+"\t"+scj_status+"\n")
    for (species,gene1,sign1,gene2,sign2,scj_status) in LOST_ADJS_DESC:
        OUTPUT_FILE_STREAM.write(ANCESTOR+" "+DESCENDANT+" "+species+"\t"+gene1+" "+gene2+" "+sign1+" "+sign2+"\t"+scj_status+"\n")
       
