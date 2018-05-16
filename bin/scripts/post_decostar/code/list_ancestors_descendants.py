__author__="Cedric Chauve and Eric Tannier"
__email__="cedric.chauve@sfu.ca eric.tannier@inria.fr"
__date__="December 18, 2016"

import sys,script_tree,string,math

# ---------------------------------------------------------------------------
USAGE="python list_ancestors_descendants.py <gene_trees_file> <species_tree_file> <output_file> [optional:<number_missing_secies>]\n"
USAGE+="  gene_trees_file:   reconciled gene trees\n"
USAGE+="                     format = 2 lines per tree: line1=>family family_id; line2=tree in Newick format\n"
USAGE+="                     reconciliation events = Spe for speciation, Extant for extant genes, Dup for duplications (losses not considered)\n"
USAGE+="  species_tree_file: species tree file in Newick format\n"
USAGE+="  output_file:       output file\n"
USAGE+="                     format = 1 line per pair of ancestor-descendant gene: species1_id species2_id gene1_name gene2_name family_id\n"
USAGE+="                     species1 = ancestor, species2 = descendant, gene names as in the gene tree file\n"
USAGE+="                     output branches are only the pre-extant or pre-speciation branches\n"
USAGE+="  number_missing_secies: if given, then only gene trees with no duplication are considered\n"
USAGE+="                         and only families absent from at most number_missing_secies are kept\n"

# ---------------------------------------------------------------------------
# Reading parameters
file_name_gene_trees   = sys.argv[1]
file_name_species_tree = sys.argv[2]
output_file_name       = sys.argv[3]
# By default all families are listed; if a fourth parameter is given, only duplication-free families 
# that are absent from a bounded missing number of species are listed
keep_all_families=True
if len(sys.argv)==5:
    nb_missing_species = int(sys.argv[4])
    keep_all_families=False
    
# ---------------------------------------------------------------------------
# Reading the species tree and extracting all species, ancestra species and branches

species_tree = script_tree.readTree(open(file_name_species_tree,"r").readline())
all_species  = script_tree.getLeavesNames(species_tree)
nodes        = script_tree.getNodes(species_tree)

ancestral_species = []
extant_species    = []
species_branches  = [] # Pair of species, (parent,desc endant)

for n in nodes:
    if script_tree.isLeaf(species_tree,n):
        species = script_tree.getName(species_tree,n).split("@")[1]
        extant_species.append(species)
    else:
        species = script_tree.getBootstrap(species_tree,n)
        ancestral_species.append(species)
    if script_tree.isRoot(species_tree,n):
        root_species = species
    else:
        species_branches.append([script_tree.getBootstrap(species_tree,script_tree.getParent(species_tree,n)),species])

nb_extant_species=len(extant_species)  # Number of extant species
        
# ---------------------------------------------------------------------------
# Reading reconciled gene trees 

file_gene_trees = open(file_name_gene_trees,"r").readlines()
reconciled_gene_trees = {} # Array of gene trees indexed by family id

# Recording the Newick format of all considered gene trees
# Comments are prefixed by \#
for line in file_gene_trees:    
    if line[0]!="#" and line[0] == ">":
        id_gene_family = line.strip().split(" ")[-1]
    elif line[0]!="#" and line[0] != ">":
        # We keep only gene trees if no filter or if they have no duplications and are not missing too many species
        if keep_all_families==True or (line.count("Dup")==0 and line.count("Extant")>=nb_extant_species-nb_missing_species):
            reconciled_gene_trees[id_gene_family] = line.strip()

# Computing and writing the output, pairs of related genes
output_file = open(output_file_name,"w")
for id_gene_family in reconciled_gene_trees.keys():
    line  = reconciled_gene_trees[id_gene_family] # Line encoding the reconciled gene tree
    tree  = script_tree.readTree(line)            # Reconciled gene tree in script_tree format
    root  = script_tree.getRoot(tree)             # Root species
    nodes = script_tree.getNodes(tree)            # List of all tree nodes
    
    for n in nodes: # Loop on nodes=genes
        # Extracting the current gene name and species and reconciliation event
        if script_tree.isLeaf(tree,n):
            name    = script_tree.getName(tree,n)     
            species = name.split("|")[0].split("@")[0]
        else:
            name = script_tree.getBootstrap(tree,n)
        species = name.split("|")[2] # Species of the gene
        event   = name.split("|")[1] # Reconciliation event
        if event == "Extant":
            name = name.split("|")[0]
        elif event == "Spe":
            name = id_gene_family+"|"+name.split("|")[0]

        # Outputing the branch if it is a gene we want to consider, extant or speciation
        if event == "Extant" or event == "Spe":
            current = n     # Current node
            found   = False # Boolean indicating if we did find the corresponding branch to its ancestor
            while not found:
                if script_tree.isRoot(tree,current):
                    found = True
                else:
                    current         = script_tree.getParent(tree,current)
                    name_current    = script_tree.getBootstrap(tree,current)
                    event_current   = name_current.split("|")[1]
                    species_current = name_current.split("|")[2]
                    if event_current == "Spe" and [species_current,species] in species_branches:
                        output_file.write(species_current+" "+species+" "+str(id_gene_family)+"|"+name_current.split("|")[0]+" "+name+" "+str(id_gene_family)+"\n")
                        found = True
            

