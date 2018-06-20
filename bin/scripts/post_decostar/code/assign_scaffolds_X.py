import sys,script_tree,string,script_graphe,math


scaffold_file = sys.argv[1]
file_name_gene_trees = sys.argv[2]
file_species_tree = sys.argv[3]
genes_file = sys.argv[4]
sep = sys.argv[5]

MAXIMUM_NUMBER = 1000000000

# read species tree
# print "read species tree"
species_tree = script_tree.readTree(open(file_species_tree,"r").readline())
all_species = script_tree.getLeavesNames(species_tree)
ancestral_species = []
nodes = script_tree.getNodes(species_tree)
species_couples = []
for n in nodes:
    if script_tree.isLeaf(species_tree,n):
        species = script_tree.getName(species_tree,n).split(sep)[1]
    else:
        species = script_tree.getBootstrap(species_tree,n)
        ancestral_species.append(species)
    if script_tree.isRoot(species_tree,n):
        root_species = species


# read ancestral and extant scaffolds >= 2
# print "read ancestral and extant scaffolds"
file_scaffolds = open(scaffold_file,"r").readlines()
gene_coordinates = {}
chromosomes = {}
for line in file_scaffolds:
    if line[0] != "#":
        words = line.split()
        species = words[0]
        scaffold_name = words[1]
        name = words[2]
        direction = words[3]
        if name.find(sep) >= 0:
            chromosome = "EXT"+species+"_"+scaffold_name
            name = name.split(sep)[1]
        else:
            chromosome = "ANC"+species+"_"+scaffold_name
        gene_coordinates[name] = [species,chromosome,direction]
        if not chromosomes.has_key(chromosome):
           chromosomes[chromosome] = []
        chromosomes[chromosome].append(name)

# FILL GENE COORDINATES (USEFUL TO HAVE ACCESS TO THE CHROMOSOME FOR EACH GENE)
file_coordinates = open(genes_file,"r").readlines()
# print "read reference scaffolds"
reference_scaffold = {}
sumX = 0
sum2L = 0
sum2R = 0
sum3L = 0
sum3R = 0
total = sumX+sum2L+sum2R+sum3L+sum3R
for line in file_coordinates:
    if line[0] != '#':    
        words = line.split()
        species = words[0]
        chromosome = words[1]
        family = words[2]
        name = words[3]
        direction = words[4]
        start = int(words[5])
        stop = words[6]
        reference_scaffold[name] = [species,chromosome]
        if chromosome == "X":
            sumX = sumX + 1
        if chromosome == "2L":
              sum2L = sum2L + 1
        if chromosome == "2R":
            sum2R = sum2R + 1
        if chromosome == "3L":
            sum3L = sum3L + 1
        if chromosome == "3R":
            sum3R = sum3R + 1
total = sumX+sum2L+sum2R+sum3L+sum3R

# read reconciled trees 
# print "read reconciled trees"
reconciled_gene_trees = {}
file_gene_trees = open(file_name_gene_trees,"r").readlines()
for line in file_gene_trees:
    if line[0] == ">":
        id_gene_family = line.strip().split(" ")[-1]
    else:
        reconciled_gene_trees[id_gene_family] = line.strip()
        tree = script_tree.readTree(line)
        root = script_tree.getRoot(tree)
        leaves = script_tree.getLeaves(tree,root)
        nodes = script_tree.getNodes(tree)

# FIND GAMBIA ORTHOLOG
#print species_couples
# print "find orthologs in reconciled gene trees"
orthologs = {}
for id_gene_family in reconciled_gene_trees.keys():
    #print "tree",id_gene_family  
    line = reconciled_gene_trees[id_gene_family]
    tree = script_tree.readTree(line)
    root = script_tree.getRoot(tree)
    leaves = script_tree.getLeaves(tree,root)
    nodes = script_tree.getNodes(tree)
    gambia = []
    # first list gambia genes in this tree
    for l in leaves:
    event = script_tree.getName(tree,l).split("|")[1]
    if event != "Loss":
        #print script_tree.getName(tree,l)
        name_leaf = script_tree.getName(tree,l).split("|")[0].split(sep)[1]
        #print name_leaf
        if gene_coordinates.has_key(name_leaf):
            gene_coordinates[name_leaf].append(id_gene_family)
            #print reference_scaffold[name_leaf]
            if reference_scaffold[name_leaf][0] == "Anopheles_gambiae" and (reference_scaffold[name_leaf][1] == "X" or len(reference_scaffold[name_leaf][1]) == 2):
                gambia.append(l)
    # then for all nodes Extant, Spe, find orthologs
    for n in nodes:
        if script_tree.isLeaf(tree,n):
            name = script_tree.getName(tree,n)
            species_name = name.split("|")[0].split(sep)[0]
        else:
            name = script_tree.getBootstrap(tree,n)
        species = name.split("|")[2]
        event = name.split("|")[1]
        if event == "Extant" or event == "Spe":
            if event == "Extant":
                name = name.split("|")[0].split(sep)[1]
            if event == "Spe":
                name = id_gene_family+"|"+name.split("|")[0]
            orthologs[name] = []
            if not gene_coordinates.has_key(name):
                if event == "Extant":
                    contig_name = "EXT"+species+"_GENE"+name
                else:
                    contig_name = "ANC"+species+"_GENE"+name
                gene_coordinates[name] = [species,contig_name,"+"]
                chromosomes[contig_name] = [name]
            for l in gambia:
                common_ancestor = script_tree.lastCommonAncestor(tree,n,l)
                #print common_ancestor,script_tree.getBootstrap(tree,common_ancestor)
                if script_tree.isLeaf(tree,common_ancestor) or (script_tree.getBootstrap(tree,common_ancestor).split("|")[1] == "Spe"):
                    orthologs[name].append(script_tree.getName(tree,l).split("|")[0].split(sep)[1])
         

   
#for k in orthologs.keys():
    #print k,orthologs[k]
  
# distribution for each chromosome function of orthologs
# print "assign scaffolds to chromosomes"
chr_names = chromosomes.keys()
distrib = {}
for c in chr_names:
    #print c,len(chromosomes[c])
    distrib[c] = {}
    nb_ortholog = 0
    #if len(chromosomes[c]) >= 1:
    for g in chromosomes[c]:
        #print g,orthologs[g]
        for o in orthologs[g]:
            nb_ortholog = nb_ortholog + 1
            chromosome = reference_scaffold[o][1]
            #print c,g,o,chromosome
            if not distrib[c].has_key(chromosome):
                distrib[c][chromosome] = 0
            distrib[c][chromosome] = distrib[c][chromosome] + 1.0/len(orthologs[g])
    somme = sum(distrib[c].values())
    for i in distrib[c].keys():
        distrib[c][i] = distrib[c][i]/somme
    if distrib[c] == {}:
        distrib[c] = {"X":sumX/float(total),"2L":sum2L/float(total),"2R":sum2R/float(total),"3L":sum3L/float(total),"3R":sum3R/float(total)}
    for k in ["X","2L","2R","3L","3R"]:
      if not distrib[c].has_key(k):
        distrib[c][k] = 0.0


# assign genes to chromosomes according to their belonging to the scaffold
# print "assign the genes to chromosomes" #,len(gene_coordinates.keys())
for g in gene_coordinates.keys():
    c = gene_coordinates[g][1]
    distrib[g] = distrib[c]

nb_genes = {"X":0.0,"2L":0.0,"2R":0.0,"3L":0.0,"3R":0.0}
for k in distrib.keys():
    #print k,distrib[k]
    if gene_coordinates.has_key(k):
        for chrom in distrib[k].keys():
            nb_genes[chrom] = nb_genes[chrom] + distrib[k][chrom]
print nb_genes,sum(nb_genes.values()),str(100*nb_genes["X"]/sum(nb_genes.values()))+"% in X"

assignment_file = open(scaffold_file+"_assignment","w")
assignment_file.write("#scaffold_name scaffold_length prob_X prob_2L prob_2R prob_3L prob_3R\n")
for c in chr_names:
    assignment_file.write(c+" "+str(len(chromosomes[c]))+" "+str(distrib[c]["X"])+" "+str(distrib[c]["2L"])+" "+str(distrib[c]["2R"])+" "+str(distrib[c]["3L"])+" "+str(distrib[c]["3R"])+"\n")

