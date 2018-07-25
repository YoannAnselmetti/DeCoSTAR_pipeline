import sys,script_tree
## The goal of this script is to provide a branchwise rearrangement rate in X, plus duplications, losses and gene moves in X and autosomes
#
# Input: 
# - reconciled gene trees
# - orthogroups
# - scaffolds
# - assignments
# - all_scaffolds_M1_assignment contains, for each scaffold (name in 1st column) its probability to be affected to X
# - all_scaffolds_M1 the correspondance between genes and scaffolds
# - SCJ file, to know to which chromosome to assign rearragnements.
#
#
# 0/ Each gene, ancestral or extant, is assigned to a scaffold (from file scaffolds)
#    and each scaffold, ancestral or extant, has a probability to be on the X (from file assignment)
#    as a consequence each gene has also a probability to be on the X (except if is it alone on its scaffold)
#
# 1/ computes gene moves
#     from orthogroups, which is a liste of ancestor-descendants, outputs a gene move if
#           the ancestor has prob to be on the X < 0.2, and the descendant >0.8, or the opposite     
#
# 2/ computing duplications and losses
#     from reconciled trees, examine all duplications and losses.
#       2.1 Let N be the node with dup or loss, and P the closest ancestor with Spe (if there are none, NA)
#       2.2 If P does not have a probability to be on the X, NA
#       2.3 Else, add its probability to dupX or LossX, and 1-its probability to dupA or LossA
#
# 3/ counting SCJ 
#     from file scj, for each line, examine the two genes involve.
#     If they both have probability >0.6 assign the SCJ to X, if both <0.4, assign the SCJ to Aut
#       The mean SCJ rate is the number of SCJ/number of genes used to compute SCJ on each branch.
#
# The tree at the end has branch lengths number of SCJ and label % of SCJ on X
#
#

file_trees = open(sys.argv[1],"r").readlines()
file_orthogroups = open(sys.argv[2],"r").readlines()
file_scaffolds = open(sys.argv[3],"r").readlines()
file_assignments = open(sys.argv[4],"r").readlines()
file_scj = open(sys.argv[5],"r").readlines()
species_tree = script_tree.readTree(open(sys.argv[6],"r").readline())

statistics_per_node = {}
Sroot=script_tree.getRoot(species_tree)
Snodes = script_tree.getNodes(species_tree)
extant_species = []
for n in Snodes:
    if script_tree.isLeaf(species_tree,n):
        name = script_tree.getName(species_tree,n).split("@")[1]
        extant_species.append(name)
    else:
        name = script_tree.getBootstrap(species_tree,n)
    statistics_per_node[name] = {"#genes":0,"#genesX":0,"#genesA":0,"#genes_alone":0,"#scaffolds":0,"#scaffoldsX":0,"#scaffoldsA":0,"#SCJ":0,"#SCJX":0,"#SCJA":0}
    
#present_scaffold={} # Number of genes on scaffolds
#for line in file_assignments:
    #if line[0] != "#":
        #words = line.split()
        #species = words[0][3:].split("_")[0]
        #scaffold = words[0][3:].split("_")[1]
        #present_scaffold[(species,scaffold)]=0
        
#for line in file_scaffolds:
    #if line[0] != "#":
        #words = line.split()
        #species = words[0]
        #scaffold = words[1]
        #present_scaffold[(species,scaffold)]+=1
        
distrib = {}
for line in file_assignments:
    if line[0] != "#":
        words = line.split()
        distrib[words[0][3:]] = float(words[2])
        species = words[0][3:].split("_")[0]
        scaffold = words[0][3:].split("_")[1]
        #nb_genes = present_scaffold[(species,scaffold)]
        nb_genes = int(words[1])
        statistics_per_node[species]["#genes"] = statistics_per_node[species]["#genes"] + nb_genes
        statistics_per_node[species]["#genesX"] = statistics_per_node[species]["#genesX"] + nb_genes*float(words[2])
        statistics_per_node[species]["#genesA"] = statistics_per_node[species]["#genesA"] + nb_genes*(1.0-float(words[2]))
        statistics_per_node[species]["#scaffolds"] = statistics_per_node[species]["#scaffolds"] + 1
        statistics_per_node[species]["#scaffoldsX"] = statistics_per_node[species]["#scaffoldsX"] + float(words[2])
        statistics_per_node[species]["#scaffoldsA"] = statistics_per_node[species]["#scaffoldsA"] + 1.0-float(words[2])
        if words[1] == "1":
            statistics_per_node[species]["#genes_alone"] = statistics_per_node[species]["#genes_alone"] + 1

gene_coordinates = {}
for line in file_scaffolds:
    if line[0] != "#":
        words = line.split()
        #print words
        gene_coordinates[words[2]] = [words[0],words[1],distrib[words[0]+"_"+words[1]]]

        
# trees for duplications, losses, movements
print "counting duplications and losses"
events = {"Loss":{"A":0.0,"X":0.0,"NA":0.0},
          "Dup":{"A":0.0,"X":0.0,"NA":0.0}}
dup_per_branch = {}
in_unicopy_family = {}
for line in file_trees:
    if line[0] == ">":
        id_family = line[1:].strip().split()[1]
    else:
        tree = script_tree.readTree(line)
        leaves = script_tree.getLeavesNames(tree)
        unicopy = (line.find("Dup") < 0) and len(leaves) == 18
        nodes = script_tree.getNodes(tree)
        for n in nodes:
                if script_tree.isLeaf(tree,n):
                    name = script_tree.getName(tree,n)
                    in_unicopy_family[name.split("|")[0]] = unicopy
                else:
                    name = script_tree.getBootstrap(tree,n)
                    in_unicopy_family[id_family+"|"+name.split("|")[0]] = unicopy
                event = name.split("|")[1]
                if event == "Loss" or event == "Dup":
                    species_event = name.split("|")[2]
                    if not dup_per_branch.has_key(species_event):
                        dup_per_branch[species_event] = 0
                    if event == "Dup":
                        dup_per_branch[species_event] = dup_per_branch[species_event] + 1
                    if not script_tree.isRoot(tree,n):
                        current = script_tree.getParent(tree,n)
                        while (not script_tree.isRoot(tree,current)) and (script_tree.getBootstrap(tree,current).split("|")[1] != "Spe"):
                            current = script_tree.getParent(tree,current)
                        if script_tree.getBootstrap(tree,current).split("|")[1] == "Spe":
                            gene = id_family+"|"+script_tree.getBootstrap(tree,current).split("|")[0]
                            if gene_coordinates.has_key(gene):
                                events[event]["X"] = events[event]["X"] + distrib[gene_coordinates[gene][0]+"_"+gene_coordinates[gene][1]]
                                events[event]["A"] = events[event]["A"] + 1 - distrib[gene_coordinates[gene][0]+"_"+gene_coordinates[gene][1]]
                            else:
                                events[event]["NA"] = events[event]["NA"] + 1
                        else:
                            events[event]["NA"] = events[event]["NA"] + 1
                    else:
                        events[event]["NA"] = events[event]["NA"] + 1
#print events

print "couting single gene moves"
moves = {"AtoX":0.0,"XtoA":0.0,"AtoXu":0.0,"XtoAu":0.0}
for line in file_orthogroups:
    words = line.split()
    if gene_coordinates.has_key(words[2]) and gene_coordinates.has_key(words[3]):
        distrib1 = distrib[gene_coordinates[words[2]][0]+"_"+gene_coordinates[words[2]][1]]
        distrib2 = distrib[gene_coordinates[words[3]][0]+"_"+gene_coordinates[words[3]][1]]
        if distrib1 < 0.2 and distrib2 > 0.8:
            moves["AtoX"] = moves["AtoX"] + 1
            if in_unicopy_family[words[2]]:
                moves["AtoXu"] = moves["AtoXu"] + 1
        if distrib2 < 0.2 and distrib1 > 0.8:
            moves["XtoA"] = moves["XtoA"] + 1
            if in_unicopy_family[words[2]]:
                moves["XtoAu"] = moves["XtoAu"] + 1
#print moves


print "counting rearrangements"
branches = {}
Pbranches = {}
for line in file_scj:
    if line[0] != "#":
        words = line.strip().split()
        if words[7] == "SCJ":
            if gene_coordinates[words[3]][2] > 0.6 and gene_coordinates[words[4]][2] > 0.6:
                branches[words[0]+"_"+words[1]]["X"] = branches[words[0]+"_"+words[1]]["X"] + 1
            elif gene_coordinates[words[3]][2] < 0.4 and gene_coordinates[words[4]][2] < 0.4:
                branches[words[0]+"_"+words[1]]["A"] = branches[words[0]+"_"+words[1]]["A"] + 1
            else:
                branches[words[0]+"_"+words[1]]["NA"] = branches[words[0]+"_"+words[1]]["NA"] + 1
        if words[7] == "PSCJ":
            if gene_coordinates[words[3]][2] > 0.6 and gene_coordinates[words[4]][2] > 0.6:
                Pbranches[words[0]+"_"+words[1]]["X"] = Pbranches[words[0]+"_"+words[1]]["X"] + 1
            elif gene_coordinates[words[3]][2] < 0.4 and gene_coordinates[words[4]][2] < 0.4:
                Pbranches[words[0]+"_"+words[1]]["A"] = Pbranches[words[0]+"_"+words[1]]["A"] + 1
            else:
                Pbranches[words[0]+"_"+words[1]]["NA"] = Pbranches[words[0]+"_"+words[1]]["NA"] + 1
    elif len(line.split()) == 4:
        words = line[1:].strip().split()
        branches[words[0]+"_"+words[1]] = {"X":0,"A":0,"NA":0,"number":int(words[2])}
        Pbranches[words[0]+"_"+words[1]] = {"X":0,"A":0,"NA":0,"number":int(words[2])}                                                  

print "writing results"
total_X = 0.0
total = 0.0
Ptotal_X = 0.0
Ptotal = 0.0

for b in branches.keys():
    #print b
    branch_length = (branches[b]["X"] + branches[b]["A"] + branches[b]["NA"])
    Plength = (Pbranches[b]["X"] + Pbranches[b]["A"] + Pbranches[b]["NA"])
    ancestor = b.split("_")[0]
    descendant = b.split("_")[1]
    for n in Snodes:
        #print n
        if script_tree.isLeaf(species_tree,n):
            #print script_tree.getName(species_tree,n)
            if len(script_tree.getName(species_tree,n).split("@")) > 1:
                name = script_tree.getName(species_tree,n).split("@")[1]
            else:
                name = -1
        else:
            name = script_tree.getBootstrap(species_tree,n)
        if name == descendant:
            #print n
            #branch_length = length
            if branch_length == 0:
                branch_length = 1e-6
            script_tree.setLength(species_tree,n,branch_length)
            prop = int(100*float(branches[b]["X"])/branch_length)
            if script_tree.isLeaf(species_tree,n):
                script_tree.setName(species_tree,n,script_tree.getName(species_tree,n).split("@")[0]+"[&label="+str(prop)+"]")
            else:
                script_tree.writeBootstrap(species_tree,n,"[&label="+str(prop)+"]")
            statistics_per_node[name]["#SCJ"] = branch_length
            statistics_per_node[name]["#SCJX"] = branches[b]["X"]
            statistics_per_node[name]["#SCJA"] = branches[b]["A"]
            statistics_per_node[name]["#SCJNA"] = branches[b]["NA"]
    if branch_length > 0:
        total_X = total_X + float(branches[b]["X"])
        total = total + branch_length
        Ptotal_X = Ptotal_X + float(Pbranches[b]["X"])
        Ptotal = Ptotal + Plength
        
#for n in nodes:
    #if script_tree.isLeaf(species_tree,n):
        #script_tree.setName(species_tree,n,script_tree.getName(species_tree,n).split("@")[0])
print script_tree.writeTree(species_tree,Sroot,False)
print "per branch statistics"
print "#species #genes #genesX #genesA #genes_alone #scaffolds(>=1) #scaffoldsX #scaffoldsA #SCJ #SCJX #SCJA #dup"
scj_X = 0.0
scj_A = 0.0
adj_X = 0.0
adj_A = 0.0
adj_extant = 0.0
for n in statistics_per_node.keys():
    if n in extant_species:
        print "E"+n,
        print statistics_per_node[n]["#genes"],
        print statistics_per_node[n]["#genesX"],
        print statistics_per_node[n]["#genesA"],
        print statistics_per_node[n]["#genes_alone"],
        print statistics_per_node[n]["#scaffolds"],
        print statistics_per_node[n]["#scaffoldsX"],
        print statistics_per_node[n]["#scaffoldsA"],
        print statistics_per_node[n]["#SCJ"],
        print statistics_per_node[n]["#SCJX"],
        print statistics_per_node[n]["#SCJA"],
        print dup_per_branch[n]
        scj_X = scj_X +  statistics_per_node[n]["#SCJX"]
        scj_A = scj_A +  statistics_per_node[n]["#SCJA"]
        adj_X = adj_X +  (statistics_per_node[n]["#genesX"]-statistics_per_node[n]["#scaffoldsX"])
        adj_A = adj_A +  (statistics_per_node[n]["#genesA"]-statistics_per_node[n]["#scaffoldsA"])
        adj_extant = adj_extant + (statistics_per_node[n]["#genesX"]-statistics_per_node[n]["#scaffoldsX"]) + (statistics_per_node[n]["#genesA"]-statistics_per_node[n]["#scaffoldsA"])
for n in statistics_per_node.keys():
    if not n in extant_species:
        print "A"+n,
        print statistics_per_node[n]["#genes"],
        print statistics_per_node[n]["#genesX"],
        print statistics_per_node[n]["#genesA"],
        print statistics_per_node[n]["#genes_alone"],
        print statistics_per_node[n]["#scaffolds"],
        print statistics_per_node[n]["#scaffoldsX"],
        print statistics_per_node[n]["#scaffoldsA"],
        print statistics_per_node[n]["#SCJ"],
        print statistics_per_node[n]["#SCJX"],
        print statistics_per_node[n]["#SCJA"],
        print dup_per_branch[n]
        scj_X = scj_X +  statistics_per_node[n]["#SCJX"]
        scj_A = scj_A +  statistics_per_node[n]["#SCJA"]
        adj_X = adj_X +  (statistics_per_node[n]["#genesX"]-statistics_per_node[n]["#scaffoldsX"])
        adj_A = adj_A +  (statistics_per_node[n]["#genesA"]-statistics_per_node[n]["#scaffoldsA"])
   
print
print "global statistics"
print "adj adjE adjX adjA SCJ SCJX fraction SCJ_rate SCJ_rate_X SCJ_rate_A Loss LossX LossA LossNA Dup DupX DupA DupNA XtoA AtoX XtoAu AtoXu"
print int(adj_X+adj_A),int(adj_extant),int(adj_X),int(adj_A),
print total,total_X,total_X/total,
#total+Ptotal,total_X+Ptotal_X,(total_X+Ptotal_X)/(total+Ptotal),
print (scj_X+scj_A)/(adj_X+adj_A),scj_X/adj_X,scj_A/adj_A,
print events["Loss"]["X"]+events["Loss"]["A"]+events["Loss"]["NA"],events["Loss"]["X"],events["Loss"]["A"],events["Loss"]["NA"],
print events["Dup"]["X"]+events["Dup"]["A"]+events["Dup"]["NA"],events["Dup"]["X"],events["Dup"]["A"],events["Dup"]["NA"],
print moves["XtoA"],moves["AtoX"],moves["XtoAu"],moves["AtoXu"]

    
