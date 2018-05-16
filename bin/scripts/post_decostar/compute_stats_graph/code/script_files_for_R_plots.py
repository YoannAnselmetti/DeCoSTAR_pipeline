#! /usr/bin/env python
# -*- coding: utf-8 -*-


from os import close, path, makedirs
from re import search
import errno

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise

def get_SCJ_STATS(INPUT,tag,output,column,bool_tag):
	input_file=open(INPUT,"r")
	for line in input_file:
		regexp=search("^[AE][0-9][0-9]? .*\n",line)
		if regexp:
			speID=line.split(" ")[0]
			TAG=tag
			if bool_tag:
				if "A" in speID:
					TAG=tag+" anc"
				elif "E" in speID:
					TAG=tag+" ext"
				else:
					exit("ERROR, species ID in file  \""+INPUT+"\"  should begin by \"A\" or \"E\" !!!")

			stats=line.split(" ")[int(column)-1]
			output.write(TAG+"\t"+stats+"\n")

	input_file.close()



Xmax=5

PREFIX_Xtopo_RAW = "results/decostar/Xtopo_RAW/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1"
PREFIX_Xtopo_withSCAFF = "results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1"
PREFIX_Xtopo_withoutSCAFF = "results/decostar/Xtopo-scaff/DeCoSTAR_Anopheles_Xtopo-scaff_Boltz_0.1"
PREFIX_WGtopo_withSCAFF = "results/decostar/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1"
outputDIR = "figures/R_plots/files"

mkdir_p(outputDIR)

extant_genes = {}
ancestral_raw = {}
ancestral_pnj = {}

file_genes = open(PREFIX_Xtopo_RAW+".genes.txt","r").readlines()
for line in file_genes:
    words = line.split()
    species = words[0]
    ident = words[1]
    if len(words) == 2:
        if not extant_genes.has_key(species):
            extant_genes[species] = 0
        extant_genes[species] = extant_genes[species] + 1
    else:
        if not ancestral_raw.has_key(species):
            ancestral_raw[species] = 0
        ancestral_raw[species] = ancestral_raw[species] + 1

file_genes = open(PREFIX_Xtopo_withSCAFF+".genes.txt","r").readlines()
for line in file_genes:
    words = line.split()
    species = words[0]
    ident = words[1]
    if len(words) > 2:
        if not ancestral_pnj.has_key(species):
            ancestral_pnj[species] = 0
        ancestral_pnj[species] = ancestral_pnj[species] + 1

#print len(extant_genes.keys()),len(ancestral_pnj.keys()),len(ancestral_raw.keys())

ekeys = extant_genes.keys()
rkeys = ancestral_raw.keys()
pkeys = ancestral_pnj.keys()

output = open(outputDIR+"/content.txt","w")
for i in range(len(ekeys)):
    output.write(str(extant_genes[ekeys[i]])+" ")
    if i<len(rkeys):
        output.write(str(ancestral_raw[rkeys[i]])+" ")
    else:
        output.write("NA ")
    if i<len(pkeys):
        output.write(str(ancestral_pnj[pkeys[i]])+"\n")
    else:
        output.write("NA\n")
        
#print "extant"
#for s in ekeys:
    #print s,extant_genes[s]
    
#print "raw"
#for s in rkeys:
    #print s,ancestral_raw[s]
 
#print "pnj"
#for s in pkeys:
    #print s,ancestral_pnj[s]
    
output = open(outputDIR+"/degrees.txt","w")
degree = {"PNJ":{},"RAW":{}}
file_adjacencies = open(PREFIX_Xtopo_withSCAFF+".adjacencies.txt","r").readlines()
for line in file_adjacencies:
    words = line.split()
    species = words[0]
    gene1 = words[1]
    gene2 = words[2]
    score = float(words[6])
    if line.find("@") < 0:
        if not degree["PNJ"].has_key(gene1):
            degree["PNJ"][gene1] = 0.0
        if not degree["PNJ"].has_key(gene2):
            degree["PNJ"][gene2] = 0.0
        degree["PNJ"][gene1] = degree["PNJ"][gene1] + score
        degree["PNJ"][gene2] = degree["PNJ"][gene2] + score
file_adjacencies = open(PREFIX_Xtopo_RAW+".adjacencies.txt","r").readlines()
for line in file_adjacencies:
    words = line.split()
    species = words[0]
    gene1 = words[1]
    gene2 = words[2]
    score = float(words[6])
    if line.find("@") < 0:
        if not degree["RAW"].has_key(gene1):
            degree["RAW"][gene1] = 0.0
        if not degree["RAW"].has_key(gene2):
            degree["RAW"][gene2] = 0.0
        degree["RAW"][gene1] = degree["RAW"][gene1] + score
        degree["RAW"][gene2] = degree["RAW"][gene2] + score

max_value = int(max(max(degree["PNJ"].values()),max(degree["RAW"].values())))
nb_categories = int(max_value)
tab = {"PNJ":[0]*(nb_categories+1),"RAW":[0]*(nb_categories+1)}
for gene in degree["RAW"].keys():
    value = int(nb_categories*degree["RAW"][gene]/max_value)
    #print value
    tab["RAW"][value] = tab["RAW"][value] + 1
for gene in degree["PNJ"].keys():
    value = int(nb_categories*degree["PNJ"][gene]/max_value)
    #print degree["PNJ"][gene],value
    tab["PNJ"][value] = tab["PNJ"][value] + 1
for i in range(nb_categories):
    if i<=Xmax:
        if i==2:
            output.write(str(i*max_value/nb_categories)+" "+str(tab["RAW"][i]/float(sum(ancestral_raw.values())))+" "+str(tab["PNJ"][i]/float(sum(ancestral_pnj.values())))+" 1\n")
        else:
            output.write(str(i*max_value/nb_categories)+" "+str(tab["RAW"][i]/float(sum(ancestral_raw.values())))+" "+str(tab["PNJ"][i]/float(sum(ancestral_pnj.values())))+" 0\n")





scj_stats_Xtopo_withSCAFF=PREFIX_Xtopo_withSCAFF+"_0.1_M1_scj_stats"
scj_stats_Xtopo_withoutSCAFF=PREFIX_Xtopo_withoutSCAFF+"_0.1_M1_scj_stats"
scj_stats_WGtopo_withSCAFF=PREFIX_WGtopo_withSCAFF+"_0.1_M1_scj_stats"

# get scaffolds number for each species to produce boxplot on scaffolds number
output = open(outputDIR+"/boxplot_scaffolds.txt","w")
column = 6

tag="XNS"
get_SCJ_STATS(scj_stats_Xtopo_withoutSCAFF,tag,output,column,True)
tag="X"
get_SCJ_STATS(scj_stats_Xtopo_withSCAFF,tag,output,column,True)
tag="WG"
get_SCJ_STATS(scj_stats_WGtopo_withSCAFF,tag,output,column,True)

output.close()


# get rearrangements/SCJ number for each species to produce boxplot on rearrangements/SCJ number
output = open(outputDIR+"/boxplot_rearrangements.txt","w")
column = 9

tag="X phylogeny, non scaffold"
get_SCJ_STATS(scj_stats_Xtopo_withoutSCAFF,tag,output,column,False)
tag="X phylogeny, scaffold"
get_SCJ_STATS(scj_stats_Xtopo_withSCAFF,tag,output,column,False)

output.close()