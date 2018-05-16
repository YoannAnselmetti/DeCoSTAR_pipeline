__author__="Cedric Chauve"
__email__="cedric.chauve@sfu.ca"
__date__="December 18, 2016"

import sys

# ---------------------------------------------------------------------------
USAGE="python filter_scaffolds_by_outgroups.py <scaffolds_file> <outgroups_file> <output_file>\n"
USAGE+="  scaffolds_file: scaffolds\n"
USAGE+="                  format = comments prefixed by #\n"
USAGE+="                  format = one line per gene: species_id scaffold_id gene_name orientation(+/-)\n"
USAGE+="  outgroups_file: pairs of genes in extant or pre-sepciation branches\n"
USAGE+="                     format = 1 line per pair of ancestor-descendant gene: species1_id species2_id gene1_name gene2_name family_id\n"
USAGE+="                     species1 = ancestor, species2 = descendant, gene names as in the gene tree file\n"
USAGE+="                     output branches are only the pre-extant or pre-speciation branches\n"
USAGE+="  output_file: scaffold file that keeps the names of existing scaffolds and remove genes that are not in the orthogroups\n"
USAGE+="               scaffolds that become empty are discarded\n"

# ---------------------------------------------------------------------------
# Reading parameters
file_name_scaffolds   = sys.argv[1]
file_name_orthogroups = sys.argv[2]
output_file_name      = sys.argv[3]
    
# ---------------------------------------------------------------------------
# Reading the scaffold file to create an index of genes that are all, by default assumed to be absent at first

scaffolds_stream = open(file_name_scaffolds,"r").readlines()
all_genes = {}
for line in scaffolds_stream:
    if line[0]!="#" and line[0]!=">":
        line1=line.split()
        (species,gene)=(line1[0],line1[2])
        all_genes[(species,gene)]=False

# ---------------------------------------------------------------------------
# Reading the orthogroups to update the all_genes array

orthogroups_stream = open(file_name_orthogroups,"r").readlines()
for line in orthogroups_stream:
    if line[0]!="#":
        line1=line.split()
        (species1,gene1)=(line1[0],line1[2])
        (species2,gene2)=(line1[1],line1[3])
        all_genes[(species1,gene1)]=True
        all_genes[(species2,gene2)]=True    

# ---------------------------------------------------------------------------
# Reading the scaffold file again to output only the genes seen in the orthogroups

scaffolds_stream = open(file_name_scaffolds,"r").readlines()
output_stream    = open(output_file_name,"w")

for line in scaffolds_stream:
    if line[0]!="#" and line[0]!=">":
        line1=line.split()
        (species,gene)=(line1[0],line1[2])
        if all_genes[(species,gene)]==True:
            output_stream.write(line)
    else:
        output_stream.write(line)
output_stream.close()
