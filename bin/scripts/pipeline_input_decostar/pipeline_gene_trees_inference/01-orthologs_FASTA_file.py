#!/usr/bin/env python
# -*- coding: utf-8 -*-
###
###   Goal:
###      Take file with gene families OR tree and directory containing FASTA files with CDS and produce 1 file / gene family                   
###                                                                        
###   INPUT:                                                               
###      1- Gene families file                                             
###         (data/GENE_TREES/unrooted_trees-ASTEI+filt.nwk)
###         (data/INPUT_DATA/ODBMOZ2_Anophelinae_tabtext)
###      2- INPUT directory path where CDS FASTA files are available       
###         (data/INPUT_DATA/FASTA/CDS)             
###      3- OUTPUT directory path where GF FASTA file will be stored       
###         (data/FASTA/GT_FASTA/CDS/Rob_Gene_TREES/)    
###         (data/FASTA/GT_FASTA/CDS/Rob_Gene_FAMILIES-ASTEI/)
###      4- OUTPUT directory path where gene trees of size 2 or 3 will be stored
###         (data/GENE_TREES/CDS/bootstrap_support/Rob_Gene_TREES/profileNJ/X_topo/PROC_GENE_TREES_100)
###         (data/GENE_TREES/CDS/bootstrap_support/Rob_Gene_TREES/profileNJ/WG_topo/PROC_GENE_TREES_100)
###      5- gene_TAG-species_name association file                         
###         (data/INPUT_DATA/name_geneID_18anopheles)    
###      6- prefix/postfix boolean                                         
###         (prefix or postfix)                                            
###      7- Species tree file                                              
###         (data/INPUT_DATA/Anopheles_species_tree_X_topology.nwk)
###         (data/INPUT_DATA/Anopheles_species_tree_WG_topology.nwk)
###      8- Directory containing RAW gene trees of Rob Waterhouse (to name gene trees)
###         (data/INPUT_DATA/OG_CDS_newtrees)            
###      9- Character separator between species name and gene ID           
###         (@)                                                            
###                                                                        
###   OUTPUT:                                                              
###      - OUTPUT directory containing 1 FASTA file/Gene family            
###                                                                        
###   Name: 01-orthologs_FASTA_file.py        Author: Yoann Anselmetti     
###   Creation date: 2016/03/09               Last modification: 2016/11/28
###


from sys import argv, path, exit
from re import search, sub, match
from os import close, path, makedirs, rename, listdir
from shutil import rmtree
from datetime import datetime
import errno
import subprocess
from ete3 import Tree

def mkdir_p(dir_path):
   try:
      makedirs(dir_path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST and path.isdir(dir_path):
         pass
      else: raise

def replaceStringInFile(filePath,old_pattern,new_pattern):
   "replaces all string by a regex substitution"
   tempName = filePath+'~~~'
   inputFile = open(filePath)
   outputFile = open(tempName,'w')
   fContent = unicode(inputFile.read(), "utf-8")

   outText = sub(old_pattern, new_pattern, fContent)
   outputFile.write((outText.encode("utf-8")))

   outputFile.close()
   inputFile.close()

   rename(tempName, filePath)

# def root_gene_tree(tree_file,notung_path,species_tree):
#    output_dir=path.dirname(path.abspath(tree_file))
#    command_line="java -jar "+notung_path+" --speciestag nhx -g "+tree_file+" -s "+species_tree+" --root --costdup 2 --costloss 1 --nolosses --maxtrees 5 --allopt --treeoutput newick --outputdir "+output_dir+"; mv "+tree_file+".rooting.0 "+tree_file
#    print "\t"+command_line
#    process = subprocess.Popen(command_line.split(), stdout=subprocess.PIPE)
#    output, error = process.communicate()



def add_species_name_to_geneID(gene,dict_geneID_speciesName,order_bool):
   species=""
   # Find species corresponding to the gene
   for geneID in dict_geneID_speciesName:
      if geneID in gene:
         print "\t"+gene
         species=dict_geneID_speciesName[geneID]
         break
   # Add species name to gene ID with separator
   spe_gene=""
   if order_bool=="prefix":
      spe_gene=species+sep+gene
   elif order_bool=="postfix":
      spe_gene=gene+sep+species
   else:
      exit("ERROR, parameter 6 should be equal to \"postfix\" or \"prefix\" !!!")

   return spe_gene


def write_output(list_genes,OUTPUT_dir_MSA,GF_ID,OUTPUT_dir_GT,dict_geneID_speciesName,order_bool,sep,species_tree,gf_size_one):
   # If gene families >= 3:
   if len(list_genes)>3:
      # Open FASTA file for the current gene tree
      output_file=open(OUTPUT_dir_MSA+"/"+GF_ID,'w')   
      # print "(GFsize>3):"
      # Browse list_genes to create FASTA file of the current gene tree
      for gene in list_genes:
         if gene in dict_ID_seq:
            # print "\t"+gene
            output_file.write(dict_ID_seq[gene]+"\n")
      output_file.close()

   # If gene families size == 2 or 3:
   elif (len(list_genes)==2 or len(list_genes)==3):
      tree_str=""
      if len(list_genes)==2:
         # print "(GFsize==2):"
         gene1=add_species_name_to_geneID(list_genes[0],dict_geneID_speciesName,order_bool)
         gene2=add_species_name_to_geneID(list_genes[1],dict_geneID_speciesName,order_bool)
         tree_str="("+gene1+","+gene2+");"
      else:
         # print "(GFsize==3):"
         gene1=add_species_name_to_geneID(list_genes[0],dict_geneID_speciesName,order_bool)
         gene2=add_species_name_to_geneID(list_genes[1],dict_geneID_speciesName,order_bool)
         gene3=add_species_name_to_geneID(list_genes[2],dict_geneID_speciesName,order_bool)
         tree_str="("+gene1+","+gene2+","+gene3+");"

      # Load current gene family in Tree object "tree" and write it in OUTPUT_tree file
      # print tree_str
      tree=Tree(tree_str)
      # print tree
      OUTPUT_tree=OUTPUT_dir_GT+"/"+GF_ID+".nwk"
      tree.write(format=9,outfile=OUTPUT_tree)


   else:
      gf_size_one.append(GF_ID)
      print "Gene family "+GF_ID+" is empty or has a size of 1 gene!!!"


if __name__ == '__main__':

   start_time = datetime.now()

   # Recovery of input parameters
   GF_file=argv[1]
   INPUT_dir_CDS=argv[2]
   OUTPUT_dir_MSA=argv[3]
   OUTPUT_dir_GT=argv[4]
   speciesName_geneID_file=open(argv[5],"r")
   order_bool=argv[6]
   species_tree=argv[7]
   GT_dir=argv[8]
   sep=argv[9]

   check_gene_in_multi_GF=False


   if not (order_bool=="prefix" or order_bool=="postfix"):
      exit("ERROR, parameter 6 should be equal to \"prefix\" or \"postfix\" and not: "+order_bool+"")

   # To be sure than directory have no "/" to the end of the path
   INPUT_dir_CDS=path.normpath(INPUT_dir_CDS)

   # To be sure than directory have no "/" to the end of the path
   OUTPUT_dir_MSA=path.normpath(OUTPUT_dir_MSA)

   # To be sure than directory have no "/" to the end of the path
   OUTPUT_dir_GT=path.normpath(OUTPUT_dir_GT)

   # Remove OUTPUT_dir_MSA if existing
   if not path.exists(OUTPUT_dir_MSA):
      mkdir_p(OUTPUT_dir_MSA)

   # Remove OUTPUT_dir_MSA if existing
   if not path.exists(OUTPUT_dir_GT):
      mkdir_p(OUTPUT_dir_GT)



   # Get list of file of RAW gene trees   
   raw_gene_trees = [f for f in listdir(GT_dir) if path.isfile(path.join(GT_dir,f))]

   # Store association gene ID with gene family ID in "dict_geneID_gfID"
   dict_geneID_gfID={}
   for gt in raw_gene_trees:
      tree_file=GT_dir+"/"+gt
      str_tree=""

      raw_tree_file=open(tree_file,'r')
      str_tree=raw_tree_file.read()

      if str_tree!="();":
         gfID=gt.split(".")[1]
         tree=Tree(tree_file)
         # Get list of extant genes in current gene tree
         for gene in tree.get_leaf_names():
            dict_geneID_gfID[gene]=gfID



   # Store association Gene-Species in "dict_gene_species"
   dict_geneID_speciesName={}
   for gene in speciesName_geneID_file:
      r=search("^([^\t ]*)[\t ]*([^\t ]*)\n$",gene)
      if r:
         speciesName=r.group(1)
         geneID=r.group(2)
         dict_geneID_speciesName[geneID]=speciesName
   speciesName_geneID_file.close()


##########################################################################
### GET ALL CDS SEQUENCE AND STORE IT BY GENE_ID, THEN BY SPECIES NAME ###
##########################################################################
   print "Indexation of gene_ID and associated CDS sequence by species name:"
   dict_ID_seq={}
   for CDS in sorted(listdir(INPUT_dir_CDS)):
      species=""
      r_cds=search("^([^-]*)-([^-]*)-.*\.fa$",CDS)
      r_cds2=search("^[A-Z0-9]*\.([^\.]*)\.cds.fas$",CDS)
      r_pep=search("^[A-Z0-9]*\.([^\.]*)\.pep.fas$",CDS)
      if r_cds:
         genus=r_cds.group(1)
         spe=r_cds.group(2)
         species=genus+"_"+spe
         print "\t"+species
      elif r_cds2:
         species=r_cds2.group(1)
         print "\t"+species
      elif r_pep:
         species=r_pep.group(1)
         print "\t"+species
      else:
         exit("ERROR, name of file "+INPUT_dir_CDS+"/"+CDS+" is bad written!!!")

      # Browse CDS FASTA file to fill dict_ID_seq
      with open(INPUT_dir_CDS+"/"+CDS,'r') as FASTA_file:
         gene_ID=""
         sequence=""
         for line in FASTA_file:
            r_ID=search("^>[^ ]* (.*)\n$",line)
            r_ID2=search("^>([^\n]*)\n$",line)
            if r_ID:
               if gene_ID:
                  dict_ID_seq[gene_ID]=sequence
               sequence=line
               gene_ID=r_ID.group(1)
            elif r_ID2:
               if gene_ID:
                  dict_ID_seq[gene_ID]=sequence
               sequence=line
               gene_ID=r_ID2.group(1)
            else:
               r_dna=search("^([ATCGNatcgn]*)\n$",line)
               r_prot=search("^([RHKDESTNQCUGPAVILMFYWXrhkdestnqcugpavilmfywx]*)\n$",line)
               if r_dna:
                  sequence+=r_dna.group(1)
               elif r_prot:
                  sequence+=r_prot.group(1)
               else:
                  exit("ERROR, the line "+line+" is bad written!!!")
      # To take the last gene of the FASTA file
      dict_ID_seq[gene_ID]=sequence



##################################################################################################################################################
### READING GENE FAMILIES/TREES FILE(S) STORE THEM AND WRITE GENE FAMILIES FASTA FILES (SIZE>=3 GENES) OR WIRTE GENE TREE (SIZE==2 OR 3 GENES) ###
##################################################################################################################################################
   print "Reading gene families file, store them and write them."
   bool_file=""
   input_file=open(GF_file,"r")
   list_check_genes=list()
   list_genes_multi=list()
   dict_gfID_geneslist=dict()
   GF_size_one=list()
   i=0
   # Browse gene trees file line by line
   for line in input_file:
      if not bool_file:
         r=search("^\(",line)
         if r:
            bool_file="GT"
            print "File containing gene cluster is a gene TREES file"
         else:
            bool_file="GF"
            print "File containing gene cluster is a gene FAMILIES file"

      if bool_file=="GF":
         r=search("^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t\n]*)\n$",line)
         if r:
            ODBMOZ2_Level=r.group(1)
            ODBMOZ2_OG_ID=r.group(2)
            Protein_id=r.group(3)
            Gene_ID=r.group(4)
            Organism=r.group(5)
            UniProt_Species=r.group(6)
            UniProt_ACC=r.group(7)
            UniProt_Description=r.group(8)
            InterPro_domains=r.group(9)

            if ODBMOZ2_Level!="ODBMOZ2_Level":
               gene_present=False
               # Check if species of the gene is present in "speciesName_geneID_file". Else we have to prune current gene tree to filter this gene
               for geneID in dict_geneID_speciesName:
                  if geneID in Gene_ID:
                     gene_present=True
                     break

               # if not gene_present:
               #    exit("ERROR, species of gene "+Gene_ID+" is not present in "+str(argv[5])+" !!!!")

               if gene_present:
                  # If want to check if gene(s) is/are present in several gene families
                  if check_gene_in_multi_GF:
                     if not Gene_ID in list_check_genes:
                        list_check_genes.append(Gene_ID)
                     else:
                        list_genes_multi.append(Gene_ID)
                        print "!!! "+Gene_ID+" is present in several gene families !!!"

                  if not ODBMOZ2_OG_ID in dict_gfID_geneslist:
                     dict_gfID_geneslist[ODBMOZ2_OG_ID]=list()
                  dict_gfID_geneslist[ODBMOZ2_OG_ID].append(Gene_ID)

      elif bool_file=="GT":
         list_genes=list()
         tree_str=line.replace("\n","")
         # print tree_str
         tree=Tree(tree_str)
         # Get list of extant genes in current gene tree
         for gene in tree.get_leaf_names():
            # If want to check if gene(s) is/are present in several gene families 
            if check_gene_in_multi_GF:
               if not gene in list_check_genes:
                  list_check_genes.append(gene)
               else:
                  list_genes_multi.append(gene)
                  print "!!! "+gene+" is present in several gene families !!!"

            gene_present=False
            # Check if species of the gene is present in "speciesName_geneID_file". Else we have to prune current gene tree to filter this gene
            for geneID in dict_geneID_speciesName:
               if geneID in gene:
                  gene_present=True
                  list_genes.append(gene)
                  break
            if not gene_present:
               exit("ERROR, species of gene "+gene+" is not present in "+str(argv[5])+" !!!!")

         gfID=dict_geneID_gfID[gene]
         write_output(list_genes,OUTPUT_dir_MSA,gfID,OUTPUT_dir_GT,dict_geneID_speciesName,order_bool,sep,species_tree,GF_size_one)
   input_file.close()
   del list_check_genes[:]


   if bool_file=="GF":
      for gfID in dict_gfID_geneslist:
         list_genes=dict_gfID_geneslist[gfID]
         print 
         write_output(list_genes,OUTPUT_dir_MSA,gfID,OUTPUT_dir_GT,dict_geneID_speciesName,order_bool,sep,species_tree,GF_size_one)

   if check_gene_in_multi_GF:
      print "There are "+str(len(list_genes_multi))+" genes that are present in several gene families/trees:"
      for gene in list_genes_multi:
         print "\t"+gene

   print "There are "+str(len(GF_size_one))+" gene families of size <=1 gene !!!:"
   for gf in GF_size_one:
      print "\t"+gf


   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
