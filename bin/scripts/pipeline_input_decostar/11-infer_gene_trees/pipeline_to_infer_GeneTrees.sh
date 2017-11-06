#!/bin/bash
###                                                                         
###   Goal:                                                                 
###      Script to obtain gene trees from CDS FASTA file and gene clustering
###      file (Here, gene trees available on WaterHouse website)            
###   INPUT:                                                                
###      1- DIRECTORY containing scripts                                    
###         (bin/scripts/pipeline_input_decostar/11-infer_gene_trees) 
###      2- Gene TREES or FAMILIES file                                     
###         (data/GENE_TREES/unrooted_trees-ASTEI+filt.nwk)
###      3- DIRECTORY containing CDS FASTA files                            
###         (data/INPUT_DATA/FASTA/CDS)               
###      4- OUTPUT DIRECTORY containing GT FASTA files (1 FASTA file contains
###         all CDS from the same gene cluster/gene tree)                   
###         (data/FASTA/GF_FASTA/CDS)                     
###      5- OUTPUT DIRECTORY containing MSA files (From MUSCLE and Gblocks) 
###         (data/FASTA/MSA/CDS)                          
###      6- OUTPUT DIRECTORY containing gene trees files                    
###         (data/GENE_TREES/CDS)                         
###      7- Species-Gene association file                                   
###         (data/INPUT_DATA/name_geneID_18anopheles)   
###      8- Species tree file                                               
###         (data/INPUT_DATA/Anopheles_species_tree_X_topology.nwk)
###         (data/INPUT_DATA/Anopheles_species_tree_WG_topology.nwk)
###      9- Max jobs number allowed
###         (Example: 100)                                                  
###      10- PROTEIN or CDS boolean                                         
###          (C)                                      
###      11- Bootstrap or SH-like support compute with RAxML boolean        
###          (bootstrap)                                          
###      12- Minimum support allow on gene trees (used in profileNJ step)   
###          (100)                                                          
###      13- Species position compare to the separator with gene ID         
###          (prefix)                                            
###      14- Directory containing RAW gene trees of Rob Waterhouse          
###         (data/INPUT_DATA/OG_CDS_newtrees)             
###      15- Character separator between species name and gene ID           
###         (@)                                                             
###      16- Directory for species tree topology used for gene trees inference
###         (X_topo / WG_topo)                                              
###      17- Gene cluster indication for OUTPUT directory name                                   
###         (Rob_Gene_TREES)
###      18- OUTPUT gene trees file
###         (data/GENE_TREES/trees_DECOSTAR_Xtopo.nwk  /  data/GENE_TREES/trees_DECOSTAR_WGtopo.nwk)
###                                                                         
###   OUTPUT:                                                               
###      - DIRECTORY containing gene trees improved with RAxML & profileNJ  
###                                                                         
###   Name: 11-pipeline_to_infer_GeneTrees.sh	Author: Yoann Anselmetti  
###   Creation date: 2016/01/11						Last modification: 2017/03/29
###                                                                         

MUSCLE=$HOME/Software/MSA/muscle3.8.425/muscle3.8.425_i86linux64
GBLOCKS=$HOME/Software/TRIMMING_ALIGN/Gblocks_0.91b/Gblocks
RAXML=$HOME/Software/TREE_INFERENCE/RAxML-8.2.8/raxmlHPC-SSE3 # AVX compilation doesn't work on the cluster => WRONG some nodes have CPUs with AVX or AVX2 compilation available (If want to use them need to give the node during job submission)
PROFILENJ=$HOME/Software/TREE_INFERENCE/profileNJ/bin/profileNJ
# NOTUNG=$HOME/Software/TREE_INFERENCE/Notung-2.8.1.7/Notung-2.8.1.7.jar

boolTopo="X"

#################
### PARAMETERS
#################
scriptDIR="bin/scripts/pipeline_input_decostar/11-infer_gene_trees"
geneTrees="data/GENE_TREES/unrooted_trees-ASTEI+filt.nwk"
CDSfasta="data/INPUT_DATA/FASTA/CDS"
GTfasta="data/FASTA/GF_FASTA/CDS"
MSAfasta="data/FASTA/MSA/CDS"
GTout="data/GENE_TREES/CDS"
nameGeneID="data/INPUT_DATA/name_geneID_18anopheles"
speciesTree=""
if [[ $boolTopo == "X" ]]; then
   speciesTree="data/INPUT_DATA/Anopheles_species_tree_X_topology.nwk"
elif [[ $boolTopo == "WG" ]]; then
   speciesTree="data/INPUT_DATA/Anopheles_species_tree_WG_topology.nwk"
else
   echo -e "\n!!! boolTopo should be equal to \"X\" or \"WG\" and not: $boolTopo !!!"
   exit
fi
JobMax=100
boolSeq="C"                             
boolSupp="bootstrap"
pNJseuil=100
boolGeneID="prefix"
rawGT="data/INPUT_DATA/OG_CDS_newtrees"
sep="@"
speciesTreeTopo=""
if [[ $boolTopo == "X" ]]; then
   speciesTreeTopo="X_topo"
elif [[ $boolTopo == "WG" ]]; then
   speciesTreeTopo="WG_topo"
else
   echo -e "\n!!! boolTopo should be equal to \"X\" or \"WG\" and not: $boolTopo !!!"
   exit
fi
outputTAG="Rob_Gene_TREES"
outputTrees=""
if [[ $boolTopo == "X" ]]; then
   outputTrees="data/GENE_TREES/trees_DeCoSTAR_Xtopo.nwk"
elif [[ $boolTopo == "WG" ]]; then
   outputTrees="data/GENE_TREES/trees_DeCoSTAR_WGtopo.nwk"
else
   echo -e "\n!!! boolTopo should be equal to \"X\" or \"WG\" and not: $boolTopo !!!"
   exit
fi


TREESdir="$GTout/${boolSupp}_support/$outputTAG"
procGTdir="$TREESdir/profileNJ/$speciesTreeTopo/PROC_GENE_TREES_$pNJseuil"
unrootGTdir="$TREESdir/profileNJ/UNROOTED_GENE_TREES"



re='^[0-9]+$'
if ! [[ $JobMax =~ $re ]]; then
   echo "ERROR: Parameter 9 (\#jobs limit) is not a number: $JobMax" >&2;
   exit 1
fi

if [[ $boolSeq != "P" && $boolSeq != "p" && $boolSeq != "C" && $boolSeq != "c" ]]; then
   echo -e "\n!!! ERROR: Parameter 10 for type of sequence used should be \"p\" (for Protein sequences) or \"c\" (for Nucleotide/CDS sequences) and not: $boolSeq" >&2
   exit
fi

if [[ $boolSupp != "bootstrap" && $boolSupp != "SH-like" ]]; then
   echo -e "\n!!! ERROR: Parameter 11 for node support computed with RAxML should be \"SH-like\" (for FAST computation) or \"bootstrap\" (for CONFIDENT support) and not: $boolSupp" >&2
   exit
fi

if ! [[ $pNJseuil =~ $re ]]; then
   echo "ERROR: parameter 12 (minimum support for profileNJ) is not a number: $pNJseuil" >&2;
   exit 1
else
   if (( "$pNJseuil" > "100" )); then 
      echo "ERROR: Threshold value: $pNJseuil > 100" >&2;
      exit 1
   fi
fi

if ! [[ $boolGeneID == "prefix" || $boolGeneID == "postfix" ]]; then
    echo "ERROR: prefix/postfix boolean should be equal to \"prefix\" or \"postfix\" and not: $boolGeneID" >&2;
    exit 1
fi






if [ -d $procGTdir ]; then
   rm -rf $procGTdir/*
else
   mkdir -p $procGTdir
fi


# 1- Write CDS seq from a same gene cluster/tree/family in 1 FASTA file with "orthologs_FASTA_file.py" script
if [ -d "$GTfasta/$outputTAG" ]; then
   rm -rf $GTfasta/$outputTAG/*
else
   mkdir -p $GTfasta/$outputTAG
fi
echo ""
echo "#############################################"
echo "### STEP 1: Write Gene Family FASTA files ###"
echo "#############################################"
echo ""
echo "$scriptDIR/01-orthologs_FASTA_file.py $geneTrees $CDSfasta $GTfasta/$outputTAG $procGTdir $nameGeneID $boolGeneID $speciesTree $rawGT $sep"

$scriptDIR/01-orthologs_FASTA_file.py $geneTrees $CDSfasta $GTfasta/$outputTAG $procGTdir $nameGeneID $boolGeneID $speciesTree $rawGT $sep







# 2- Align CDS seq from the same GT_FASTA file for all GT_FASTA files with MUSCLE
if [ -d "$MSAfasta/$outputTAG/MUSCLE" ]; then
   rm -rf $MSAfasta/$outputTAG/MUSCLE/*
else
   mkdir -p $MSAfasta/$outputTAG/MUSCLE
fi
echo ""
echo "#######################################################"
echo "### STEP 2: Multiple Sequence Alignment with Muscle ###"
echo "#######################################################"
echo ""

startTime=`date +%s`

list_GF=(`ls $GTfasta/$outputTAG`)
nb_tot=$(ls $GTfasta/$outputTAG | wc -l)
x=0
while [ $x -lt $nb_tot ];  # until the ending one
do
   jobs_nb=$(qstat | grep $USER | wc -l)
   if [ $jobs_nb -le $JobMax ];
   then
      echo "qsub -V -N \"muscle_yanselmetti\" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y \"$scriptDIR/02-MUSCLE.sh $MUSCLE $GTfasta/$outputTAG/${list_GF[$x]} $MSAfasta/$outputTAG/MUSCLE\""
      qsub -V -N "muscle_yanselmetti" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y "$scriptDIR/02-MUSCLE.sh $MUSCLE $GTfasta/$outputTAG/${list_GF[$x]} $MSAfasta/$outputTAG/MUSCLE"
      x=$(($x+1))
      sleep 0.2
   else
      echo "No spot yet"
      sleep 10
   fi
done

# Loop until MUSCLE jobs are over
muscle_nb=$(qstat | grep $USER | grep "muscle_yan" | wc -l)
while [ $muscle_nb -ne 0 ];
do
   sleep 10
   muscle_nb=$(qstat | grep $USER | grep "muscle_yan" | wc -l)
done

# Compute excution time of Step 2
endTime=`date +%s`
runtime=$(($endTime-$startTime))

secondes=""
minutes=""
hours=$(($runtime/3600))
modH=$(($runtime%3600))
if [ $modH == 0 ]; then
    minutes="00"
    secondes="00"
else
   minutes=$(($modH/60))
   modM=$(($modH%60))
    if [ $modM == 0 ]; then
       secondes="00"
   else
      secondes=$modM
   fi
fi

echo "Execution time: "$hours":"$minutes":"$secondes








# 3- Restrict MSA to blocks select by Gblocks
if [ -d "$MSAfasta/$outputTAG/Gblocks" ]; then
   rm -rf $MSAfasta/$outputTAG/Gblocks/*
else
   mkdir -p $MSAfasta/$outputTAG/Gblocks
fi
echo ""
echo "######################################################################"
echo "### STEP 3: Select Multiple Sequence Alignment blocks with Gblocks ###"
echo "######################################################################"
echo ""

startTime=`date +%s`

# Get list of MSA files in MUSCLE repertory to apply Gblocks on them
list_GF=(`ls $MSAfasta/$outputTAG/MUSCLE`)
nb_tot=$(ls $MSAfasta/$outputTAG/MUSCLE | wc -l)
x=0
while [ $x -lt $nb_tot ];  # until the ending one
do
   jobs_nb=$(qstat | grep $USER | wc -l)
   if [ $jobs_nb -le $JobMax ];
   then
      echo "qsub -V -N \"gblocks_yanselmetti\" -hold_jid \"muscle_yanselmetti\" --V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y \"$scriptDIR/03-Gblocks.sh $GBLOCKS $MSAfasta/$outputTAG/MUSCLE/${list_GF[$x]} $MSAfasta/$outputTAG/Gblocks $boolSeq\""
      qsub -V -N "gblocks_yanselmetti" -hold_jid "muscle_yanselmetti" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y "$scriptDIR/03-Gblocks.sh $GBLOCKS $MSAfasta/$outputTAG/MUSCLE/${list_GF[$x]} $MSAfasta/$outputTAG/Gblocks $boolSeq"
      x=$(($x+1))
      sleep 0.2
   else
      echo "No spot yet"
      sleep 10
   fi
done

# Loop until Gblocks jobs are over
gblocks_nb=$(qstat | grep $USER | grep "gblocks_ya" | wc -l)
while [ $gblocks_nb -ne 0 ];
do
   sleep 10
   gblocks_nb=$(qstat | grep $USER | grep "gblocks_ya" | wc -l)
done

# Compute excution time of Step 3
endTime=`date +%s`
runtime=$(($endTime-$startTime))

secondes=""
minutes=""
hours=$(($runtime/3600))
modH=$(($runtime%3600))
if [ $modH == 0 ]; then
    minutes="00"
    secondes="00"
else
   minutes=$(($modH/60))
   modM=$(($modH%60))
    if [ $modM == 0 ]; then
       secondes="00"
   else
      secondes=$modM
   fi
fi

echo "Execution time: "$hours":"$minutes":"$secondes








# 4- Infer gene trees from MSA (Gblocks) with RAxML
if [ -d "$TREESdir/RAxML" ]; then
   rm -rf $TREESdir/RAxML/*
else
   mkdir -p $TREESdir/RAxML
fi
echo ""
echo "##########################################"
echo "## STEP 4: Infer Gene Trees with RAxML ###"
echo "##########################################"
echo ""

startTime=`date +%s`

# Get list of MSA files in Gblocks repertory to compute Best ML RAxML tree on them
list_GF=(`ls $MSAfasta/$outputTAG/Gblocks`)
nb_tot=$(ls $MSAfasta/$outputTAG/Gblocks | wc -l)
x=0
while [ $x -lt $nb_tot ];  # until the ending one
do
   jobs_nb=$(qstat | grep $USER | wc -l)
   if [ $jobs_nb -le $JobMax ];
   then
      echo "qsub -V -N \"raxml_yanselmetti\" -hold_jid \"gblocks_yanselmetti\" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y \"$scriptDIR/04-RAxML.sh $RAXML $MSAfasta/$outputTAG/Gblocks/${list_GF[$x]} $TREESdir/RAxML\" $boolSeq $boolSupp"
      qsub -V -N "raxml_yanselmetti" -hold_jid "gblocks_yanselmetti" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y "$scriptDIR/04-RAxML.sh $RAXML $MSAfasta/$outputTAG/Gblocks/${list_GF[$x]} $TREESdir/RAxML $boolSeq $boolSupp"
      x=$(($x+1))
      sleep 0.2
   else
      echo "No spot yet"
      sleep 10
   fi
done


# Remove ".reduced" files created by RAxML
echo -e "\tqsub -V -N \"rm_reduced\" -hold_jid \"raxml_yanselmetti\" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y \"rm $MSAfasta/$outputTAG/Gblocks/*.reduced\""
qsub -V -N "rm_reduced" -hold_jid "raxml_yanselmetti" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y "rm $MSAfasta/$outputTAG/Gblocks/*.reduced"



# Loop until RAxML jobs are over
raxml_nb=$(qstat | grep $USER | grep "raxml_yans" | wc -l)
while [ $raxml_nb -ne 0 ];
do
   sleep 10
   raxml_nb=$(qstat | grep $USER | grep "raxml_yans" | wc -l)
done

# Compute excution time of Step 4
endTime=`date +%s`
runtime=$(($endTime-$startTime))

secondes=""
minutes=""
hours=$(($runtime/3600))
modH=$(($runtime%3600))
if [ $modH == 0 ]; then
    minutes="00"
    secondes="00"
else
   minutes=$(($modH/60))
   modM=$(($modH%60))
    if [ $modM == 0 ]; then
       secondes="00"
   else
      secondes=$modM
   fi
fi

echo "Execution time: "$hours":"$minutes":"$secondes








# 5- Process gene trees with profileNJ to contract weak branches and resolve multifurcation produced by min(GDup and GLos)
if [ -d $unrootGTdir ]; then
   rm -rf $unrootGTdir/*
else
   mkdir -p $unrootGTdir
fi

echo ""
echo "################################################"
echo "### STEP 5: Refine Gene trees with profileNJ ###"
echo "################################################"
   
# 5a- Compute distance matrix for gene trees give as INPUT of profileNJ
echo ""
echo -e "\t################################################################################"
echo -e "\t### STEP 5a: Modify gene trees inferred by RAxML to be readable by profileNJ ###"
echo -e "\t################################################################################"
echo ""

startTime=`date +%s`

# Get list of Gene Family ID repertory containing Best ML RAxML tree to compute distance matrix on them
list_GF=(`ls $TREESdir/RAxML`)
nb_tot=$(ls $TREESdir/RAxML | wc -l)
x=0
while [ $x -lt $nb_tot ];  # until the ending one
do
   jobs_nb=$(qstat | grep $USER | wc -l)
   if [ $jobs_nb -le $JobMax ];
   then
      echo ${list_GF[$x]}
      echo -e "\tqsub -V -N \"modifGT_yanselmetti\" -hold_jid \"raxml_yanselmetti\" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y \"$scriptDIR/05a-modifGT.py $TREESdir/RAxML/${list_GF[$x]} $unrootGTdir $boolSupp $nameGeneID $boolGeneID\""
      qsub -V -N "modifGT_yanselmetti" -hold_jid "raxml_yanselmetti" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y "$scriptDIR/05a-modifGT.py $TREESdir/RAxML/${list_GF[$x]} $unrootGTdir $boolSupp $nameGeneID $boolGeneID"
      x=$(($x+1))
      sleep 0.2
   else
      echo "No spot yet"
      sleep 10
   fi
done

# Loop until Gene trees modification jobs are over
modifGT_nb=$(qstat | grep $USER | grep "modifGT_ya" | wc -l)
while [ ${modifGT_nb} -ne 0 ];
do
   sleep 10
   modifGT_nb=$(qstat | grep $USER | grep "modifGT_ya" | wc -l)
done



# 5b- Contract weak branches and resolve multifurcation produced by min(GDup and GLos number) with profileNJ
echo ""
echo -e "\t##################################################################################"
echo -e "\t### STEP 5b: Refine Gene Trees with profileNJ by contracting weak branches and ###"
echo -e "\t### resolving multifurcation produced by minimizing Gene Duplication and Loss  ###"
echo -e "\t##################################################################################"
echo ""

# Get list of gene trees in UNROOTED_GENE_TREES repertory to apply profileNJ on them
list_GF=(`ls $unrootGTdir`)
nb_tot=$(ls $unrootGTdir | wc -l)
x=0
while [ $x -lt $nb_tot ];  # until the ending one
do
   jobs_nb=$(qstat | grep $USER | wc -l)
   if [ $jobs_nb -le $JobMax ];
   then
      echo "qsub -V -N \"profileNJ_yanselmetti\" -hold_jid \"modifGT_yanselmetti\" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y \"$scriptDIR/05b-profileNJ.sh $PROFILENJ $speciesTree $unrootGTdir/${list_GF[$x]} $procGTdir $pNJseuil $boolGeneID\""
      qsub -V -N "profileNJ_yanselmetti" -hold_jid "modifGT_yanselmetti" -V -cwd -l h_rt=1000:00:00,s_rt=1000:00:00 -b y "$scriptDIR/05b-profileNJ.sh $PROFILENJ $speciesTree $unrootGTdir/${list_GF[$x]} $procGTdir $pNJseuil $boolGeneID"
      x=$(($x+1))
      sleep 0.2
   else
      echo "No spot yet"
      sleep 10
   fi
done


profileNJ_nb=$(qstat | grep $USER | grep "profileNJ_ya" | wc -l)
while [ ${profileNJ_nb} -ne 0 ];
do
   sleep 10
   profileNJ_nb=$(qstat | grep $USER | grep "profileNJ_ya" | wc -l)
done

sleep 10


#####
### Get the last solution for all gene trees and store it in resulting file
tail -n1 $procGTdir/* |grep \( > $outputTrees


# Compute excution time of Step 5
endTime=`date +%s`
runtime=$(($endTime-$startTime))

secondes=""
minutes=""
hours=$(($runtime/3600))
modH=$(($runtime%3600))
if [ $modH == 0 ]; then
    minutes="00"
    secondes="00"
else
   minutes=$(($modH/60))
   modM=$(($modH%60))
    if [ $modM == 0 ]; then
       secondes="00"
   else
      secondes=$modM
   fi
fi

echo "Execution time: "$hours":"$minutes":"$secondes