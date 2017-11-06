#!/bin/bash
###                                                                    
###   Goal:                                                            
###      Compute best ML-tree with RAxML from MSA compute with MUSCLE  
###      and Gblocks                                                   
###                                                                    
###   INPUT:                                                           
###      1- PATH of "RAxML" executive file                             
###         (Software/RAxML/raxmlHPC-AVX)                              
###      2- INPUT Gblocks MSA file                                     
###         (DATA/DATA_WaterHouse/FASTA/MSA/PROT/Gblocks/GF0000001-gb) 
###      3- OUTPUT directory that will contain ML trees with support   
###         (DATA/DATA_WaterHouse/GENE_TREES/PROT/bootstrap/RAxML)     
###      4- PROTEIN or RNA boolean (Protein: P/p | CDS: C/c)           
###      5- Boolean for SH-like or bootstrap support on ML trees       
###         (SH-like | bootstrap)                                      
###                                                                    
###   OUTPUT:                                                          
###      - 1 directory per gene tree with files produced by RAxML run  
###                                                                    
###   Name: RAxML.sh                      Author: Yoann Anselmetti     
###   Creation date: 2016/01/06           Last modification: 2017/03/29
###                                                                    

startTime=`date +%s`



if [ ! -d "$3" ]; then
	mkdir -p $3
fi

# Get full path of OUTPUT directory
OUTPUT_dir=$(cd $(dirname $3) && pwd -P)/$(basename $3)

# Get Gene Family ID
GF_ID=$(basename $2 | sed 's/-gb//g')
if [ -d "$OUTPUT_dir/$GF_ID" ]; then
	rm -rf $OUTPUT_dir/$GF_ID/*
else
	mkdir -p $OUTPUT_dir/$GF_ID
fi
echo "Processing Gene Family $GF_ID:"

# IF CDS FASTA FILE ARE PROTEIN SEQUENCES
if [[ $4 == "P" || $4 == "p" ]]; then

	# If SH-like support (FAST computation)
   if [[ $5 == "SH-like" ]]; then
		# COMPUTE the BEST ML-tree for the alignment with default algo of RAxML (-f d: new rapid hill­ climbing) and with model of AA substitution with GAMMA distribution and the GTR AA substitution matrix with ML estimation of AA frequencies ("X") => -N 3: test 3 starting trees to find bestML tree
		echo -e "\t$1 -f d -m PROTGAMMAGTRX -p 12345 -N 3 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_tree
		$1 -f d -m PROTGAMMAGTRX -p 12345 -N 3 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_tree"

		# THEN COMPUTE, SH-like support (-f J) from the BEST ML-tree (Designed from PhyML 3.0) with model of AA substitution with GAMMA distribution and the GTR AA substitution matrix with ML estimation of AA frequencies ("X")
		echo -e "\t$1 -f J -m PROTGAMMAGTRX -p 12345 -s $2 -t $OUTPUT_dir/$GF_ID/RAxML_bestTree.$GF_ID"_tree" -w $OUTPUT_dir/$GF_ID -n $GF_ID"_treeSH
		$1 -f J -m PROTGAMMAGTRX -p 12345 -s $2 -t $OUTPUT_dir/$GF_ID/RAxML_bestTree.$GF_ID"_tree" -w $OUTPUT_dir/$GF_ID -n $GF_ID"_treeSH"

	# If boostrap support (CONFIDENT computation)
   elif [[ $5 == "bootstrap" ]]; then
		# COMPUTE the BEST ML-tree for the alignment with default algo of RAxML (-f d: new rapid hill­ climbing) and with model of AA substitution with GAMMA distribution and the GTR AA substitution matrix with ML estimation of AA frequencies ("X") => -N 10: test 10 starting trees to find bestML tree
		echo -e "\t$1 -f d -m PROTGAMMAGTRX -p 12345 -N 10 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_tree
		$1 -f d -m PROTGAMMAGTRX -p 12345 -N 10 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_tree"

		# THEN, turn on boostrapping 100 times (-N 100 -b 12345) with model of AA substitution with GAMMA distribution and the GTR AA substitution matrix with ML estimation of AA frequencies ("X")
		echo -e "\t$1 -f d -m PROTGAMMAGTRX -p 12345 -N 100 -b 12345 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_BS
		$1 -f d -m PROTGAMMAGTRX -p 12345 -N 100 -b 12345 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_BS"

		# THEN, get bootstrap trees (RAxML_bootstrap.$GF_ID"_BS") and compute bootstrap support for each branch of the Best ML trees compute in first command-line
		echo -e "\t$1 -f b -m PROTGAMMAGTRX -p 12345 -t $OUTPUT_dir/$GF_ID/RAxML_bestTree.$GF_ID"_tree" -z $OUTPUT_dir/$GF_ID/RAxML_bootstrap.$GF_ID"_BS" -w $OUTPUT_dir/$GF_ID -n $GF_ID"_treeBS
		$1 -f b -m PROTGAMMAGTRX -p 12345 -t $OUTPUT_dir/$GF_ID/RAxML_bestTree.$GF_ID"_tree" -z $OUTPUT_dir/$GF_ID/RAxML_bootstrap.$GF_ID"_BS" -w $OUTPUT_dir/$GF_ID -n $GF_ID"_treeBS"
   else
		echo -e "\n!!! ERROR: Param 5 should be SH-like (for FAST computation) or bootstrap (for CONFIDENT support) and not $5 !!!"
   	exit
	fi

# IF CDS FASTA FILE ARE DNA/RNA SEQUENCES
elif [[ $4 == "C" || $4 == "c" ]]; then

	# If SH-like support (FAST computation)
   if [[ $5 == "SH-like" ]]; then
		# COMPUTE the BEST ML-tree for the alignment with default algo of RAxML (-f d: new rapid hill­ climbing) and with model of Nt substitution with GAMMA distribution and the GTR Nt substitution matrix with ML estimation of bases frequencies ("X") => -N 3: test 3 starting trees to find bestML tree
		echo -e "\t$1 -f d -m GTRGAMMAX -p 12345 -N 3 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_tree
		$1 -f d -m GTRGAMMAX -p 12345 -N 3 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_tree"

		# THEN COMPUTE, SH-like support (-f J) from the BEST ML-tree (Designed from PhyML 3.0) with model of Nt substitution with GAMMA distribution and the GTR Nt substitution matrix with ML estimation of bases frequencies ("X")
		echo -e "\t$1 -f J -m GTRGAMMAX -p 12345 -s $2 -t $OUTPUT_dir/$GF_ID/RAxML_bestTree.$GF_ID"_tree" -w $OUTPUT_dir/$GF_ID -n $GF_ID"_treeSH
		$1 -f J -m GTRGAMMAX -p 12345 -s $2 -t $OUTPUT_dir/$GF_ID/RAxML_bestTree.$GF_ID"_tree" -w $OUTPUT_dir/$GF_ID -n $GF_ID"_treeSH"

	# If boostrap support (CONFIDENT computation)
   elif [[ $5 == "bootstrap" ]]; then
		# COMPUTE the BEST ML-tree for the alignment with default algo of RAxML (-f d: new rapid hill­ climbing) and with model of Nt substitution with GAMMA distribution and the GTR AA substitution matrix with ML estimation of AA frequencies ("X") => -N 10: test 10 starting trees to find bestML tree
		echo -e "\t$1 -f d -m GTRGAMMAX -p 12345 -N 10 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_tree
		$1 -f d -m GTRGAMMAX -p 12345 -N 10 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_tree"

		# THEN, turn on boostrapping 100 times with model of AA substitution with GAMMA distribution and the GTR Nt substitution matrix with ML estimation of AA frequencies ("X")
		echo -e "\t$1 -f d -m GTRGAMMAX -p 12345 -N 100 -b 12345 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_BS
		$1 -f d -m GTRGAMMAX -p 12345 -N 100 -b 12345 -s $2 -w $OUTPUT_dir/$GF_ID -n $GF_ID"_BS"

		# THEN, get bootstrap trees and compute bootstrap support for each branch of the Best ML trees compute in first command-line
		echo -e "\t$1 -f b -m GTRGAMMAX -p 12345 -t $OUTPUT_dir/$GF_ID/RAxML_bestTree.$GF_ID"_tree" -z $OUTPUT_dir/$GF_ID/RAxML_bootstrap.$GF_ID"_BS" -w $OUTPUT_dir/$GF_ID -n $GF_ID"_treeBS
		$1 -f b -m GTRGAMMAX -p 12345 -t $OUTPUT_dir/$GF_ID/RAxML_bestTree.$GF_ID"_tree" -z $OUTPUT_dir/$GF_ID/RAxML_bootstrap.$GF_ID"_BS" -w $OUTPUT_dir/$GF_ID -n $GF_ID"_treeBS"
   else
		echo -e "\n!!! ERROR: Param 5 should be SH-like (for FAST computation) or bootstrap (for CONFIDENT support) and not $5 !!!"
   	exit
	fi
else
	echo -e "\n!!! ERROR: Param 4 should be p (for Protein sequences) or r (for RNA/DNA sequences) and not $4 !!!"
	exit
fi



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