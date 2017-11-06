#!/bin/bash
###                                                                    
###   Goal:                                                            
###      Multiple Sequence Alignment of CDS sequence of genes from same gene trees with MUSCLE                              
###                                                                    
###   INPUT:                                                           
###      1- PATH of "MUSCLE" executive file                            
###         (Software/MSA/muscle3.8.425/muscle3.8.425_i86linux64)      
###      2- INPUT FASTA file with CDS sequence (PROT/CDS) from a same Gene family                                                
###         (DECOSTAR/DATA_DeCoSTAR/FASTA/GT_FASTA/CDS/MZ22500058)     
###      3- OUTPUT directory that will contain MSA output files of MUSCLE
###         (DECOSTAR/DATA_DeCoSTAR/FASTA/MSA/CDS/MUSCLE/)             
###                                                                    
###   OUTPUT:                                                          
###      - MSA files for each gene tree to the FASTA file format       
###                                                                    
###   Name: 02-MUSCLE.sh                  Author: Yoann Anselmetti     
###   Creation date: 2016/01/06           Last modification: 2017/03/29
###                                                                    

startTime=`date +%s`



if [ ! -d "$3" ]; then
	mkdir -p $3
fi

path=$(pwd)
GF=$(basename $2)


runtime=$((end-start))

echo -e "\tProcessing Gene Family $GF:"
# If 1 sequence is too long in $GF (>32000 char): use "-maxiters 2" option
seq_too_long=false
for line in $(cat $path/$2); do
   nb_char=$(wc -c <<< $line);
   if (($nb_char >= 32000)); then
      seq_too_long=true;
      break;
   fi
done
if $seq_too_long; then
	echo -e "\t\t$1 -in $path/$2 -out $path/$3/$GF -maxiters 2"
	$1 -in $path/$2 -out $path/$3/$GF -maxiters 2
else
	echo -e "\t\t$1 -in $path/$2 -out $path/$3/$GF"
	$1 -in $path/$2 -out $path/$3/$GF
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