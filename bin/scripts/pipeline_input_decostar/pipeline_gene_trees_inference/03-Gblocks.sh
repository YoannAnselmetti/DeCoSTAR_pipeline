#!/bin/bash
###                                                                    
###   Goal:                                                            
###      Alignment of CDS sequence of gene from same gene trees        
###                                                                    
###   INPUT:                                                           
###      1- PATH of "Gblocks" executive file                           
###         (Software/TRIMMING_ALIGN/Gblocks_0.91b/Gblocks)            
###      2- INPUT MSA tool (MUSCLE) file                               
###         (DECOSTAR/DATA_DeCoSTAR/FASTA/MSA/CDS/MUSCLE/MZ22500058)   
###      3- OUTPUT directory that will contain output of Gblocks       
###         (DECOSTAR/DATA_DeCoSTAR/FASTA/MSA/CDS/Gblocks/)            
###      4- PROTEIN or CDS boolean (Protein: P/p | CDS: C/c)           
###                                                                    
###   OUTPUT:                                                          
###      - MSA files filter with Gblocks to the FASTA file format      
###                                                                    
###   Name: 03-Gblocks.sh                 Author: Yoann Anselmetti     
###   Creation date: 2016/01/06           Last modification: 2017/03/29
###                                                                    

startTime=`date +%s`



if [ ! -d "$3" ]; then
	mkdir -p $3
fi


MSA=$(basename $2)
MUSCLE_dir=$(dirname $2)

echo -e "Processing Gene Family $MSA:"
# If PROT sequence
if [[ $4 == "P" || $4 == "p" ]]; then
   echo -e "\t$1 $2 -t=p -b4=5 -b5=h"
   $1 $2 -t=p -b4=5 -b5=h
# If RNA sequence
elif [[ $4 == "C" || $4 == "c" ]]; then
   echo -e "\t$1 $2 -t=d -b4=5 -b5=h"
   $1 $2 -t=d -b4=5 -b5=h
else
   echo -e "\n!!! ERROR: Param 4 should be p (for protein) or r (for RNA) and not $4 !!!"
   exit
fi

# Move OUTPUT file in OUTPUT directory and remove html file produced by Gblocks
mv $MUSCLE_dir/$MSA"-gb" $3
rm $MUSCLE_dir/$MSA"-gb.htm"



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
