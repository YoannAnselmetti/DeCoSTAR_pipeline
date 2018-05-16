# Cedric Chauve, December 18, 2016
#!/bin/sh

# $1 = reconciled gene trees: file name from the bin directory
# $2 = species tree: same
# $3 = DeCo* adjacencies file: same
# $4 = genes coordinates file: same
# Note: $1, $2, $3 are obtained from DeCo*
# $5 = threshold for linearization
# $6 = directory of results: same
# $7 = prefix of files
# $8 = method of linearization
# $9 = python command

# Setting up the python path to find the networkx module
PWD=`pwd`
PYTHONPATH=${PWD}/bin/scripts/post_decostar/code:$PYTHONPATH

SCF=${6}/${7}_${8}_scaffolds
OG=${6}/${7}_orthogroups
SCJ=${6}/${7}_${8}_scj

codeDIR=bin/scripts/post_decostar/code

startTime=`date +%s`

echo "--> Linearizing genomes with method" $8
$9 $codeDIR/linearize_genomes.py $3 ALL $5 ${6}/${7} ${8} 1.0 0.001 > ${SCF}_log

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





startTime=`date +%s`

echo "--> Assigning scaffolds to chromosomes"
$9 $codeDIR/assign_scaffolds_X.py ${SCF} $1 $2 $4

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





startTime=`date +%s`

echo "--> Creating orthogroups from reconciled gene trees"
$9 $codeDIR/list_ancestors_descendants.py $1 $2 ${OG}

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





startTime=`date +%s`

echo "--> Detecting SCJ rearrangements"
$9 $codeDIR/scj.py ${OG} ALL ${SCF} 3 ${SCJ} > ${SCJ}_log

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





startTime=`date +%s`

echo "--> Counting SCJ per chromosome"
ASSIGNMENT=${SCF}_assignment
$9 $codeDIR/count_X_scj.py $1 ${OG} ${SCF} ${ASSIGNMENT} ${SCJ} $2 > ${SCJ}_stats

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


# for MISSING in 1 2 3 4 5 6 7 8 9 10
# do
#     echo "--> Duplication-free gene families missing in at most " ${MISSING} " genomes"
#     OG1=${OG}_DF${MISSING}
#     SCF1=${SCF}_DF${MISSING}
#     SCJ1=${SCJ}_DF${MISSING}
#     ASSIGNMENT1=${SCF1}_assignment
#     echo "----> Creating orthogroups from reconciled gene trees"
#     $9 $codeDIR/list_ancestors_descendants.py $1 $2 ${OG1} ${MISSING}
#     echo "----> Filtering scaffolds"
#     $9 $codeDIR/filter_scaffolds_by_outgroups.py ${SCF} ${OG1} ${SCF1}
#     echo "----> Detecting SCJ rearrangements"
#     $9 $codeDIR/scj.py ${OG1} ALL ${SCF} 3 ${SCJ1} > ${SCJ1}_log
#     echo "----> Counting SCJ per chromosome"
#     $9 $codeDIR/count_X_scj.py $1 ${OG1} ${SCF1} ${ASSIGNMENT} ${SCJ1} $2 > ${SCJ1}_stats1
#     echo "----> Assigning scaffolds to filtered chromosomes and recounting scj"
#     $9 $codeDIR/assign_scaffolds_X.py ${SCF1} $1 $2 $4
#     $9 $codeDIR/count_X_scj.py $1 ${OG1} ${SCF1} ${ASSIGNMENT1} ${SCJ1} $2 > ${SCJ1}_stats2
# done
