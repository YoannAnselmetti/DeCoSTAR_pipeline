#!/bin/bash
###
###   Goal:
###      Script to scaffold contigs with paired reads (PE, MP or Fosill) on reference genome
###      (1 run/species => Submitting 1 job/SRX)
###        
###   INPUT:
###      1- PATH of INPUT MAPPING directory
###         (/share/nas-isem_i/yanselmetti/BIRDS_PROJECT/DATA_SEQ/MAPPING/Bowtie2_k100/TRIM/Acanthisitta_chloris)
###      2- OUTPUT directory name (BESST-2.2.6)
###      3- FASTA file of reference genome
###         (/share/nas-isem_i/yanselmetti/BIRDS_PROJECT/INPUT_DATA/FASTA/SCAFF/Acanthisitta_chloris-GCA_000695815.1_ASM69581v1_genomic.fna)
###      4- TAG to allow severa experiments without interfering with jobs
###         of others experiments (Achl_Bowtie2_k100_TRIM)
###      5- Orientation file
###         (/share/nas-isem_i/yanselmetti/BIRDS_PROJECT/DATA_SEQ/orientation_libraries)
###      6- Multiple alignment used for MAPPING (If not given then all alignments were considered)
###         (50)
###
###   Name: 06-run_BESST_oneRUNperSPECIES.sh   Author: Yoann Anselmetti
###   Creation date: 2015/08/26                Last modification: 2017/11/13
###

mapDIR=$1
suffDIR=$2
REF=$3
TAG=$4
file_orientation=$5


job_max=200
k=50          # Nb of multiple alignments consider by Bowtie2
boolK=true    # IF false: consider all alignments with Bowtie2 if true: consider $k best alignments for Bowtie2 


re='^[0-9]+$'
if [ -z $6 ]; then
   echo "Default usage"
   echo $k" multiple alignments were considered for mapping of reads on reference genome"
else
   if [ $6 = "a" ]; then
      echo "\tAll multiple alignments were considered for mapping of reads on reference genome"
      boolK=false
   elif [[ $6 =~ $re ]] ; then
      k=$6
      echo $k" multiple alignments were considered for mapping of reads on reference genome"
      boolK=true
   else
      echo "ERROR, parameter 7 should be an integer or equal to \"a\" (if all alignments were considered during MAPPING step)!!!"
      exit
   fi
fi

besst="/home/yanselmetti/Software/SCAFFOLDING/BESST-2.2.6/runBESST"
headSGE=#\!"/bin/bash\n#$ -V\n#$ -S /bin/bash\n#$ -cwd\n#$ -l h_rt=200000:00:00\n#$ -N "


cd $PBS_O_WORKDIR
echo "Current working directory is `pwd`"
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"


declare -A SRX_ori
while read line; do
   SRX=$(echo $line | cut -d" " -f2)
   # echo $SRX
   ori=$(echo $line | cut -d" " -f3)
   # echo $ori
   insertMean=$(echo $line | cut -d" " -f4)
   # echo $insertMean
   insertSD=$(echo $line | cut -d" " -f5)
   # echo $insertSD
   SRX_ori+=([$SRX]=$ori":"$insertMean":"$insertSD)
done < ${file_orientation}


# Create OUTPUT SCAFFOLDING Directory
scaffDIR=$(echo $mapDIR | sed "s/MAPPING/SCAFFOLDING\/${suffDIR}/g")
echo -e "mkdir -p $scaffDIR"
mkdir -p $scaffDIR


strSRX=""
strORI=""
strMean=""
strSD=""
listSRX=$(ls $mapDIR)
for SRX in $listSRX; do
   ori=$(echo ${SRX_ori["$SRX"]} | cut -d: -f1)
   insertMean=$(echo ${SRX_ori["$SRX"]} | cut -d: -f2)
   insertSD=$(echo ${SRX_ori["$SRX"]} | cut -d: -f3)
   echo $insertMean" "$insertSD":"$SRX":"$ori
done |
sort -n -k1 |( # Parenthesis very important due to the fact that var strSRX and strORI will be deleted once loop While is over...
while read -r elem; do
   # echo $elem
   mean=$(echo $elem | cut -d" " -f1)
   sd=$(echo $elem | cut -d" " -f2 | cut -d: -f1)
   SRX=$(echo $elem | cut -d" " -f2 | cut -d: -f2)
   ori=$(echo $elem | cut -d" " -f2 | cut -d: -f3)
   if $boolK; then
      strSRX=$strSRX"$mapDIR/$SRX/${SRX}_k$k.sorted.bam "
   else
      strSRX=$strSRX"$mapDIR/$SRX/$SRX.sorted.bam "
   fi
   strORI=$strORI$ori" "
   strSD=$strSD$sd" "
   strMean=$strMean$mean" "
done
# echo $strSRX
# echo $strORI

commandLine="module load python2.7; module load libz-1.2.11; module load atlas; $besst -c $REF -f $strSRX -orientation $strORI -m $strMean -s $strSD -o $scaffDIR -plot --print_scores -z 10000 --separate_repeats --min_mapq 0"
jobName="BESST_$TAG"


# WAIT until a slot is available
jobs_nb=$(qstat | grep $USER | wc -l)
while [ $jobs_nb -ge ${job_max} ]; do
   echo "No spot yet"
   sleep 10
   jobs_nb=$(qstat | grep $USER | wc -l)
done


# Print current job to run
echo -e "\t\t"$commandLine"\n"
# Run current job
echo -e $headSGE$jobName"\n\n"$commandLine > "job_"$jobName".sge"

sleep 10

jobID=$(qsub -q mem.q "job_"$jobName".sge" | cut -d' ' -f 3)
echo "Your job $jobID (\"$jobName\") has been submitted"
jobWaiting=true
while $jobWaiting; do
   sleep 1
   if qstat | grep -P "^$jobID" | grep -q " r   " || qstat | grep -P "^$jobID" | grep -q " qw   "; then
      jobWaiting=false
   else
      echo "Job $jobID is still not taken into account"
   fi
done
# rm "job_"$jobName".sge"
)

echo "Finishing run at: `date`"
