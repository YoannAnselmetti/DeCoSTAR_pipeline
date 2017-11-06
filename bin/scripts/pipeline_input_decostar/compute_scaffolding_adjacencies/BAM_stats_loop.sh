#!/bin/bash
###
###   Goal:
###       Run "BAM_stats.py" script for all species for a choosen library orientation
###
###   INPUT:
###      1- Orientation library choosen (fr/rf/all)
###      2- Directory containing BAM files
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/MAPPING/Bowtie2_ALL/RAW/ALL/Anopheles_albimanus)
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/MAPPING/Bowtie2_k100/RAW/ALL/Anopheles_albimanus)
###      3- TAG to allow severa experiments without interfering with jobs of others experiments
###         (Aalb_Bowtie2_ALL_RAW_ALL)
###         (Aalb_Bowtie2_k100_RAW_ALL)
###      4- Orientation file
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/orientation_libraries)
###      5- Multiple alignment used for MAPPING (If not given then all alignments were considered)
###         (50)
###
###   OUTPUT:
###      - stats file on contigs pairs links from "BAM_stats.py" script
###
###   Name: BAM_stats_loop.sh                 Author: Yoann Anselmetti
###   Creation date: 2015/10/21               Last modification: 2017/03/15
###

TAG=$3

job_max=130
k=50          # Nb of multiple alignments consider by Bowtie2
boolK=true    # IF false: consider all alignments with Bowtie2 if true: consider $k best alignments for Bowtie2

re='^[0-9]+$'
if [ -z $5 ]; then
   echo "Default usage"
   echo $k" multiple alignments were considered for mapping of reads on reference genome"
else
   if [ $5 = "a" ]; then
      echo "All multiple alignments were considered for mapping of reads on reference genome"
      boolK=false
   elif [[ $5 =~ $re ]] ; then
      echo $k" multiple alignments were considered for mapping of reads on reference genome"
      boolK=true
      k=$5
   else
      echo "ERROR, parameter 5 should be an integer or equal to \"a\" (if all alignments were considered during MAPPING step)!!!"
      exit
   fi
fi


script_BAM_stats="/home/yanselmetti/Anopheles_project/script/cluster/cluster_ISEM/DATA_SEQ_DeCoSTAR/BAM_stats.py"


headSGE=#\!"/bin/bash\n#$ -V\n#$ -S /bin/bash\n#$ -cwd\n#$ -l h_rt=200000:00:00\n#$ -N "


declare -A SRX_ori

while read line; do
   SRX=$(echo $line | cut -d" " -f2)
#   echo $SRX
   ori=$(echo $line | cut -d" " -f3)
#   echo $ori
   if [[ $1 == "fr" ]]; then
      if [[ $ori == "fr" ]]; then
         SRX_ori+=([$SRX]=$ori)
      fi
   elif [[ $1 == "rf" ]]; then
      if [[ $ori == "rf" ]]; then
         SRX_ori+=([$SRX]=$ori)
      fi
   elif [[ $1 == "all" ]]; then
      SRX_ori+=([$SRX]=$ori)
   else
      echo "Error the 2nd paramter is incorrect: Should be \"fr\" for Paired-End or fosill library OR \"rf\" for Mate-Pairs library OR \"all\" if all orientations!!!"
      exit 1
   fi
done < $4

#for SRX in "${!SRX_ori[@]}"; do
#   echo "$SRX ${SRX_ori["$SRX"]}"
#done

function BAM_REC {
   DirList=$(find $1/* -prune)
   f1=""
   f2=""
   for DirFil in $DirList; do
      if [ -d "$DirFil" ]; then
#         echo "DIR $DirFil"
         # If directory is SRX directory
         if [ $(echo $DirFil | grep "SRX") ]; then
            DIR=$(echo $DirFil | sed "s/MAPPING/STATS\/MAPPING/g")
            SRX=$(echo $DirFil | awk 'BEGIN{FS="/"};{print $NF}')
            orientation=${SRX_ori["$SRX"]}

            # WAIT until a slot is available
            jobs_nb=$(qstat | grep $USER | wc -l)
            while [ ${jobs_nb} -ge ${job_max} ]; do
               echo "No spot yet"
               sleep 10
               jobs_nb=$(qstat | grep $USER | wc -l)
            done

            commandLine=""
            jobName=""
            if [[ $orientation != "" ]]; then
               # echo $orientation
               if $boolK; then
                  commandLine="module load python2.7;\n$script_BAM_stats $DirFil/"$SRX"_k"$k".sorted.bam $DIR $orientation"
               else
                  commandLine="module load python2.7;\n$script_BAM_stats $DirFil/"$SRX".sorted.bam $DIR $orientation"
               fi
               jobName="BAMstats_${SRX}_$TAG"
               # Print current job to run
               echo -e "\t\t"$commandLine"\n"
               # Run current job
               echo -e $headSGE$jobName"\n\n"$commandLine > "job_"$jobName".sge"
               jobID=$(qsub -V "job_"$jobName".sge" | cut -d' ' -f 3)
               echo "Your job $jobID (\"$jobName\") has been submitted"
               jobWaiting=true
               while $jobWaiting; do
                  sleep 1
                  if qstat | grep -P "^$jobID" | grep -q " r   " || qstat | grep -P "^$jobID" | grep -q " qw   "; then
                     jobWaiting=false
                  fi
               done
               # rm "job_"$jobName".sge"
            fi
         else
            BAM_REC $DirFil $2
         fi
      fi
   done
}

cd $PBS_O_WORKDIR
echo "Current working directory is `pwd`"
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

echo "BAM_REC $2 $1"
BAM_REC $2 $1

echo "Finishing run at: `date`"
