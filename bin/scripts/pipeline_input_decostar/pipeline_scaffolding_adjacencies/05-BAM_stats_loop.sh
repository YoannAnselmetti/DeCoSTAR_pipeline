#!/bin/bash
###
###   Goal:
###       Run "BAM_stats.py" script for all libraries ([ES]RX) of 1 species (Must be executed for all species)
###
###   INPUT:
###      1- Directory containing BAM files
###         (/share/nas-isem_i/yanselmetti/ANOPHELES_PROJECT/DATA/DATA_SEQ/MAPPING/Bowtie2_ALL/TRIMMOMATIC3/ALL/Anopheles_albimanus)
###      2- TAG to allow severa experiments without interfering with jobs of others experiments
###         (Aalb_Bowtie2_ALL_RAW_ALL)
###      3- Multiple alignment used for MAPPING (If not given then the 50 "best" alignments are considered)
###         (a)
###
###   OUTPUT:
###      - stats file on contigs pairs links from "BAM_stats.py" script
###
###   Name: 05-BAM_stats_loop.sh              Author: Yoann Anselmetti
###   Creation date: 2015/10/21               Last modification: 2018/06/01
###

TAG=$2

job_max=200
k=50          # Nb of multiple alignments consider by Bowtie2
boolK=true    # IF false: consider all alignments with Bowtie2 if true: consider $k best alignments for Bowtie2

re='^[0-9]+$'
if [ -z $3 ]; then
   echo "Default usage"
   echo $k" multiple alignments were considered for mapping of reads on reference genome"
else
   if [ $3 = "a" ]; then
      echo "All multiple alignments were considered for mapping of reads on reference genome"
      boolK=false
   elif [[ $3 =~ $re ]] ; then
      echo $k" multiple alignments were considered for mapping of reads on reference genome"
      boolK=true
      k=$3
   else
      echo "ERROR, parameter 3 should be an integer or equal to \"a\" (if all alignments were considered during MAPPING step)!!!"
      exit
   fi
fi


script_BAM_stats="/share/nas-isem_i/yanselmetti/BIRDS_PROJECT/script/pipeline_input-decostar/DATA_SEQ_DECOSTAR/BAM_stats.py"

headSGE=#\!"/bin/bash\n#$ -V\n#$ -S /bin/bash\n#$ -cwd\n#$ -l h_rt=200000:00:00\n#$ -N "


function BAM_REC {
   DirList=$(find $1/* -prune)
   f1=""
   f2=""
   for DirFil in $DirList; do
      if [ -d "$DirFil" ]; then
#         echo "DIR $DirFil"
         # If directory is SRX directory
         if [ $(echo $DirFil | grep "[ES]RX") ]; then
            DIR=$(echo $DirFil | sed "s/MAPPING/STATS\/MAPPING/g")
            SRX=$(echo $DirFil | awk 'BEGIN{FS="/"};{print $NF}')

            # WAIT until a slot is available
            jobs_nb=$(qstat | grep $USER | wc -l)
            while [ ${jobs_nb} -ge ${job_max} ]; do
               echo "No spot yet"
               sleep 10
               jobs_nb=$(qstat | grep $USER | wc -l)
            done

            commandLine=""
            jobName=""
            if $boolK; then
               commandLine="module load python2.7;\n$script_BAM_stats $DirFil/"$SRX"_k"$k".sorted.bam $DIR"
            else
               commandLine="module load python2.7;\n$script_BAM_stats $DirFil/"$SRX".sorted.bam $DIR"
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
         else
            BAM_REC $DirFil
         fi
      fi
   done
}

cd $PBS_O_WORKDIR
echo "Current working directory is `pwd`"
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

echo "BAM_REC $1"
BAM_REC $1

echo "Finishing run at: `date`"
