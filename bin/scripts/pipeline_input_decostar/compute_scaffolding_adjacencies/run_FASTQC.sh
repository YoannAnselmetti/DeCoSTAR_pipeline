#!/bin/bash
###
###   Goal:
###      Script to get statistics on FASTQ files with FASTQc programm
###
###   INPUT:
###      1- PATH of INPUT directory containing FASTQ files
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/FASTQ/RAW/ALL/Anopheles_albimanus)
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/FASTQ/TRIMMOMATIC3/ALL/Anopheles_albimanus)
###      2- TAG to allow severa experiments without interfering with jobs of others experiments
###         (Aalb_RAW_ALL)
###         (Aalb_TRIMMOMATIC3_ALL)
###
###   Name: run_FASTQC.sh                      Author: Yoann Anselmetti
###   Creation date: 2016/10/18                Last modification: 2016/12/09
###

FASTQ_dir=$1
TAG=$2

job_max=100

fastqc="/home/yanselmetti/Software/TRIMMING_READS/FastQC-v0.11.5/fastqc"

headSGE=#\!"/bin/bash\n#$ -V\n#$ -S /bin/bash\n#$ -cwd\n#$ -l h_rt=200000:00:00\n#$ -N "

function FASTQc_REC {
   DirList=$(ls $1)
   f1=""
   f2=""
   for DirFil in $DirList; do
      if $(echo "$DirFil" | grep -q "SRR"); then
         SRR=$DirFil
         f1=$1/$SRR/$SRR"_1.fastq.gz"
         f2=$1/$SRR/$SRR"_2.fastq.gz"
         echo -e "\tfastq_1 = "$f1
         echo -e "\tfastq_2 = "$f2

         # WAIT until a slot is available
         jobs_nb=$(qstat | grep $USER | wc -l)
         while [ ${jobs_nb} -ge ${job_max} ]; do
            # echo "No spot yet"
            sleep 10
            jobs_nb=$(qstat | grep $USER | wc -l)
         done

         commandLine=$fastqc" "$f1" "$f2
         jobName="FASTQc_${SRR}_$TAG"

         # Print current job to run
         echo -e "\t\t"$commandLine"\n"
         # Run current job
         echo -e $headSGE$jobName"\n\n"$commandLine > "job_"$jobName".sge"
         jobID=$(qsub -q long.q "job_"$jobName".sge" | cut -d' ' -f 3)
         echo "Your job $jobID (\"$jobName\") has been submitted"
         jobWaiting=true
         while $jobWaiting; do
            sleep 1
            if qstat | grep -P "^$jobID" | grep -q " r   " || qstat | grep -P "^$jobID" | grep -q " qw   "; then
               jobWaiting=false
            # else
            #    echo "Job $jobID is still not taken into account"
            fi
         done
         rm "job_"$jobName".sge"

      elif [ -d "$1/$DirFil" ]; then
         echo "DIR $1/$DirFil"
         FASTQc_REC $1/$DirFil
      fi
   done
}


echo "Current working directory is `pwd`"
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

FASTQc_REC ${FASTQ_dir}

echo "Finishing run at: `date`"