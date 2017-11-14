#!/bin/bash
###
###   Goal:
###      Script to trim reads to remove reads of low qualities on cluster
###
###   INPUT:
###      1- TRIMMING tool (ERNE-filter or TRIMMOMATIC)
###      2- PATH of INPUT directory containing FASTQ files
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/FASTQ/RAW/ALL)
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/FASTQ/RAW/50_pourc)
###      3- TAG to allow several experiments without interfering with jobs of others experiments
###         (ALL / 50pourc)
###      4- Trimmomatic paramters option 
###         (1 / 2 / 3 / 4)
###
###   Name: 01-run_TRIMMING.sh                 Author: Yoann Anselmetti
###   Creation date: 2015/09/21                Last modification: 2017/07/24
###

TRIMMING_tool=$1
FASTQ_dir=$2
TAG=$3
trimmomaticOpt=$4

job_max=200


erne_filter="/home/yanselmetti/Software/TRIMMING_READS/erne-2.1.1-linux/bin/erne-filter"
trimmomatic="/home/yanselmetti/Software/TRIMMING_READS/Trimmomatic-0.36/trimmomatic-0.36.jar"

headSGE=#\!"/bin/bash\n#$ -V\n#$ -S /bin/bash\n#$ -cwd\n#$ -l h_rt=200000:00:00\n#$ -N "


function FASTQ_REC {
   DirList=$(ls $1)
   f1=""
   f2=""
   for DirFil in $DirList; do
      if $(echo "$DirFil" | grep -q "[ES]RR"); then
         SRR=$DirFil
         f1=$1/$SRR/$SRR"_1.fastq.gz"
         f2=$1/$SRR/$SRR"_2.fastq.gz"
         echo -e "\tfastq_1 = "$f1
         echo -e "\tfastq_2 = "$f2
         DIR=""
         if [[ ${TRIMMING_tool} == "ERNE-filter" ]]; then
            DIR=$(echo $1/$SRR | sed "s/FASTQ\/RAW/FASTQ\/${TRIMMING_tool}/g")
         elif [[ ${TRIMMING_tool} == "TRIMMOMATIC" ]]; then
            DIR=$(echo $1/$SRR | sed "s/FASTQ\/RAW/FASTQ\/${TRIMMING_tool}${trimmomaticOpt}/g")
         else
            echo "ERROR, trimming tool should be \"TRIMMOMATIC\" or \"ERNE-FILTER\" !!!"
            exit 1
         fi
         echo -e "\tmkdir -p $DIR"
         mkdir -p $DIR

         # WAIT until a slot is available
         jobs_nb=$(qstat | grep $USER | wc -l)
         while [ ${jobs_nb} -ge ${job_max} ]; do
            echo "No spot yet"
            sleep 10
            jobs_nb=$(qstat | grep $USER | wc -l)
         done

         commandLine=""
         jobName=""
         if [[ ${TRIMMING_tool} == "ERNE-filter" ]]; then
            #######################
            ### ERNE-filter run ###
            #######################*
            commandLine="${erne_filter} --threads 4 --query1 $f1 --query2 $f2 --output-prefix $DIR/$SRR;\nrm $DIR/${SRR}_unpaired.fastq;\ngzip $DIR/$SRR*"
            jobName="erne-filter_${SRR}_$TAG"

         elif [[ ${TRIMMING_tool} == "TRIMMOMATIC" ]]; then
            #######################
            ### TRIMMOMATIC run ###
            #######################
            if [[ ${trimmomaticOpt} == "1" ]]; then
               # OPTION 1: BASIC FILTER parameters for TRIMMOMATIC
               commandLine="java -jar $trimmomatic PE $f1 $f2 $DIR/${SRR}_1.fastq.gz $DIR/${SRR}_1_unpaired $DIR/${SRR}_2.fastq.gz $DIR/${SRR}_2_unpaired LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36;\nrm $DIR/${SRR}_1_unpaired $DIR/${SRR}_2_unpaired;"
               jobName="trimmomatic_opt${trimmomaticOpt}_${SRR}_$TAG"

            elif [[ ${trimmomaticOpt} == "2" ]]; then
               # OPTION 2: CROP extremities of reads for TRIMMOMATIC => BEST option (See FASTQC report for Anopheles_albimanus/SRX200219/SRR606148)
               commandLine="java -jar $trimmomatic PE $f1 $f2 $DIR/${SRR}_1.fastq.gz $DIR/${SRR}_1_unpaired $DIR/${SRR}_2.fastq.gz $DIR/${SRR}_2_unpaired CROP:85 HEADCROP:10;\nrm $DIR/${SRR}_1_unpaired $DIR/${SRR}_2_unpaired;"
               jobName="trimmomatic_opt${trimmomaticOpt}_${SRR}_$TAG"

            elif [[ ${trimmomaticOpt} == "3" ]]; then
               # OPTION 3: CROP extremities + BASIC FILTER for TRIMMOMATIC
               commandLine="java -jar $trimmomatic PE $f1 $f2 $DIR/${SRR}_1.fastq.gz $DIR/${SRR}_1_unpaired $DIR/${SRR}_2.fastq.gz $DIR/${SRR}_2_unpaired CROP:85 HEADCROP:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75;\nrm $DIR/${SRR}_1_unpaired $DIR/${SRR}_2_unpaired;"
               jobName="trimmomatic_opt${trimmomaticOpt}_${SRR}_$TAG"
            elif [[ ${trimmomaticOpt} == "3" ]]; then
               # OPTION 3: CROP extremities + BASIC FILTER for TRIMMOMATIC
               commandLine="java -jar $trimmomatic PE $f1 $f2 $DIR/${SRR}_1.fastq.gz $DIR/${SRR}_1_unpaired $DIR/${SRR}_2.fastq.gz $DIR/${SRR}_2_unpaired CROP:85 HEADCROP:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75;\nrm $DIR/${SRR}_1_unpaired $DIR/${SRR}_2_unpaired;"
               jobName="trimmomatic_opt${trimmomaticOpt}_${SRR}_$TAG"
            else
               echo "ERROR, option parameter for 4 should be \"1\", \"2\", \"3\" or \"4\" and not $trimmomaticOpt"
               exit 1
            fi

         else
            echo "ERROR, trimming tool should be \"TRIMMOMATIC\" or \"ERNE-FILTER\" !!!"
            exit 1
         fi

         # Print current job to run
         echo -e "\t\t"$commandLine"\n"
         # Run current job
         echo -e $headSGE$jobName"\n\n"$commandLine > "job_"$jobName".sge"
         jobID=$(qsub -q mem.q "job_"$jobName".sge" | cut -d' ' -f 3)
         echo "Your job $jobID (\"$jobName\") has been submitted"
         jobWaiting=true
         while $jobWaiting; do
            sleep 1
            if qstat | grep -P "^$jobID" | grep -q " r   " || qstat | grep -P "^$jobID" | grep -q " qw   "; then
               jobWaiting=false
            fi
         done
         # rm "job_"$jobName".sge"

      elif [ -d "$1/$DirFil" ]; then
         echo "DIR $1/$DirFil"
         FASTQ_REC $1/$DirFil
      fi
   done
}


echo "Current working directory is `pwd`"
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

FASTQ_REC ${FASTQ_dir}

echo "Finishing run at: `date`"