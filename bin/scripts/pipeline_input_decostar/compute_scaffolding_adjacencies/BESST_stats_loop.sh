#!/bin/bash
###
###   Goal:
###       Run "BESST_stats.py" script for all species for a choosen library orientation
###
###   INPUT:
###      1- Orientation library choosen (fr/rf/all)
###      2- Directory containing BAM files
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_ALL/RAW/ALL/Anopheles_albimanus)
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/SCAFFOLDING/BESST-2.2.5/Bowtie2_k100/RAW/ALL/Anopheles_albimanus)
###      3- TAG to allow severa experiments without interfering with jobs of others experiments
###         (Aalb_Bowtie2_ALL_RAW_ALL)
###         (Aalb_Bowtie2_k100_RAW_ALL)
###
###   OUTPUT:
###      - stats file on contigs pairs links and gap size from "BESST_stats.py" script
###
###   Name: BESST_stats_loop.sh              Author: Yoann Anselmetti
###   Creation date: 2015/10/21              Last modification: 2017/01/17
###

TAG=$3

script_BESST_stats="/home/yanselmetti/Anopheles_project/script/cluster/cluster_ISEM/DATA_SEQ_DeCoSTAR/BESST_stats.py"
file_orientation="/share/nas-isem_i/yanselmetti/DATA_SEQ/orientation_libraries"

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
done < $file_orientation


#for SRX in "${!SRX_ori[@]}"; do
#   echo "$SRX ${SRX_ori["$SRX"]}"
#done

function SCAFF_REC {
  DirList=$(find $1/* -prune)
  for DirFil in $DirList; do
    if [ -d "$DirFil" ]; then
      #         echo "DIR $DirFil"
      # If directory is SRX directory
      if [ $(echo $DirFil | grep "SRX") ]; then
        echo $DirFil
        DIR=$(echo $DirFil | sed "s/SCAFFOLDING/STATS\/SCAFFOLDING/g")
        SRX=$(echo $DirFil | awk 'BEGIN{FS="/"};{print $NF}')
        orientation=""
        orientation=${SRX_ori["$SRX"]}
        if [[ $orientation != "" ]]; then
          mkdir -p $DIR
          # echo $orientation
          commandLine="module load python2.7;\n$script_BESST_stats $DirFil/BESST_output/score_file_pass_1.tsv $DIR"
          jobName="SCAFFstats_${SRX}_$TAG"
          # Print current job to run
          echo -e "\t\t"$commandLine"\n"
          # Run current job
          echo -e $headSGE$jobName"\n\n"$commandLine > "job_"$jobName".sge"
          jobID=$(qsub "job_"$jobName".sge" | cut -d' ' -f 3)
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
        SCAFF_REC $DirFil
      fi
    fi
  done
}

cd $PBS_O_WORKDIR
echo "Current working directory is `pwd`"
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

SCAFF_REC $2

echo "Finishing run at: `date`"
