#!/bin/bash
###
###   Goal:
###      Script to index ref FASTA file with "bowtie2-build" or "bwa index" before mapping with bowtie2 or bwa
###
###   INPUT:
###      1- MAPPING tool (mrsFAST, Bowtie2 or BWA)
###      2- Directory containing reference genomes FASTA files
###         (/share/nas-isem_i/yanselmetti/INPUT_DATA/FASTA/SCAFF)
###      3- INDEX of reference genome
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/INDEX/...)
###      4- TAG to allow severa experiments without interfering with jobs of others experiments
###         (blastn_50pourc)
###
###   Name: 02-index_MAPPING.sh                 Author: Yoann Anselmetti
###   Creation date: 2015/08/24                 Last modification: 2017/07/19
###


MAPPING_tool=$1
FASTA_dir=$2
INDEX_dir=$3
TAG=$4

job_max=200

mrsFAST="/home/yanselmetti/Software/MAPPING/mrsfast-3.3.9/mrsfast"
bowtie2_build="/home/yanselmetti/Software/MAPPING/bowtie2-2.2.9/bowtie2-build" # IF CHANGE NAME OF READ, THIS ONE CAN BE MAPPED TO AN OTHER POSITION IF MULTIPLE ALIGNMENT (Firtina & Alkan, 2016)
bwa="/home/yanselmetti/Software/MAPPING/bwa-0.7.15/bwa" # NOT RECOMMANDED SINCE BWA HAS BEEN SHOWN TO BE NON-DETERMINISTIC METHODS (Firtina & Alkan, 2016)

headSGE=#\!"/bin/bash\n#$ -V\n#$ -S /bin/bash\n#$ -cwd\n#$ -l h_rt=200000:00:00,mem_free=4G\n#$ -N "


function FASTA_REC {
   DirList=$(find $1/* -prune)
   for DirFil in $DirList; do
      echo $DirFil
      if [ -f "$DirFil" ]; then
#         basename $DirFil
################################################################################
### COMMAND TO RECOVER SPECIES NAME => MAKE IT USABLE FOR SEVERAL CASES !!!! ###
################################################################################
         species=$(basename $DirFil | cut -d. -f1)
         # echo $species
         mkdir -p $2/$species

         # WAIT until a slot is available
         jobs_nb=$(qstat | grep $USER | wc -l)
         while [ ${jobs_nb} -ge ${job_max} ]; do
            echo "No spot yet"
            sleep 10
            jobs_nb=$(qstat | grep $USER | wc -l)
         done

         if [[ ${MAPPING_tool} == *"BWA"* ]]; then
            commandLine="$bwa index -p $2/$species/$species $DirFil"
            jobName="Index_BWA_$3_$species"

         elif [[ ${MAPPING_tool} == *"Bowtie2"* ]]; then
            commandLine="${bowtie2_build} $DirFil $2/$species/$species"
            jobName="Index_Bowtie2_$3_$species"

         elif [[ ${MAPPING_tool} == *"mrsFAST"* ]]; then
            commandLine="$mrsFAST --index $DirFil --ws 14;\nmv $DirFil".index" $2/$species/$species".index""
            jobName="Index_mrsFAST_$3_$species"
            
         else
            echo "ERROR: mapping tool should be \"BWA\", \"Bowtie2\" or \"mrsFAST\" !!!"
            exit 1
         fi

#######################
### RUN CURRENT JOB ###
#######################
         # Print current job to run
         echo -e "\t\t"$commandLine"\n"
         # Run current job
         echo -e $headSGE$jobName"\n\n"$commandLine > "job_"$jobName".sge"

         sleep 5

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

      elif [ -d "$DirFil" ]; then
         FASTA_REC $DirFil $2 $3
      fi
   done
}

cd $PBS_O_WORKDIR
echo "Current working directory is `pwd`"
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

FASTA_REC ${FASTA_dir} ${INDEX_dir} $TAG

echo "Finishing run at: `date`"
