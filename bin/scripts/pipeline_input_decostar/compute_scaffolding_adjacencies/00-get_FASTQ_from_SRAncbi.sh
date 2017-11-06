#!/bin/bash
###
###   Prerequisites: Need to create the directory architecture of the SRX*/SRR* ID for each species
###   ex: mkdir -p Anopheles_albimanus/SRX111456/SRR389778 Anopheles_albimanus/SRX111456/SRR389781 Anopheles_albimanus/SRX111456/SRR390324 Anopheles_albimanus/SRX111456/SRR390326 Anopheles_albimanus/SRX200219/SRR606148
###
###   Goal:
###      Script to obtain FASTQ files from SRA files from NCBI website
###
###   Comments:
###      Get FASTQ directly from NCBI website without downloading SRA
###
###   INPUT:
###      1- INPUT directory containing SRA files
###      2- TAG to allow several experiments without interfering with jobs of others experiments (Aalb)
###
###   OUTPUT:
###      - OUTPUT directory containing directory for each species that will contain FASTQ files
###
###   Name: get_FASTQ_from_SRAncbi.sh        Author: Yoann Anselmetti
###   Creation date: 2015/07/30              Last modification: 2016/12/09
###

SRAdir=$1
TAG=$2

job_max=100

fastqdump='/home/yanselmetti/Software/SRA/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump.2.8.0'

headSGE=#\!"/bin/bash\n#$ -V\n#$ -S /bin/bash\n#$ -cwd\n#$ -l h_rt=200000:00:00,mem_free=4G\n#$ -N "


function DIR_REC {
	DirList=$(ls $1)
	for DirFil in $DirList; do
		DIR=$1/$DirFil
		if $(echo "$DirFil" | grep -q "SRR"); then
			SRR=$DirFil
			echo $DIR"/"

         # WAIT until a slot is available
         jobs_nb=$(qstat | grep $USER | wc -l)
         while [ ${jobs_nb} -ge ${job_max} ]; do
            echo "No spot yet"
            sleep 10
            jobs_nb=$(qstat | grep $USER | wc -l)
         done

         # Command line to run
         commandLine="$fastqdump --split-3 --gzip -O $DIR $SRR;"
         jobName="SRA_${SRR}_$TAG"

			# Print current job to run
			echo -e "\t\t"$commandLine
			# Run current job
			echo -e $headSGE$jobName"\n\n"$commandLine > "job_"$jobName".sge"
			jobID=$(qsub "job_"$jobName".sge" | cut -d' ' -f 3)
			echo -e "Your job $jobID (\"$jobName\") has been submitted\n"

			# Slow job execution to avoid file truncation!!!
			jobWaiting=true
			while $jobWaiting; do
				sleep 1
				if qstat | grep -P "^$jobID" | grep -q " r   " || qstat | grep -P "^$jobID" | grep -q " qw   "; then
					jobWaiting=false
				fi
			done
         # rm "job_"$jobName".sge"

		elif [ -d "$DIR" ]; then
			# echo $DIR"/"
			DIR_REC $DIR
		fi
	done
}


DIR_REC $SRAdir