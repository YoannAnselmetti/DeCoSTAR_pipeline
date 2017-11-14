#!/bin/bash
###
###   Goal:
###      Script to get statistics on FASTQ files with FASTQc programm
###
###   INPUT:
###      1- PATH of INPUT directory containing MAPPING STATS files
###         (/share/nas-isem_i/yanselmetti/BIRDS_PROJECT/DATA_SEQ/STATS/MAPPING/Bowtie2_k100)
###      2- script to write orientation file form MAPPING STATS file for one species
###         (/share/nas-isem_i/yanselmetti/BIRDS_PROJECT/script/pipeline_input-decostar/DATA_SEQ_DECOSTAR/write_orientation_file.py)
###      3- orientation file
###         (/share/nas-isem_i/yanselmetti/BIRDS_PROJECT/DATA_SEQ/orientation_libraries)
###
###   Name: write_orientation_file_for_ALL_species.sh   Author: Yoann Anselmetti
###   Creation date: 2017/11/13                         Last modification: 2017/11/13
###

STATS_dir=$1
script=$2
orientation_file=$3
echo "Current working directory is `pwd`"
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

for spe in $(ls ${STATS_dir}); do
    $script ${STATS_dir} $spe ${orientation_file};
done

echo "Finishing run at: `date`"
