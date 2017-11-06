#!/bin/bash
###
###   Goal:
###      Script to scaffold contigs with paired reads (PE, MP or Fosill) on reference genome
###      (1 run/species => Submitting 1 job/SRX)
###        
###   INPUT:
###      1- PATH of INPUT MAPPING directory
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/MAPPING/Bowtie2_ALL/TRIMMOMATIC3/ALL/Anopheles_albimanus)
###      2- SCAFFOLDING tool (BESST or OPERA-LG)
###      3- OUTPUT directory name (BESST-2.2.5)
###      4- FASTA file of reference genome
###         (/share/nas-isem_i/yanselmetti/INPUT_DATA/FASTA/SCAFF/Anopheles-albimanus-STECLA_SCAFFOLDS_AalbS1.fa)
###      5- TAG to allow severa experiments without interfering with jobs
###         of others experiments (Aalb_Bowtie2_ALL_TRIMMOMATIC3_ALL)
###      6- Orientation file
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/orientation_libraries)
###      7- Multiple alignment used for MAPPING (If not given then all alignments were considered)
###         (50)
###
###   Name: run_SCAFFOLDING.sh                 Author: Yoann Anselmetti
###   Creation date: 2015/08/26                Last modification: 2017/03/15
###

MAPPING_dir=$1
SCAFFOLDING_tool=$2
SCAFF_dir=$3
REF=$4
TAG=$5
file_orientation=$6

job_max=100
k=50          # Nb of multiple alignments consider by Bowtie2
boolK=true    # IF false: consider all alignments with Bowtie2 if true: consider $k best alignments for Bowtie2 

re='^[0-9]+$'
if [ -z $7 ]; then
   echo "Default usage"
   echo $k" multiple alignments were considered for mapping of reads on reference genome"
else
   if [ $7 = "a" ]; then
      echo "All multiple alignments were considered for mapping of reads on reference genome"
      boolK=false
   elif [[ $7 =~ $re ]] ; then
      echo $k" multiple alignments were considered for mapping of reads on reference genome"
      boolK=true
      k=$7
   else
      echo "ERROR, parameter 7 should be an integer or equal to \"a\" (if all alignments were considered during MAPPING step)!!!"
      exit
   fi
fi

besst="/home/yanselmetti/Software/SCAFFOLDING/BESST-2.2.6/runBESST"
operalg="/home/yanselmetti/Software/SCAFFOLDING/OPERA-LG_v2.0.5/bin/OPERA-LG"

headSGE=#\!"/bin/bash\n#$ -V\n#$ -S /bin/bash\n#$ -cwd\n#$ -l h_rt=200000:00:00\n#$ -N "


declare -A SRX_ori

while read line; do
   SRX=$(echo $line | cut -d" " -f2)
#   echo $SRX
   ori=$(echo $line | cut -d" " -f3)
#   echo $ori
   insert=$(echo $line | cut -d" " -f4)
#   echo $insert
   SRX_ori+=([$SRX]=$ori)
   # SRX_ori+=([$SRX]=$ori":"$insert)
done < ${file_orientation}


function BAM_REC {
   DirList=$(find $1/* -prune)
   for DirFil in $DirList; do
      if [ -d "$DirFil" ]; then
         echo "DIR $DirFil"
         # If directory is SRX directory
         if [ $(echo $DirFil | grep "SRX") ]; then
            SRX=$(echo $DirFil | awk 'BEGIN{FS="/"};{print $NF}')
            # echo -e "$SRX"
            DIR=$(echo $DirFil | sed "s/MAPPING/SCAFFOLDING\/${SCAFF_dir}/g")
            echo -e "\tmkdir -p $DIR"
            mkdir -p $DIR
            orientation=${SRX_ori["$SRX"]}
            # orientation=$(echo ${SRX_ori["$SRX"]} | cut -d: -f1)
            # insert=$(echo ${SRX_ori["$SRX"]} | cut -d: -f2)

            # WAIT until a slot is available
            jobs_nb=$(qstat | grep $USER | wc -l)
            while [ $jobs_nb -ge ${job_max} ]; do
               echo "No spot yet"
               sleep 10
               jobs_nb=$(qstat | grep $USER | wc -l)
            done

            commandLine=""
            jobName=""
            if [[ ${SCAFFOLDING_tool} == "BESST" ]]; then
               if $boolK; then
                  commandLine="module load python2.7; module load libz-1.2.11; module load atlas; $besst -c $2 -f $DirFil/${SRX}_k$k.sorted.bam -plot --print_scores -o $DIR -orientation $orientation -z 10000 --separate_repeats --min_mapq 0"
               else
                  commandLine="module load python2.7; module load libz-1.2.11; module load atlas; $besst -c $2 -f $DirFil/$SRX.sorted.bam -plot --print_scores -o $DIR -orientation $orientation -z 10000 --separate_repeats --min_mapq 0"
               fi
               jobName="BESST_${SRX}_$TAG"

            elif [[ ${SCAFFOLDING_tool} == "OPERA-LG" ]]; then
               ## TO DO !!! => Need to create config file for OPERA-LG and mapping data !!!
               commandLine="$operalg $2 $DirFil/${SRX}.sorted.bam $DIR"
               # commandLine="$operalg $2 $DirFil/${SRX}_k$k.sorted.bam $DIR"
               jobName="OPERA-LG_${SRX}_$TAG"

            else
               echo "!!! ERROR: SCAFFOLDING tool should be \"BESST\" or \"OPERA-LG\" and not: ${SCAFFOLDING_tool} !!!"
               exit 1
            fi

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
               else
                  echo "Job $jobID is still not taken into account"
               fi
            done

            # rm "job_"$jobName".sge"

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

BAM_REC ${MAPPING_dir} $REF

echo "Finishing run at: `date`"