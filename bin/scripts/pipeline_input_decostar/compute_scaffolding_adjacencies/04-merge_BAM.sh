#!/bin/bash
###
###   Goal:
###       Merge BAM files SRR directories in 1 BAM file of corresponding SRX directory
###
###   INPUT:
###      1- Directory containing BAM files (FULL path!!!)
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/MAPPING/Bowtie2_ALL/RAW/ALL/Anopheles_albimanus)
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/MAPPING/Bowtie2_k100/RAW/ALL/Anopheles_albimanus)
###      2- TAG to allow several experiments without interfering with jobs of others experiments
###         (Aalb_Bowtie2_ALL_RAW_ALL)
###         (Aalb_Bowtie2_k100_RAW_ALL)
###      3- Multiple alignment used for MAPPING (If not given then 50 alignments were considered)
###         (50)
###
###   OUTPUT:
###      - 1 merged BAM file/SRX DIR from all BAM files of SRR DIR
###
###   Name: merge_BAM.sh                Author: Yoann Anselmetti
###   Creation date: 2015/12/07         Last modification: 2017/03/15
###

BAMdir=$1
TAG=$2

job_max=100
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


mem=2
vmem=$(($mem*2))

samtools="/home/yanselmetti/Software/MAPPING/samtools-1.3.1/samtools"

headSGE=#\!"/bin/bash\n#$ -V\n#$ -S /bin/bash\n#$ -cwd\n#$ -l h_rt=200000:00:00,h_vmem="$vmem"G\n#$ -N "


function BAM_REC {
   DirList=$(find $1/* -prune)
   f1=""
   f2=""
   for DirFil in $DirList; do
      if [ -d "$DirFil" ]; then
         echo "DIR $DirFil"
         # If directory is SRX directory
         if [ $(echo $DirFil | grep "SRX") ]; then
            SRX=$(echo $DirFil | awk 'BEGIN{FS="/"};{print $NF}')
            list_SRR=$(ls -1 $DirFil | grep SRR)
            list_bam_to_merge=""
            j=0
            for SRR in $list_SRR; do
               if $boolK; then
                  list_bam_to_merge=${list_bam_to_merge}$DirFil"/"$SRR"/"$SRR"_k"$k".bam "
               else
                  list_bam_to_merge=${list_bam_to_merge}$DirFil"/"$SRR"/$SRR.bam "
               fi
               j=$(($j+1))
            done

            # WAIT until a slot is available
            jobs_nb=$(qstat | grep $USER | wc -l)
            while [ $jobs_nb -ge ${job_max} ]; do
               echo "No spot yet"
               sleep 10
               jobs_nb=$(qstat | grep $USER | wc -l)
            done

            commandLine=""
            jobName=""
	         if [[ $j = 1 ]]; then
               echo -e "\tSRR: $SRR"
               #################################################################################################
               ### If only 1 SRR just copy sorted.bam and sorted.bam.bai from SRR directory in SRX directory ###
               #################################################################################################

               if $boolK; then
                  commandLine="$samtools sort -m "$mem"G -o $DirFil/${SRX}_k$k.sorted.bam $DirFil/$SRR/${SRR}_k$k.bam;\n$samtools index $DirFil/${SRX}_k$k.sorted.bam;\nsleep 10;\nif [[ -f $DirFil/${SRX}_k$k.sorted.bam.bai ]]; then\n\trm -rf $DirFil/SRR*;\nelse\n\techo -e \"ERROR, file:\\n\\t$DirFil/${SRX}_k$k.sorted.bam.bai\\nis NOT present!!!\";\nfi"
               else
                  commandLine="$samtools sort -m "$mem"G -o $DirFil/$SRX.sorted.bam $DirFil/$SRR/$SRR.bam;\n$samtools index $DirFil/$SRX.sorted.bam;\nsleep 10;\nif [[ -f $DirFil/$SRX.sorted.bam.bai ]]; then\n\trm -rf $DirFil/SRR*;\nelse\n\techo -e \"ERROR, file:\\n\\t$DirFil/$SRX.sorted.bam.bai\\nis NOT present!!!\";\nfi"
               fi
               jobName="merge_BAM_${SRX}_$TAG"

            elif [[ $j -gt 1 ]]; then
               ################################################ ###################################################################
               ### IF more than 1 SRR: MERGE BAM files from SRR directories in SRX directory and sort and index merged BAM file ###
               ####################################################################################################################

               if $boolK; then
                  commandLine="$samtools merge -f $DirFil/${SRX}_k$k.bam ${list_bam_to_merge};\n$samtools sort -m "$mem"G -o $DirFil/${SRX}_k$k.sorted.bam $DirFil/${SRX}_k$k.bam;\nsleep 10;\nif [[ -f $DirFil/${SRX}_k$k.sorted.bam ]]; then\n\trm -rf $DirFil/SRR*;\nelse\n\techo -e \"ERROR, file:\\n\\t$DirFil/${SRX}_k$k.sorted.bam\\nis NOT present!!!\";\nfi;\n$samtools index $DirFil/${SRX}_k$k.sorted.bam;\nif [[ -f $DirFil/${SRX}_k$k.sorted.bam ]]; then\n\trm -rf $DirFil/${SRX}_k$k.bam;\nelse\n\techo -e \"ERROR, file:\\n\\t$DirFil/${SRX}_k$k.sorted.bam\\nis NOT present!!!\";\nfi"
               else
                  commandLine="$samtools merge -f $DirFil/$SRX.bam ${list_bam_to_merge};\n$samtools sort -m "$mem"G -o $DirFil/$SRX.sorted.bam $DirFil/$SRX.bam;\nsleep 10;\nif [[ -f $DirFil/$SRX.sorted.bam ]]; then\n\trm -rf $DirFil/SRR*;\nelse\n\techo -e \"ERROR, file:\\n\\t$DirFil/$SRX.sorted.bam\\nis NOT present!!!\";\nfi;\n$samtools index $DirFil/$SRX.sorted.bam;\nif [[ -f $DirFil/$SRX.sorted.bam ]]; then\n\trm -rf $DirFil/$SRX.bam;\nelse\n\techo -e \"ERROR, file:\\n\\t$DirFil/$SRX.sorted.bam\\nis NOT present!!!\";\nfi"
               fi
               jobName="merge_BAM_${SRX}_$TAG"

            else
               echo "NO SRR directory in SRX directory: "$SRX
            fi
            

            if [ ! -z $jobName ]; then
               # Print current job to run
               echo -e "\t\tCommand line:\n"$commandLine"\n"
               # Run current job
               echo -e $headSGE$jobName"\n\n"$commandLine > "job_"$jobName".sge"

               sleep 5

               jobID=$(qsub "job_"$jobName".sge" | cut -d' ' -f 3)
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
            fi

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

BAM_REC $BAMdir

echo "Finishing run at: `date`"