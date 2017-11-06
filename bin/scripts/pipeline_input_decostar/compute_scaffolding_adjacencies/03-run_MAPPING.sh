#!/bin/bash
###
###   Goal:
###      Script to map reads on genome reference on cluster with qsub
###      and transform SAM files in BAM files (indexed and compressed format)
###
###   INPUT:
###      1- MAPPING tool (Bowtie2, BWA or mrsFAST)
###      2- PATH of INPUT directory containing FASTQ files
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/FASTQ/TRIMMOMATIC3/ALL/Anopheles_albimanus)
###      3- OUTPUT directory name (not complete PATH)
###         (Bowtie2_ALL / Bowtie2_k100)
###         (BWA_ALL / BWA_k100)
###         (mrsFAST_ALL / mrsFAST_k100)
###      4- INDEX of reference genome
###         (/share/nas-isem_i/yanselmetti/DATA_SEQ/INDEX/Bowtie2/...)
###      5- Species to process
###         (Anopheles_albimanus)
###      6- TAG to allow severa experiments without interfering with jobs
###         of others experiments (Aalb_TRIMMOMATIC3_ALL)
###      7- Multiple alignment allowed for MAPPING (If not given then 50 alignments are considered)
###         (50)
###
###   Name: run_MAPPING.sh                     Author: Yoann Anselmetti
###   Creation date: 2015/08/26                Last modification: 2017/03/15
###

MAPPING_tool=$1
FASTQ_dir=$2
MAP_dir=$3
INDEX=$4
species=$5
TAG=$6

job_max=100
k=50          # Nb of multiple alignments consider by Bowtie2
boolK=true    # IF false: consider all alignments with Bowtie2 if true: consider $k best alignments for Bowtie2 

re='^[0-9]+$'
if [ -z $7 ]; then
   echo "Default usage"
   echo $k" multiple alignments will be considered for mapping of reads on reference genome"
else
   if [ $7 = "a" ]; then
      echo "All multiple alignments will be considered for mapping of reads on reference genome"
      boolK=false
   elif [[ $7 =~ $re ]]; then
      echo $k" multiple alignments will be considered for mapping of reads on reference genome"
      boolK=true
      k=$7
   else
      echo "ERROR, parameter 7 should be an integer or equal to \"a\" (if want to consider all alignments)!!!"
      exit
   fi
fi


file_orientation="/share/nas-isem_i/yanselmetti/DATA_SEQ/orientation_libraries"
mrsFAST="/home/yanselmetti/Software/MAPPING/mrsfast-3.3.9/mrsfast"
bowtie2="/home/yanselmetti/Software/MAPPING/bowtie2-2.2.9/bowtie2" # IF CHANGE NAME OF READ, THIS ONE CAN BE MAPPED TO AN OTHER POSITION IF MULTIPLE ALIGNMENT (Firtiea & Alkan, 2016)
bwa="/home/yanselmetti/Software/MAPPING/bwa-0.7.15/bwa" # NOT RECOMMANDED SINCE BWA HEVE BEEN SHOWN TO BE NON-DETERMINISTIC METHODS (Firtiea & Alkan, 2016)
samtools="/home/yanselmetti/Software/MAPPING/samtools-1.3.1/samtools"

headSGE=#\!"/bin/bash\n#$ -V\n#$ -S /bin/bash\n#$ -cwd\n#$ -l h_rt=200000:00:00,mem_free=8G\n#$ -N "


declare -A SRX_ori
declare -A SRX_insert

while read line; do
   SRX=$(echo $line | cut -d" " -f2)
#   echo $SRX
   ori=$(echo $line | cut -d" " -f3)
#   echo $ori
   insert=$(echo $line | cut -d" " -f4)
#   echo $insert
   SRX_ori+=([$SRX]=$ori)
   SRX_insert+=([$SRX]=$insert)
   # SRX_ori+=([$SRX]=$ori":"$insert)
done < ${file_orientation}



function FASTQ_REC {
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
         DIR=$(echo $1/$SRR | sed "s/FASTQ/MAPPING\/${MAP_dir}/g")
         echo -e "\tmkdir -p $DIR"
         mkdir -p $DIR

         SRX=$(echo $1 | awk -F/ '{print $(NF)}')
         insertSize=${SRX_insert["$SRX"]}
         echo -e "\tInsert size = "$insertSize

         # WAIT until a slot is available
         jobs_nb=$(qstat | grep $USER | wc -l)
         while [ ${jobs_nb} -ge ${job_max} ]; do
            echo "No spot yet"
            sleep 10
            jobs_nb=$(qstat | grep $USER | wc -l)
         done

         commandLine=""
         jobName=""
         if [[ ${MAPPING_tool} == *"Bowtie2"* ]]; then
            ###################
            ### Bowtie2 run ###
            ###################
            # Get all alignments of each reads on genome assembly (option -a) => !!! TAKE LOT OF TIME!!! -> Bowtie2 is not designed or this
            # If want less use option -k $k where $k is the number of best alignment taken into account by Bowtie2

            if $boolK; then
               commandLine="$bowtie2 -t -x $2/$3/$3 -k $k -1 $f1 -2 $f2 -S $DIR/${SRR}_k$k.sam;\n$samtools view -bS $DIR/${SRR}_k$k.sam -o $DIR/${SRR}_k$k.bam;\nsleep 10;\nif [[ -f $DIR/${SRR}_k$k.bam ]]; then rm $DIR/${SRR}_k$k.sam fi;"
            else
               commandLine="$bowtie2 -t -x $2/$3/$3 -a -1 $f1 -2 $f2 -S $DIR/$SRR.sam;\n$samtools view -bS $DIR/$SRR.sam -o $DIR/$SRR.bam;\nsleep 10;\nif [[ -f $DIR/$SRR.bam ]]; then\n\trm $DIR/$SRR.sam\nfi;"
            fi
            jobName="${MAP_dir}_${SRR}_$TAG"

         elif [[ ${MAPPING_tool} == *"mrsFAST"* ]]; then
            ###################
            ### mrsFAST run ###
            ###################

            if [[ $insertSize == "180" ]]; then
               commandLine="$mrsFAST --search $2/$3/$3 --mem 4 --threads 1 --pe --seqcomp --min 0 --max 500 --seq1 $f1 --seq2 $f2 --disable-nohits -o $DIR/$SRR.sam;\n$samtools view -bS $DIR/$SRR.sam -o $DIR/$SRR.bam;\nsleep 1;\nif [[ -f $DIR/$SRR.bam ]]; then\n\trm $DIR/$SRR.sam\nfi;"
            elif [[ $insertSize == "1500" ]]; then
               commandLine="$mrsFAST --search $2/$3/$3 --mem 4 --threads 1 --pe --seqcomp --min 0 --max 6000 --seq1 $f1 --seq2 $f2 --disable-nohits -o $DIR/$SRR.sam;\n$samtools view -bS $DIR/$SRR.sam -o $DIR/$SRR.bam;\nsleep 1;\nif [[ -f $DIR/$SRR.bam ]]; then\n\trm $DIR/$SRR.sam\nfi;"
            elif [[ $insertSize == "38000" ]]; then
               commandLine="$mrsFAST --search $2/$3/$3 --mem 4 --threads 1 --pe --seqcomp --min 10000 --max 80000 --seq1 $f1 --seq2 $f2 --disable-nohits -o $DIR/$SRR.sam;\n$samtools view -bS $DIR/$SRR.sam -o $DIR/$SRR.bam;\nsleep 1;\nif [[ -f $DIR/$SRR.bam ]]; then\n\trm $DIR/$SRR.sam\nfi;"
            else
               echo "ERROR, insert size: "$insertSize" doesn't correspond to known value."
               exit 1
            fi
            jobName="${MAP_dir}_${SRR}_$TAG"

         elif [[ ${MAPPING_tool} == *"BWA"* ]]; then
            ###############
            ### BWA run ###
            ###############

            commandLine="$bwa mem $2/$3/$3 $f1 $f2 > $DIR/$SRR.sam;\n$samtools view -bS $DIR/$SRR.sam -o $DIR/$SRR.bam;\nsleep 1;\nif [[ -f $DIR/$SRR.bam ]]; then\n\trm $DIR/$SRR.sam\nfi;"
            jobName="${MAP_dir}_${SRR}_$TAG"

         else
            echo "ERROR: mapping tool should be \"mrsFAST\", \"BWA\" or \"Bowtie2\" !!!"
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
            fi
         done

         # rm "job_"$jobName".sge"

      elif [ -d "$1/$DirFil" ]; then
         echo "DIR $1/$DirFil"
         FASTQ_REC $1/$DirFil $2 $3
      fi
   done
}

cd $PBS_O_WORKDIR
echo "Current working directory is `pwd`"
echo "Running on hostname `hostname`"
echo "Starting run at: `date`"

FASTQ_REC ${FASTQ_dir} $INDEX $species

echo "Finishing run at: `date`"