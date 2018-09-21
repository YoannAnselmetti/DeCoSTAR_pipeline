#/bin/bash

initDIR=`pwd`
echo $initDIR


################
### INSTALL SOFTWARE: DeCoSTAR, Treerecs, ...
################
softlibDIR=$initDIR/bin/software_libraries
echo $softlibDIR
if [ ! -d $softlibDIR/$softlibDIR ]; then
	mkdir -p $softlibDIR
fi


############
### SOFTWARE WITH GITHUB REPOSITORY
############


#######
### Installation of DeCoSTAR software (GitHub repository)
#######
echo -e "\n\nINSTALL DeCoSTAR:\n"
cd $softlibDIR
if [ -d $softlibDIR/DeCoSTAR ]; then
	rm -rf $softlibDIR/DeCoSTAR
fi
git clone https://github.com/WandrilleD/DeCoSTAR.git
cd DeCoSTAR
# sed -i 's/BPP_INCLUDE=$(HOME)\/local\/bpp\/dev\/include/BPP_INCLUDE=$bppDIR\/include/g' makefile
# sed -i 's/BPP_LIB=$(HOME)\/local\/bpp\/dev\/lib/BPP_INCLUDE=$bppDIR\/lib/g' makefile
make bin/DeCoSTAR


#######
### Installation of Treerecs software (GitHub repository)
#######
echo -e "\n\nINSTALL Treerecs:\n"
cd $softlibDIR
if [ -d $softlibDIR/Treerecs ]; then
	rm -rf $softlibDIR/Treerecs
fi
git clone https://gitlab.inria.fr/Phylophile/Treerecs.git
cd Treerecs
cmake -DCMAKE_BUILD_TYPE=Release .
make


#######
### Installation of BESST software (GitHub repository)
#######
# Requirements: pysam, networkx, numpy and scipy
echo -e "\n\nINSTALL BESST:\n"
cd $softlibDIR
if [ -d $softlibDIR/BESST ]; then
	rm -rf $softlibDIR/BESST
fi
git clone https://github.com/ksahlin/BESST.git
# cd BESST
# ./runBESST


#######
### Installation of SRA Toolkit (GitHub repository)
#######
# ERRO: ngs-sdk not found
echo -e "\n\nINSTALL SRA toolkit:\n"
cd $softlibDIR
if [ -d $softlibDIR/sra-tools ]; then
	rm -rf $softlibDIR/sra-tools
fi
git clone https://github.com/ncbi/sra-tools.git
cd sra-tools
./configure



#######
### Installation of Bowtie2 software (GitHub repository)
#######
echo -e "\n\nINSTALL Bowtie2:\n"
cd $softlibDIR
if [ -d $softlibDIR/bowtie2 ]; then
	rm -rf $softlibDIR/bowtie2
fi
git clone https://github.com/BenLangmead/bowtie2.git
cd bowtie2
make NO_TBB=1


#######
### Installation of RAxML software (GitHub repository)
#######
echo -e "\n\nINSTALL RAxML:\n"
cd $softlibDIR
if [ -d $softlibDIR/standard-RAxML ]; then
	rm -rf $softlibDIR/standard-RAxML
fi
git clone https://github.com/stamatak/standard-RAxML.git
cd standard-RAxML
make -f Makefile.gcc
make -f Makefile.AVX2.gcc
make -f Makefile.SSE3.gcc
rm *.o


#######
### Installation of htslib library and SAMtools software (GitHub repository)
#######
echo -e "\n\nINSTALL htslib library:\n"
cd $softlibDIR
if [ -d $softlibDIR/htslib ]; then
	rm -rf $softlibDIR/htslib
fi
git clone https://github.com/samtools/htslib.git
cd htslib
autoheader
autoconf -Wno-synta
./configure
make
cd $softlibDIR
if [ -d $softlibDIR/samtools ]; then
	rm -rf $softlibDIR/samtools
fi
git clone https://github.com/samtools/samtools.git
cd samtools
autoheader
autoconf -Wno-synta
./configure
make



############
### SOFTWARE WITHOUT GITHUB REPOSITORY
############


#######
### Installation of MUSCLE software (v3.8.1551)
#######
echo -e "\n\nINSTALL MUSCLE:\n"
cd $softlibDIR
if [ -d $softlibDIR/muscle_src_3.8.1551 ]; then
	rm -rf $softlibDIR/muscle_src_3.8.1551
fi
mkdir -p muscle_src_3.8.1551
cd muscle_src_3.8.1551
wget https://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz
tar -zxvf muscle_src_3.8.1551.tar.gz
rm muscle_src_3.8.1551.tar.gz
make


#######
### Installation of Trimmomatic software (v0.36)
#######
echo -e "\n\nINSTALL Trimmomatic:\n"
cd $softlibDIR
if [ -d $softlibDIR/Trimmomatic-Src-0.36 ]; then
	rm -rf $softlibDIR/Trimmomatic-Src-0.36
fi
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-Src-0.36.zip
unzip Trimmomatic-Src-0.36.zip
rm Trimmomatic-Src-0.36.zip


#######
### Installation of Gblocks software (v0.91b)
#######
echo -e "\n\nINSTALL GBlocks:\n"
cd $softlibDIR
if [ -d $softlibDIR/Gblocks_0.91b ]; then
	rm -rf $softlibDIR/Gblocks_0.91b
fi
wget http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z
uncompress Gblocks_Linux64_0.91b.tar.Z
tar -xvf Gblocks_Linux64_0.91b.tar
rm Gblocks_Linux64_0.91b.tar





############
### OTHERS POTENTIAL SOFTWARE
############

### ALIGNMENT TRIMMING
# HMMCleaner
# BMGE
# GUIDANCE
# FAS
# ZORRO
# PSAR
# Noisy
# trimAl

### MULTIPLE SEQURENCE ALIGNMENT
# MACSE
# MAFFT

### TREE INFERENCE
# IQ-TREE


### READS TRIMMING
# Skewer
# Trim Galore

### MAPPING
# mrsFAST-Ultra (required all reads have the same size) / mrFAST => different tools!!!
# BWA-MEM
# RazerS3

### SCAFFOLDING
# OPERA-LG





################
### INSTALL LIBRARIES: Bio++, Boost, ...
################
# libDIR=$initDIR/bin/libraries
# echo $libDIR
# mkdir -p $libDIR

# #######
# ### Installation of Bio++ library
# #######
# bppDIR=$libDIR/bpp/dev
# rm -rf $bppDIR
# mkdir -p $bppDIR/sources
# cd $bppDIR/sources
# echo -e "\ngit clone https://github.com/BioPP/bpp-core"
# git clone https://github.com/BioPP/bpp-core
# echo -e "\ngit clone https://github.com/BioPP/bpp-seq"
# git clone https://github.com/BioPP/bpp-seq
# echo -e "\ngit clone https://github.com/BioPP/bpp-popgen"
# git clone https://github.com/BioPP/bpp-popgen
# echo -e "\ngit clone https://github.com/BioPP/bpp-phyl"
# git clone https://github.com/BioPP/bpp-phyl

# for elem in `ls $bppDIR/sources`; do
# 	echo -e "\nBuild $elem library:"
# 	cd $bppDIR/sources/$elem
# 	mkdir build
# 	cd build
# 	cmake -DCMAKE_INSTALL_PREFIX=$bppDIR .. # prepare compilation
# 	make # compile
# 	make install # move files to the installation directory
# done


#######
### Installation of Boost library
#######
# NOT SUPPORTED

