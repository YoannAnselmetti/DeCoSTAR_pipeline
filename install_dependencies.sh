#/bin/bash

initDIR=`pwd`
echo $initDIR


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





################
### INSTALL SOFTWARES: Bio++, Boost, ...
################
softDIR=$initDIR/bin/softwares
echo $softDIR
mkdir -p $softDIR

#######
### Installation of DeCoSTAR software (GitHub repository)
#######
cd $initDIR/bin/softwares
rm -rf DeCoSTAR
git clone https://github.com/WandrilleD/DeCoSTAR.git
cd DeCoSTAR
# sed -i 's/BPP_INCLUDE=$(HOME)\/local\/bpp\/dev\/include/BPP_INCLUDE=$bppDIR\/include/g' makefile
# sed -i 's/BPP_LIB=$(HOME)\/local\/bpp\/dev\/lib/BPP_INCLUDE=$bppDIR\/lib/g' makefile
make bin/DeCoSTAR


#######
### Installation of Treerecs software (GitHub repository)
#######
cd $initDIR/bin/softwares
rm -rf Treerecs
git clone https://gitlab.inria.fr/Phylophile/Treerecs.git
cd Treerecs
cmake -DCMAKE_BUILD_TYPE=MinSizeRel .
make Treerecs


#######
### Installation of BESST software (GitHub repository)
#######
cd $initDIR/bin/softwares
rm -rf BESST
git clone https://github.com/ksahlin/BESST.git
cd BESST
./runBESST


#######
### Installation of SRA Toolkit (GitHub repository)
#######
cd $initDIR/bin/softwares
rm -rf sra-tools
git clone https://github.com/ncbi/sra-tools.git
cd sra-tools
rm Trimmomatic-Src-0.36.zip


#######
### Installation of Bowtie2 software (GitHub repository)
#######
cd $initDIR/bin/softwares
git clone https://github.com/BenLangmead/bowtie2.git
cd bowtie2
make


#######
### Installation of MUSCLE software (v3.8.1551)
#######
cd $initDIR/bin/softwares
mkdir -p muscle_src_3.8.1551
cd muscle_src_3.8.1551
wget https://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz
tar -zxvf muscle_src_3.8.1551.tar.gz
rm muscle_src_3.8.1551.tar.gz
make


#######
### Installation of Trimmomatic software (v0.36)
#######
cd $initDIR/bin/softwares
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-Src-0.36.zip
unzip Trimmomatic-Src-0.36.zip
rm Trimmomatic-Src-0.36.zip







