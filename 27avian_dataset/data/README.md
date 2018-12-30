# README content:
This README contains informations on the content of this directory that corresponds to input data used to produce the results of DeCoSTAR for scaffolding improvment of 27 avian dataset described in Leroy et al., in prep.
It contains also some command lines necessary to download genome assemblies FASTA files and preprocess data to the DeCoSTAR input files format to reproduce DeCoSTAR scaffolding of the 27 avian dataset. 



# Directory content:
This directory contains all input data files necessary to the execution of DeCoSTAR for scaffolding improvment of 26 avian reference genomes (see list below) and the newly sequenced reference genome of Zosterops borbonicus.  

List of the 26 species name present in this study (in addition to Zosterops borbonicus) with the corresponding RefSeq assembly accesion and assemblies ID used:
Acanthisitta_chloris	GCF_000695815.1	ASM69581v1
Anas_platyrhynchos	GCF_000355885.1	BGI_duck_1.0
Aptenodytes_forsteri	GCF_000699145.1	ASM69914v1
Calypte_anna	GCF_000699085.1	ASM69908v1
Chaetura_pelagica	GCF_000747805.1	ChaPel_1.0
Charadrius_vociferus	GCF_000708025.1	ASM70802v2
Columba_livia	GCF_000337935.1	Cliv_1.0
Corvus_brachyrhynchos	GCF_000691975.1	ASM69197v1
Corvus_cornix	GCF_000738735.2	ASM73873v2
Cuculus_canorus	GCF_000709325.1	ASM70932v1
Egretta_garzetta	GCF_000687185.1	ASM68718v1
Falco_peregrinus	GCF_000337955.1	F_peregrinus_v1.0
Ficedula_albicollis	GCF_000247815.1	FicAlb1,5
Gallus_gallus	GCF_000002315.5	GRCg6a
Geospiza_fortis	GCF_000277835.1	GeoFor_1.0
Haliaeetus_leucocephalus	GCF_000737465.1	Haliaeetus_leucocephalus-4.0
Manacus_vitellinus	GCF_001715985.1	ASM171598v1
Melopsittacus_undulatus	GCF_000238935.1	Melopsittacus_undulatus_6.3
Nipponia_nippon	GCF_000708225.1	ASM70822v1
Opisthocomus_hoazin	GCF_000692075.1	ASM69207v1
Picoides_pubescens	GCF_000699005.1	ASM69900v1
Pygoscelis_adeliae	GCF_000699105.1	ASM69910v1
Struthio_camelus_australis	GCF_000698965.1	ASM69896v1
Taeniopygia_guttata	GCF_000151805.1	Taeniopygia_guttata-3.2.4
Tinamus_guttatus	GCF_000705375.1	ASM70537v2
Zonotrichia_albicollis	GCF_000385455.1	Zonotrichia_albicollis-1.0.1



## Directory architecture:
.
├── DATA_SEQ/
├── INPUT_DATA/
│   ├── 27avian_species_tree.nwk
│   ├── 27avian_species.txt
│   ├── 27avian_trees.nwk
│   ├── chromosome_map/Zosterops_borbonicus/correspondance_syntenicblocks_Zobov2_ZeFi_final_version_assigned.agp
│   ├── FASTA/SCAFF/
│   └── RAW_GFF/
└── README.md



### "DATA_SEQ/" directory
This directory contains scored scaffolding links obtained with BESST software (Sahlin et al., 2014) to produce scaffolding adjacencies that will be used as input of DeCoSTAR to bring additional adjacencies with paired-sequencing data support.
These data have been produced with the following pipeline:
1. Reads trimming with Trimmomatic
2. Mapping of trimmed reads on reference genomes with Bowtie2
3. Execute BESST on reads mapping to produce scored scaffolding links between scaffolds of reference genomes 

List of species for which NO paired-sequencing data were available:
- Taeniopygia guttata



### "INPUT_DATA/" directory
This directory contains 4 files and 1 directory:
1. 27avian_species_tree.nwk	-> species tree file in newick format
2. 27avian_species.txt		-> association between species name and chromosome number (haplome)
3. 27avian_trees.nwk		-> gene trees in newick format
4. chromosome_map/Zosterops_borbonicus/correspondance_syntenicblocks_Zobov2_ZeFi_final_version_assigned.agp	-> assignment of Zosterops borbonicus scaffolds on Taeniopygia guttata chromosome by pairwise alignment (LASTZ)
5. FASTA/SCAFF/				-> directory that will contain the 27 avian genome assembly FASTA files
6. RAW_GFF/					-> directory containing GFF files of the 26 avian species corresponding to the genome assembly ID used + newly generated GFF file of Zosterops borbonicus


#### "FASTA/SCAFF/" directory
In this directory will stored the reference genome assembly FASTA files for the 27 avian species.
You first need to download the Zosterops borbonicus genome assembly FASTA file in this directory (Linux):
````
wget LINK_TO_THE_ZOBO_FASTA_FILE FASTA/SCAFF/
````
For the 26 additional species you have to run the following command lines (Linux) from the directory where this README file is stored:
````
cd INPUT_DATA/FASTA/SCAFF
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/695/815/GCF_000695815.1_ASM69581v1/GCF_000695815.1_ASM69581v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/355/885/GCF_000355885.1_BGI_duck_1.0/GCF_000355885.1_BGI_duck_1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/699/145/GCF_000699145.1_ASM69914v1/GCF_000699145.1_ASM69914v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/699/085/GCF_000699085.1_ASM69908v1/GCF_000699085.1_ASM69908v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/747/805/GCF_000747805.1_ChaPel_1.0/GCF_000747805.1_ChaPel_1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/708/025/GCF_000708025.1_ASM70802v2/GCF_000708025.1_ASM70802v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/337/935/GCF_000337935.1_Cliv_1.0/GCF_000337935.1_Cliv_1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/691/975/GCF_000691975.1_ASM69197v1/GCF_000691975.1_ASM69197v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/738/735/GCF_000738735.2_ASM73873v2/GCF_000738735.2_ASM73873v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/709/325/GCF_000709325.1_ASM70932v1/GCF_000709325.1_ASM70932v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/687/185/GCF_000687185.1_ASM68718v1/GCF_000687185.1_ASM68718v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/337/955/GCF_000337955.1_F_peregrinus_v1.0/GCF_000337955.1_F_peregrinus_v1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/247/815/GCF_000247815.1_FicAlb1.5/GCF_000247815.1_FicAlb1.5_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.5_GRCg6a/GCF_000002315.5_GRCg6a_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/277/835/GCF_000277835.1_GeoFor_1.0/GCF_000277835.1_GeoFor_1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/737/465/GCF_000737465.1_Haliaeetus_leucocephalus-4.0/GCF_000737465.1_Haliaeetus_leucocephalus-4.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/715/985/GCF_001715985.1_ASM171598v1/GCF_001715985.1_ASM171598v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/935/GCF_000238935.1_Melopsittacus_undulatus_6.3/GCF_000238935.1_Melopsittacus_undulatus_6.3_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/708/225/GCF_000708225.1_ASM70822v1/GCF_000708225.1_ASM70822v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/692/075/GCF_000692075.1_ASM69207v1/GCF_000692075.1_ASM69207v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/699/005/GCF_000699005.1_ASM69900v1/GCF_000699005.1_ASM69900v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/699/105/GCF_000699105.1_ASM69910v1/GCF_000699105.1_ASM69910v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/698/965/GCF_000698965.1_ASM69896v1/GCF_000698965.1_ASM69896v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/151/805/GCF_000151805.1_Taeniopygia_guttata-3.2.4/GCF_000151805.1_Taeniopygia_guttata-3.2.4_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/705/375/GCF_000705375.1_ASM70537v2/GCF_000705375.1_ASM70537v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/385/455/GCF_000385455.1_Zonotrichia_albicollis-1.0.1/GCF_000385455.1_Zonotrichia_albicollis-1.0.1_genomic.fna.gz
gunzip *gz
cd ../../..
````


#### "RAW_GFF/" directory
To adapt GFF files to the format used by DeCoSTAR, you need to execute the following command lines (Linux): 
````
indir="INPUT_DATA/RAW_GFF"
outdir="INPUT_DATA/GFF"
gunzip $indir"/"*
mkdir -p $outdir
for gff in $(ls $indir); do
	species=$(echo $gff|cut -d. -f1)
	echo $species
	grep -P "\tCDS\t" $indir"/"$gff | grep "Name" | awk -v spe=$species 'BEGIN {OFS="\t"} {print $1,spe,$3,$4,$5,$6,$7,$8}' > buffer1
	grep -P "\tCDS\t" $indir"/"$gff | grep "Name" | cut -d";" -f4|awk -F= '{print $NF}' > buffer2
	paste buffer1 buffer2 > buffer3
	awk '{if ($3=="CDS") print $0}' buffer3 > $outdir"/"$gff
	rm buffer1 buffer2 buffer3
done;

gff="INPUT_DATA/RAW_GFF/Zosterops_borbonicus.gff"
species="Zosterops_borbonicus"
echo $species
grep CDS $gff |awk -v spe=$species 'BEGIN {OFS="\t"} {print $1,spe,"CDS",$4,$5,$6,$7,$8}' > buffer1
grep CDS $gff |awk -F= '{print $2}'|awk -F\; '{print $1}' > buffer2
paste buffer1 buffer2 > $outdir"/Zosterops_borbonicus.gff"
rm buffer1 buffer2
````

Now you can go to the head of this GitHub repository.
Check that the 4 snakefiles upload the 27avian dataset snakemake configuration file, by uncommenting line:
````
configfile: "config_files/snakemake/config_27avian.yaml"
````
Then execute the 4 snakefile scripts:
```
snakemake --snakefile preprocessing.snakefile
snakemake --snakefile input_decostar.snakefile
snakemake --snakefile run_decostar.snakefile
snakemake --snakefile create_adjacencies_graph.snakefile
```
