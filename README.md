Snakemake pipeline to generate input data of the DeCoSTAR software and 
=====

Pipeline to produce input data for DeCoSTAR software [1] to apply ARt-DeCo, ADseq and DeClone algorithms on a dataset from input files whose format is described in the section **INPUT**.



# Requirements
* [Python2](https://www.python.org/.) (≥2.6):
	* [Biopython](http://biopython.org/) 
	* [ETE toolkit](http://etetoolkit.org/)
	* [Matplotlib](https://matplotlib.org/)
	* [Numpy](http://www.numpy.org/)
* [Python3](https://www.python.org/.) (≥3.5):
	* [Snakemake](http://snakemake.readthedocs.io/en/stable/)
* C++:
	* [Bio++](http://biopp.univ-montp2.fr/) ([installation instructions](http://biopp.univ-montp2.fr/wiki/index.php/Installation))
	* [Boost](www.boost.org)




# Libraries/Softwares used (installation with script "install_dependencies.sh"):
* [SRAtoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Samtools](http://samtools.sourceforge.net/)
* [BESST](https://github.com/ksahlin/BESST)
* [MUSCLE](https://www.drive5.com/muscle/)
* [GBlocks](http://molevol.cmima.csic.es/castresana/Gblocks.html)
* [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html)
* [Treerecs](https://gitlab.inria.fr/Phylophile/Treerecs) (=> reimplementation of [profileNJ](https://github.com/maclandrol/profileNJ) in C++ with Bio++ library)
* [DeCoSTAR](http://pbil.univ-lyon1.fr/software/DeCoSTAR/index.html) ([gitHub repository](https://github.com/WandrilleD/DeCoSTAR))



# Data reproduction of the article [2]
The repository contains a dataset composed of 18 Anopheles genomes corresponding to the dataset used in [2].
For reproduction of input data of DeCoSTAR used in [2], go to [commit 572d5a5](https://github.com/YoannAnselmetti/DeCoSTAR_pipeline/tree/572d5a50248fa7e0f22c5a8b8dfc52a9fc78275c) of this repository.



# Quickstart
After downloading the snakemake pipeline with the command line:
```
git clone https://github.com/YoannAnselmetti/DeCoSTAR_pipeline.git
cd DeCoSTAR_pipeline
```
First, step will consist to install all dependencies/softwares necessary to execute the different steps of the pipeline:
```
install_dependencies.sh **(IN PROGRESS)**
```

## Adapt pipeline to your dataset:
1. Provide files from your dataset described in section "INPUT".
2. Create snakemake configuration file for your dataset from "config_files/config_example.yaml"
3. Set the correct configuration file path at the top of each "\*.snakefile" files (configfile: "config_files/your_config_file.yaml").


## Command lines of the different snakemake pipeline steps:
```
snakemake --snakefile preprocessing.snakefile -j <N>
<!-- snakemake --snakefile pipeline_trees_inference.snakefile -j <N> **(TO DO)** -->
<!-- snakemake --snakefile pipeline_scaffolding.snakefile -j <N> **(TO DO)** -->
snakemake --snakefile input_decostar.snakefile -j <N>
snakemake --snakefile run_decostar.snakefile -j <N> **(TO DO)**
<!-- snakemake --snakefile linearization_and_stats.snakefile -j <N> **(TO DO)** -->
<!-- snakemake --snakefile stats_graphs.snakefile -j <N> **(TO DO)** -->
```



# Pipeline
The pipeline to execute DeCoSTAR on a dataset is divided in 5 parts:

1. Pre-process data to produce GENE file and discard gene families/trees containing included genes (script **preprocessing.snakefile**)
2. Gene trees inference with MUSCLE, GBlocks, RAXML and profileNJ (script **pipeline_trees_inference.snakefile** - TO DO: facultative step)
3. Generate scaffolding adjacencies with BESST (script **pipeline_scaffolding.snakefile** - TO DO: facultative step) 
4. Generation of input files for DeCoSTAR (script **input_decostar.snakefile**)
5. Execution of DeCoSTAR (script **run_decostar.snakefile** - TO DO)
<!-- 6. Linearization and stats generation on DeCoSTAR results -->
<!-- 7. Stats graphics generation -->



# INPUT
* Species tree file in newick format
* Tab-separated file composed of 2 columns:
	1. Species_name
	2. Expected chromosome number
* Gene trees file in newick or NHX format || Gene families tab-seperated file composed of 3 columns:
	1. Gene family ID
	2. Gene ID
	3. Species name
* GFF files containing exons positions on reference genome assembly of all species present in species tree (see [GFF file format](https://www.ensembl.org/info/website/upload/gff.html))  
=> Column 9 ('attribute' in GFF file format) corresponds to gene ID (must match with gene ID present in gene trees/families)
* Reference genome assemblies of extant species in FASTA file format
* If scaffolding data are available, user has to provide directory with SRA architecture, ex:
```
FASTQ/RAW/Anopheles_albimanus/
├── SRX084279
│   ├── SRR314655
│   ├── SRR314656
│   └── SRR314659
├── SRX111456
│   ├── SRR389778
│   ├── SRR389781
│   ├── SRR390324
│   └── SRR390326
└── SRX200219
    └── SRR606148
```


# References
[1] Duchemin, W., Anselmetti, Y., Patterson, M., Ponty, Y., Bérard, S., Chauve, C., Scornavacca, C., Daubin, V., Tannier, E. (2017). DeCoSTAR: Reconstructing the Ancestral Organization of Genes or Genomes Using Reconciled Phylogenies. Genome Biology and Evolution, 9(5), 1312–1319.

[2] Anselmetti, Y., Duchemin, W., Tannier, E., Chauve, C., Bérard S. (2018). Phylogenetic signal from rearrangements in 18 Anopheles species by joint scaffolding extant and ancestral genomes. The Sixteenth Asia Pacific Bioinformatics Conference (APBC 2018), Yokohama, Japan, 15-17 January 2018 (to appear in BMC Genomics).
