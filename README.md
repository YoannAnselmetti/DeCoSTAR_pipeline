Snakemake pipeline for a friendly-use of the DeCoSTAR software
=====

This repository contains a pipeline to produce input data for the DeCoSTAR software [1] in order to apply ARt-DeCo [2], ADseq [3] and DeClone [4] algorithms on a dataset from input files whose format is described in the section **INPUT**.



# Requirements
* [Python2](https://www.python.org/.) (≥v.2.6) -> v.2.7:
	* [Biopython](http://biopython.org/) -> v.1.70
	* [ETE toolkit](http://etetoolkit.org/) -> v.3.1.1
	* [Matplotlib](https://matplotlib.org/) -> v.2.1.1
	* [Numpy](http://www.numpy.org/) -> v.1.14.0
	* [NetworkX](https://networkx.github.io/) -> v.1.11
	* [PyGraphviz](https://pygraphviz.github.io/) -> v.1.3.1
* [Python3](https://www.python.org/.) (≥v.3.5) -> v.3.5.2:
	* [Snakemake](http://snakemake.readthedocs.io/en/stable/) -> v.3.5.5
* C++ (≥c++11):
	* [Bio++](http://biopp.univ-montp2.fr/) ([installation instructions](http://biopp.univ-montp2.fr/wiki/index.php/Installation)) -> v.2.4.0
	* [Boost](www.boost.org) -> v.1.58


# Software/Tools
* Gene trees inference pipeline:
	* [SRAtoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/) ([GitHub repository](https://github.com/ncbi/sra-tools))
	* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
	* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) ([GitHub repository](https://github.com/BenLangmead/bowtie2))
	* [Samtools](http://samtools.sourceforge.net/) ([GitHub repository](https://github.com/samtools/samtools))
	* [BESST](https://github.com/ksahlin/BESST)
* Scaffolding adjacencies pipeline:	
	* [MUSCLE](https://www.drive5.com/muscle/)
	* [GBlocks](http://molevol.cmima.csic.es/castresana/Gblocks.html)
	* [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html) ([GitHub repository](https://github.com/stamatak/standard-RAxML))
	* [Treerecs](https://gitlab.inria.fr/Phylophile/Treerecs) (C++ reimplementation of [profileNJ](https://github.com/maclandrol/profileNJ) with Bio++ library)
* [DeCoSTAR](http://pbil.univ-lyon1.fr/software/DeCoSTAR/index.html) ([GitHub repository](https://github.com/WandrilleD/DeCoSTAR))



# Quickstart
After downloading the Snakemake pipeline with the command line:
```
git clone https://github.com/YoannAnselmetti/DeCoSTAR_pipeline.git
cd DeCoSTAR_pipeline
```
First, step will consist to install all software necessary to execute the different steps of the pipeline:
```
./install_dependencies.sh
```


## Dataset example on 18 Anopheles species
The repository contains a dataset composed of 18 Anopheles genomes corresponding to the dataset used in [3].
Input data files for this dataset are present in the directory [18Anopheles_dataset/](18Anopheles_dataset) and 2 configurations files are present in directory [config_files/snakemake](config_files/snakemake) ("config_18Anopheles_WGtopo.yaml" and "config_18Anopheles_Xtopo.yaml") to apply the pipeline on the 18 Anopheles dataset with 2 species tree topologies.
For reproduction of input data of DeCoSTAR used in [3], go to [commit 572d5a5](https://github.com/YoannAnselmetti/DeCoSTAR_pipeline/tree/572d5a50248fa7e0f22c5a8b8dfc52a9fc78275c) of this repository.


## Adapt pipeline to your dataset:
1. Produce data files described in section **INPUT** from your dataset 
2. Create Snakemake configuration file adapted to your dataset from the [configuration file example](config_files/snakemake/config_example.yaml)
3. Set the correct configuration file path at the top of each "\*.snakefile" files (configfile: "config_files/snakemake/your_config_file.yaml")


## Command lines for the different steps of the Snakemake pipeline (<N>: #CPUs):
```
snakemake --snakefile preprocessing.snakefile -j <N>
snakemake --snakefile input_decostar.snakefile -j <N>
snakemake --snakefile run_decostar.snakefile -j <N>
snakemake --snakefile create_adjacencies_graph.snakefile -j <N>
```



# Pipeline
The pipeline to execute DeCoSTAR on a dataset is divided in 6 parts:
1. Input data preprocessing to produce GENE file and discard gene families/trees containing included genes (script **preprocessing.snakefile**)
2. Gene trees inference with MUSCLE, GBlocks, RAXML and profileNJ (script **pipeline_trees_inference.snakefile** - TO DO: optional step)
3. Generate scaffolding adjacencies with BESST (script **pipeline_scaffolding.snakefile** - TO DO: optional step) 
4. Generation of input files for DeCoSTAR (script **input_decostar.snakefile**)
5. Execution of DeCoSTAR and linearization of adjacencies prediction (script **run_decostar.snakefile**)
6. Generation of adjacencies graph (script **create_adjacencies_graph.snakefile**)

Steps 2 and 3 are optional. If you don't use the step 2 (pipeline for gene trees inference), you need to provide a gene trees file in the Snakemake configuration file (**families: path_to_your_gene_trees_file**) and to remove the line: **gene_trees: path_to_gene_trees_after_inference_pipeline**



# INPUT
The pipeline required 6 input:
* A pecies tree file in newick format
* A tab-separated file composed of 2 columns (chromosome file):
	1. Species name
	2. Expected chromosome number
* A gene trees file in newick or NHX format || A gene families tab-separated file composed of 2 columns (required to use step 2 - script in progress):
	1. Gene family ID
	2. Gene ID
* GFF files containing exons positions on reference genome assembly of all species present in species tree (see [GFF file format](https://www.ensembl.org/info/website/upload/gff.html))  
=> Column 9 ('attribute' in GFF file format) corresponds to gene ID (must match with gene ID present in gene trees/families)
* The eference genome assemblies of extant species in FASTA file format (name file format: $(species_name)\.\*)
* If scaffolding data are available (used in step 3: optional), user has to provide directory with SRA architecture, ex:
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

It is important that in gene trees/families and gene sequences FASTA file the gene ID is present in the format: $(species_name)$separator$(gene_ID). By default we use '@' as separator character between the species name and the gene ID. The choice of the separator character used is given in the Snakemake configuration file (line: **separator: '@'**).
Species name should not contain space and is commonly represented with the format: $genus_$species.



# References
[1] Duchemin, W., Anselmetti, Y., Patterson, M., Ponty, Y., Bérard, S., Chauve, C., Scornavacca, C., Daubin, V., Tannier, E. (2017). DeCoSTAR: Reconstructing the Ancestral Organization of Genes or Genomes Using Reconciled Phylogenies. Genome Biology and Evolution, 9(5), 1312–1319.

[2] Anselmetti, Y., Berry, V., Chauve, C., Chateau, A., Tannier, E., & Bérard, S. (2015). Ancestral gene synteny reconstruction improves extant species scaffolding. BMC Genomics, 16(Suppl 10), S11.

[3] Anselmetti, Y., Duchemin, W., Tannier, E., Chauve, C., & Bérard, S. (2018). Phylogenetic signal from rearrangements in 18 Anopheles species by joint scaffolding extant and ancestral genomes. BMC Genomics, 19(S2), 96.

[4] Chauve, C., Ponty, Y., & Zanetti, J. P. P. (2014). Evolution of genes neighborhood within reconciled phylogenies: an ensemble approach. In Lecture Notes in Computer Science (Advances in Bioinformatics and Computational Biology) (Vol. 8826 LNBI, pp. 49–56).
