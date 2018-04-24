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
* C++ (≥c++11):
	* [Bio++](http://biopp.univ-montp2.fr/) ([installation instructions](http://biopp.univ-montp2.fr/wiki/index.php/Installation))
	* [Boost](www.boost.org)


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
	* [Treerecs](https://gitlab.inria.fr/Phylophile/Treerecs) (C++ reimplementation of [profileNJ](https://github.com/maclandrol/profileNJ)with Bio++ library)
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
The repository contains a dataset composed of 18 Anopheles genomes corresponding to the dataset used in [2].
Input data files for this dataset are present in the directory [18Anopheles_dataset/](18Anopheles_dataset) and 2 configurations files are present in directory [config_files](config_files) () to apply the pipeline on the 18 Anopheles dataset.
For reproduction of input data of DeCoSTAR used in [2], go to [commit 572d5a5](https://github.com/YoannAnselmetti/DeCoSTAR_pipeline/tree/572d5a50248fa7e0f22c5a8b8dfc52a9fc78275c) of this repository.


## Adapt pipeline to your dataset:
1. Produce data files described in section **INPUT** from your dataset 
2. Create Snakemake configuration file adapted to your dataset from the [configuration file example](config_files/config_example.yaml)
3. Set the correct configuration file path at the top of each "\*.snakefile" files (configfile: "config_files/your_config_file.yaml")


## Command lines for the different steps of the Snakemake pipeline:
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
1. Input data preprocessing to produce GENE file and discard gene families/trees containing included genes (script **preprocessing.snakefile**)
2. Gene trees inference with MUSCLE, GBlocks, RAXML and profileNJ (script **pipeline_trees_inference.snakefile** - TO DO: optional step)
3. Generate scaffolding adjacencies with BESST (script **pipeline_scaffolding.snakefile** - TO DO: optional step) 
4. Generation of input files for DeCoSTAR (script **input_decostar.snakefile**)
5. Execution of DeCoSTAR (script **run_decostar.snakefile** - TO DO)
<!-- 6. Linearisation and stats generation on DeCoSTAR results -->
<!-- 7. Stats graphics generation -->

Steps 2 and 3 are optional. If you don't use the step 2 (pipeline for gene trees inference), you need to provide a gene trees file in the Snakemake configuration file (**families: path_to_your_gene_trees_file**) and to remove the line: **gene_trees: path_to_gene_trees_after_inference_pipeline**



# INPUT
* Species tree file in newick format
* Tab-separated file composed of 2 columns:
	1. Species name
	2. Expected chromosome number
* Gene trees file in newick or NHX format || Gene families tab-separated file composed of 2 columns:
	1. Gene family ID
	2. Gene ID
* GFF files containing exons positions on reference genome assembly of all species present in species tree (see [GFF file format](https://www.ensembl.org/info/website/upload/gff.html))  
=> Column 9 ('attribute' in GFF file format) corresponds to gene ID (must match with gene ID present in gene trees/families)
* Reference genome assemblies of extant species in FASTA file format (name file format: $(species_name)\.\*)
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

It is important that in gene trees/families and gene sequences FASTA file the gene ID is present in the format $(species_name)$separator$(gene_ID). By default we use the '@' as the separator character between the species name and the gene ID. The choice of the separator character used is given in the Snakemake configuration file (line: **separator: '@'**).
Species name should not contain space and is commonly represented with the format: $genus_$species.



# References
[1] Duchemin, W., Anselmetti, Y., Patterson, M., Ponty, Y., Bérard, S., Chauve, C., Scornavacca, C., Daubin, V., Tannier, E. (2017). DeCoSTAR: Reconstructing the Ancestral Organization of Genes or Genomes Using Reconciled Phylogenies. Genome Biology and Evolution, 9(5), 1312–1319.

[2] Anselmetti, Y., Duchemin, W., Tannier, E., Chauve, C., Bérard S. (2018). Phylogenetic signal from rearrangements in 18 Anopheles species by joint scaffolding extant and ancestral genomes. The Sixteenth Asia Pacific Bioinformatics Conference (APBC 2018), Yokohama, Japan, 15-17 January 2018 (to appear in BMC Genomics).
