
Snakemake pipeline to generate input data of the DeCoSTAR software
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



# Quickstart
The current repository version was created to reproduce the input data of DeCoSTAR for the 18 Anopheles dataset used in [2].
This current version don't allow to reproduce the inference of gene trees with "MUSCLE->GBlocks->RAXML->profileNJ" and to generate scaffolding adjacencies with "Trimmomatic->Bowtie2->BESST". However, the scripts to generate these data are respectively available in the directories "bin/scripts/pipeline_input_decostar/pipeline_gene_trees_inference/" and "bin/scripts/pipeline_input_decostar/pipeline_scaffolding_adjacencies/".

Below are the command lines to execute the two steps of the pipeline to produce the input data of the 18 Anopheles dataset for [2]: 

```
snakemake --snakefile preprocessing.snakefile -j <N>
snakemake --snakefile input_decostar.snakefile -j <N>
```

After execution of the snakemake pipeline for the 2 configurations files:
* "config_files/config_input_18Anopheles_WGtopo.yaml"
* "config_files/config_input_18Anopheles_Xtopo.yaml"

In the directory "data/data_DeCoSTAR", DeCoSTAR configuration files have to be handwritten and content of the directory has to be manually reorganised to correspond at the architecture of "data/data_DeCoSTAR" directory in the GitHub repository [ADseq-Anopheles-APBC2018
](https://github.com/YoannAnselmetti/ADseq-Anopheles-APBC2018) (repository of [2]).



# Pipeline
The pipeline to produce input data for DeCoSTAR is divided in 2 steps:
1. Pre-process data to produce GENE file and discard gene families/trees containing included genes (script **preprocessing.snakefile**)
2. Generation of input files for DeCoSTAR after gene trees inference and scaffolding adjacencies pipelines (script **input_decostar.snakefile**)



# INPUT
* Species tree file in newick format
* File with species name and expected chromosome number file (format: "species_name\t#chromosome")
* Gene trees file in newick or NHX format
* File linking species name and gene ID prefix (see file: "data/INPUT_DATA/name_geneID_18Anopheles")
* GFF file containing exons positions on reference genome assemblies (see [GFF file format](https://www.ensembl.org/info/website/upload/gff.html))  
=> Column 9 ('attribute' in GFF file format) corresponds to gene ID (must match with gene ID present in gene trees)
* Reference genome assemblies of extant genomes in FASTA file format



# References
[1] Duchemin, W., Anselmetti, Y., Patterson, M., Ponty, Y., Bérard, S., Chauve, C., Scornavacca, C., Daubin, V., Tannier, E. (2017). DeCoSTAR: Reconstructing the Ancestral Organization of Genes or Genomes Using Reconciled Phylogenies. Genome Biology and Evolution, 9(5), 1312–1319.

[2] Anselmetti, Y., Duchemin, W., Tannier, E., Chauve, C., Bérard S. (2018). Phylogenetic signal from rearrangements in 18 Anopheles species by joint scaffolding extant and ancestral genomes. The Sixteenth Asia Pacific Bioinformatics Conference (APBC 2018), Yokohama, Japan, 15-17 January 2018 (to appear in BMC Genomics).
