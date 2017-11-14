Pipeline for a friendly use of the DeCoSTAR software from standard format files
=====

Pipeline to produce input data for DeCoSTAR software to apply ARt-DeCo, ADseq and DeClone algorithms on a dataset composed of files to the standard format described in the section **INPUT**.


# Requirements
* [Python2](https://www.python.org/.) (≥2.6):
	* [Biopython](http://biopython.org/) 
	* [ETE toolkit](http://etetoolkit.org/)
	* [Matplotlib](https://matplotlib.org/)
	* [Numpy](http://www.numpy.org/)
* [Python3](https://www.python.org/.) (≥3.5):
	* [Snakemake](http://snakemake.readthedocs.io/en/stable/)
* C++:
	* [Bio++](http://biopp.univ-montp2.fr/wiki/index.php/Main_Page) ([installation intructions](http://biopp.univ-montp2.fr/wiki/index.php/Installation))
	* [Boost](www.boost.org)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [Samtools](http://samtools.sourceforge.net/)
* [DeCoSTAR](http://pbil.univ-lyon1.fr/software/DeCoSTAR/index.html) ([gitHub repository](https://github.com/WandrilleD/DeCoSTAR))
* [profileNJ](https://github.com/maclandrol/profileNJ)
* [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html)



# Quickstart
```
<!--install_DeCoSTAR.sh **(TO DO)**-->
snakemake --snakefile preprocessing.snakefile -j <N>
snakemake --snakefile pipeline_trees_inference.snakefile -j <N> **(In progress)**
snakemake --snakefile pipeline_scaffolding.snakefile -j <N> **(In progress)**
snakemake --snakefile input_decostar.snakefile -j <N>
snakemake --snakefile decostar.snakefile -j <N>
snakemake --snakefile linearization_and_stats.snakefile -j <N>
snakemake --snakefile stats_graphs.snakefile -j <N>
```



# Pipeline
The pipeline to execute DeCoSTAR on a dataset is divided in 7 parts:

1. Preprocess data to produce GENE file and discard geen trees containing included genes (script **preprocessing.snakefileconfig file**)
2. Pipeline inference gene trees with profileNJ (script development in progress)
3. Pipeline to produce scaffolding adjacencies with BESST (facultative step and script development in progress) 
4. Production of input files for DeCoSTAR (script **input_decostar.yaml**)
5. Execution of DeCoSTAR
6. Linearization and stats generation on DeCoSTAR results
7. Stats graphics generation



# INPUT
* species tree file in newick format
* file with species name and expected chromosome number file (format: "species_name\t#chromosome")
* gene trees/families file in newick or NHX format (future developement => only need to give a file with gene cluster)
* GFF file giving informations on gene position on reference genome assemblies (see [GFF file format](https://www.ensembl.org/info/website/upload/gff.html))
* Reference genome assemblies of extant genomes in FASTA file format (to compute scaffolding adjacencies)
* If scaffolding data are available, user has to provide directory with architecture 



<!--# Usage-->

<!-- The directory **bin/scripts/pipeline_input_decostar** contains scripts developed to process raw input data of the 18 Anopheles dataset (located in directory *data/INPUT_DATA*) into input data for the software DeCoSTAR to execute the ARt-DeCo and ADseq algorithms. All the data produced by these scripts are located in the directory *data/*. In this directory, the script names are prefixed by a number corresponding to their execution order from 01 to 16, except the directory *bin/scripts/pipeline_input_decostar/compute_scaffolding_adjacencies* containing scripts to produce scaffolding adjacencies with paired sequencing data and the scaffolding tool BESST. The scripts to scaffold with BESST have been developed for a computing cluster based on Sun Grid Engine (SGE), they have to be executed for the 18 Anopheles species and before the step 14 of the directory *bin/scripts/pipeline_input_decostar/compute_scaffolding_adjacencies*.-->-->

<!--For each script the content is headed by a description explaining its input, output and goal. The pipeline to produce input data for DeCoSTAR is divided in 16 steps (+scripts to compute scaffolding adjacencies), among which step 11 has been developed to be run on the computing cluster with SGE scripts and have to be adapted to run on other cluster architecture.-->

<!--To reproduce the results of our paper, the scripts have been written without need to give parameters to reproduce the architecture of this repository and global scripts have been produced to chained several scripts in one executive file. For more details on these scripts present in directory *bin/scripts/pipeline_input_decostar* see README file in this directory.-->



# References
[1] Duchemin, W., Anselmetti, Y., Patterson, M., Ponty, Y., Bérard, S., Chauve, C., Scornavacca, C., Daubin, V., Tannier, E. (2017). DeCoSTAR: Reconstructing the Ancestral Organization of Genes or Genomes Using Reconciled Phylogenies. Genome Biology and Evolution, 9(5), 1312–1319.

[2] Anselmetti, Y., Duchemin, W., Tannier, E., Chauve, C., Bérard S. (to appear). Phylogenetic signal from rearrangements in 18 {Anopheles} species by joint scaffolding extant and ancestral genomes. BMC Genomics.
