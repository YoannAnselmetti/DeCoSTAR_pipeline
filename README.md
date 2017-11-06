Pipeline for a friendly use of the DeCoSTAR software from standard format files
=====

Pipeline to produce input data for DeCoSTAR software to apply ARt-DeCo, ADseq and DeClone algorithms on a dataset composed of files to the standard format described in the part XXXX.


### Requirements

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

### Quickstart ###

```
install_DeCoSTAR.sh
snakemake --snakefile preprocessing.snakefile -j <N>


```

## Pipeline
The pipeline to execute DeCoSTAR on a dataset is divided in 4 main parts:
1) Subpipeline to produce input data for the DeCoSTAR software to use ARt-DeCo or ADseq, and DeClone. 
2) 
3)
4) 


** The directory *bin/scripts/pipeline_input_decostar* contains scripts developed to process raw input data of the 18 Anopheles dataset (located in directory *data/INPUT_DATA*) into input data for the software DeCoSTAR to execute the ARt-DeCo and ADseq algorithms. All the data produced by these scripts are located in the directory *data/*. In this directory, the script names are prefixed by a number corresponding to their execution order from 01 to 16, except the directory *bin/scripts/pipeline_input_decostar/compute_scaffolding_adjacencies* containing scripts to produce scaffolding adjacencies with paired sequencing data and the scaffolding tool BESST. The scripts to scaffold with BESST have been developed for a computing cluster based on Sun Grid Engine (SGE), they have to be executed for the 18 Anopheles species and before the step 14 of the directory *bin/scripts/pipeline_input_decostar/compute_scaffolding_adjacencies*.

For each script the content is headed by a description explaining its input, output and goal. The pipeline to produce input data for DeCoSTAR is divided in 16 steps (+scripts to compute scaffolding adjacencies), among which step 11 has been developed to be run on the computing cluster with SGE scripts and have to be adapted to run on other cluster architecture.

To reproduce the results of our paper, the scripts have been written without need to give parameters to reproduce the architecture of this repository and global scripts have been produced to chained several scripts in one executive file. For more details on these scripts present in directory *bin/scripts/pipeline_input_decostar* see README file in this directory.




### INPUT
* species tree file in newick format
* species and expected chromosome number file 
* gene trees/families file in newick or NHX format
* 
*








### Preprocessing






### References
[1] Duchemin, W., Anselmetti, Y., Patterson, M., Ponty, Y., Bérard, S., Chauve, C., Scornavacca, C., Daubin, V., Tannier, E. (2017). DeCoSTAR: Reconstructing the Ancestral Organization of Genes or Genomes Using Reconciled Phylogenies. Genome Biology and Evolution, 9(5), 1312–1319.

[2] Anselmetti, Y., Duchemin, W., Tannier, E., Chauve, C., Bérard S. (to appear). Phylogenetic signal from rearrangements in 18 {Anopheles} species by joint scaffolding extant and ancestral genomes. BMC Genomics.