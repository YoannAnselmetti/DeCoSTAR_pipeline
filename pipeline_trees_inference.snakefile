configfile: "config_input_18Anopheles_Xtopo.yaml"


snakemake.utils.makedirs(config["outputdir"]+"/data/GENE_TREES")


rule all:
	input:
		config["outputdir"]+"/data/GFF_to_GENE_files/sorted_GFF",
		config["outputdir"]+"/data/GFF_to_GENE_files/GENE",
		config["outputdir"]+"/data/GFF_to_GENE_files/sorted_GENE",
		config["outputdir"]+"/data/GFF_to_GENE_files/filtered_GENE",
		config["outputdir"]+"/data/GFF_to_GENE_files/with_filter",
		config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file",
		config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_Inclusion_file",
		config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF",
		config["outputdir"]+"/data/GENE_TREES/unrooted_trees_filtered.nwk"



rule filter_GENE_with_gene_trees:
	input:
		input_GENE=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF"
        gene_trees=config["output_dir"]+"/"+config["gene_trees"]
	output:
		output_GENE=config["outputdir"]+"/data/GFF_to_GENE_files/sorted_GFF"
        
	log:
		config["outputdir"]+"/logs/GFF_to_GENE_files/sort_GFF.log"
	shell:
		"bin/scripts/pipeline_input_decostar/sort_GFF.sh {input.gff} {output.sorted_gff}"

