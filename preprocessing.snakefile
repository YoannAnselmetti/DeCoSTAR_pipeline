configfile: "config_files/config_18Anopheles_Xtopo+scaff.yaml"


# Create output directories
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



rule sort_GFF:
	input:
		gff=config["outputdir"]+"/"+config["GFF_dir"]
	output:
		sorted_gff=config["outputdir"]+"/data/GFF_to_GENE_files/sorted_GFF"
	# log:
	# 	config["outputdir"]+"/logs/GFF_to_GENE_files/sort_GFF.log"
	shell:
		"bin/scripts/pipeline_input_decostar/sort_GFF.sh {input.gff} {output.sorted_gff}"



rule from_GFF_to_GENE:
	input:
		sorted_gff=config["outputdir"]+"/data/GFF_to_GENE_files/sorted_GFF"
	output:
		gene=config["outputdir"]+"/data/GFF_to_GENE_files/GENE",
		graph_gff=config["outputdir"]+"/data/GFF_to_GENE_files/GRAPH_GFF"
	# log:
	# 	config["outputdir"]+"/logs/GFF_to_GENE_files/from_GFF_to_GENE.log"
	shell:
		"bin/scripts/pipeline_input_decostar/from_GFF_to_GENE.py {input.sorted_gff} {output.gene} {output.graph_gff}"



rule sort_GENE:
	input:
		gene=config["outputdir"]+"/data/GFF_to_GENE_files/GENE"
	output:
		sorted_gene=config["outputdir"]+"/data/GFF_to_GENE_files/sorted_GENE"
	# log:
	# 	config["outputdir"]+"/logs/GFF_to_GENE_files/sort_GENE.log"
	shell:
		"bin/scripts/pipeline_input_decostar/sort_GENE.sh {input} {output}"



rule filter_GENE_with_families:
	input:
		families=config["outputdir"]+"/"+config["families"],
		sorted_gene=config["outputdir"]+"/data/GFF_to_GENE_files/sorted_GENE"
	output:
		filtered_gene=config["outputdir"]+"/data/GFF_to_GENE_files/filtered_GENE"
	# log:
	# 	config["outputdir"]+"/logs/GFF_to_GENE_files/filter_GENE_with_families.log"
	shell:
		"bin/scripts/pipeline_input_decostar/filter_GENE_with_families.py {input.families} {input.sorted_gene} {output.filtered_gene} "+config["SEP"]



rule detect_includedGenes:
	input: 
		gene=config["outputdir"]+"/data/GFF_to_GENE_files/filtered_GENE"
	output:
		overlapDir=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter",
		ouput_gene=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file",
		inclusion=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_Inclusion_file"
	# log:
	# 	config["outputdir"]+"/logs/GFF_to_GENE_files/detect_includedGenes.log"
	shell:
		"bin/scripts/pipeline_input_decostar/detect_includedGenes.py {input.gene} {output.overlapDir}"



rule add_geneFamilyID:
	input: 
		gene=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file",
		families=config["outputdir"]+"/"+config["families"]
	output:
		geneWithGF=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF"
	# log:
	# 	config["outputdir"]+"/logs/GFF_to_GENE_files/add_geneFamilyID.log"
	shell:
		"bin/scripts/pipeline_input_decostar/add_geneFamilyID.py {input.gene} {input.families} {output.geneWithGF} "+config["SEP"]



# THIS STEP HAS TO BE IMPROVED TO TAKE INTO ACCOUNT GENE FAMILIES FILE
rule filter_geneTrees_with_includedGenes:
	input: 
		families=config["outputdir"]+"/"+config["families"],
		inclusion=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_Inclusion_file"
	output:
		filtered_families=config["outputdir"]+"/data/GENE_TREES/unrooted_trees_filtered.nwk"
	# log:
	# 	config["outputdir"]+"/logs/GFF_to_GENE_files/add_geneFamilyID.log"
	shell:
		"bin/scripts/pipeline_input_decostar/filter_geneTrees_with_includedGenes.sh {input.families} {input.inclusion} {output.filtered_families}"
