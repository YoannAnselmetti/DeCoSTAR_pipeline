configfile: "config_files/config_input_18Anopheles_Xtopo.yaml"


snakemake.utils.makedirs(config["outputdir"]+"/logs/GFF_to_GENE_files")
snakemake.utils.makedirs(config["outputdir"]+"/data/GENE_TREES")


rule all:
	input:




rule filter_GENE_with_geneTrees:
	input:
		input_GENE=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF"
        gene_trees=config["output_dir"]+"/"+config["gene_trees"]
        sep=config["separator"]
	output:
		output_GENE=config["outputdir"]+"/data_DeCoSTAR/GENE_file"
	shell:
		"bin/scripts/pipeline_input_decostar/filter_GENE_with_geneTrees.py {input.input_GENE} {input.gene_trees} {output.output_GENE} {input.separator} no"


rule create_CTG_file:
	input:
		GENE=config["outputdir"]+"/data_DeCoSTAR/GENE_file"
        SCAFF=config["outputdir"]+"/"+config["SCAFF_dir"]
	output:
		CTG=config["outputdir"]+"/data_DeCoSTAR/CTG_file"
	shell:
		"bin/scripts/pipeline_input_decostar/create_CTG_file.py {input.GENE} {input.SCAFF} {output.CTG}"



rule create_scaff_adj_prefile:



rule create_scaff_adj_file_final:



rule create_ADJfile_for_DeCoSTAR:





rule write_1_tree_per_file:
	input:
		gene_trees=config["outputdir"]+"/"+config["gene_trees"]
        SCAFF=config["outputdir"]+"/"+config["SCAFF_dir"]
	output:
		CTG=config["outputdir"]+"/data_DeCoSTAR/CTG_file"
	shell:
		"bin/scripts/pipeline_input_decostar/create_CTG_file.py {input.GENE} {input.SCAFF} {output.CTG}"
