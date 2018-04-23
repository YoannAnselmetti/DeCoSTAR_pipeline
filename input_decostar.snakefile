# configfile: "config_files/config_input_18Anopheles_Xtopo.yaml"
configfile: "config_files/config_input_18Anopheles_WGtopo.yaml"


# Create output directory
snakemake.utils.makedirs(config["outputdir"]+"/data/data_DeCoSTAR/decostar")


rule all:
	input:
		config["outputdir"]+"/data/data_DeCoSTAR/GENE_file",
		config["outputdir"]+"/data/data_DeCoSTAR/CTG_file",
		config["outputdir"]+"/data/data_DeCoSTAR/scaff_BESST_ALL_3_TRIMMOMATIC3",
		config["outputdir"]+"/data/data_DeCoSTAR/scaff_BESST_DeCoSTAR",
		config["outputdir"]+"/data/data_DeCoSTAR/decostar/adjacencies.txt",
		config["outputdir"]+"/data/data_DeCoSTAR/decostar/adjacencies-scaff.txt",
		config["outputdir"]+"/data/data_DeCoSTAR/decostar/"+config["decostar_gene_trees"],
		config["outputdir"]+"/data/data_DeCoSTAR/decostar/distrib_"+config["decostar_gene_trees"]+".txt"



rule filter_GENE_with_geneTrees:
	input:
		input_GENE=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF",
		gene_trees=config["outputdir"]+"/"+config["gene_trees"]
	output:
		output_GENE=config["outputdir"]+"/data/data_DeCoSTAR/GENE_file"
	shell:
		"bin/scripts/pipeline_input_decostar/filter_GENE_with_geneTrees.py {input.input_GENE} {input.gene_trees} {output.output_GENE} "+config["SEP"]+" no"



rule create_CTG_file:
	input:
		GENE=config["outputdir"]+"/data/data_DeCoSTAR/GENE_file",
		SCAFF=config["outputdir"]+"/"+config["ASSEMBLY_dir"]
	output:
		CTG=config["outputdir"]+"/data/data_DeCoSTAR/CTG_file"
	shell:
		"bin/scripts/pipeline_input_decostar/create_CTG_file.py {input.GENE} {input.SCAFF} {output.CTG}"



rule create_scaff_adj_prefile:
	input:
		BESST_dir=config["outputdir"]+"/"+config["BESST_dir"]
	output:
		prescaff_file=config["outputdir"]+"/data/data_DeCoSTAR/scaff_BESST_ALL_3_TRIMMOMATIC3"
	shell:
		"bin/scripts/pipeline_input_decostar/create_scaff_adj_prefile.py {input.BESST_dir} 1000000000 3 {output.prescaff_file}"



rule create_scaff_adj_file_final:
	input:
		prescaff_file=config["outputdir"]+"/data/data_DeCoSTAR/scaff_BESST_ALL_3_TRIMMOMATIC3",
		CTG=config["outputdir"]+"/data/data_DeCoSTAR/CTG_file"
	output:
		scaff_file=config["outputdir"]+"/data/data_DeCoSTAR/scaff_BESST_DeCoSTAR"
	shell:
		"bin/scripts/pipeline_input_decostar/create_scaff_adj_file_final.py {input.prescaff_file} {input.CTG} {output.scaff_file}"



rule create_ADJfile_for_DeCoSTAR:
	input:
		GENE=config["outputdir"]+"/data/data_DeCoSTAR/GENE_file",
		scaff_file=config["outputdir"]+"/data/data_DeCoSTAR/scaff_BESST_DeCoSTAR"
	output:
		adj_file=config["outputdir"]+"/data/data_DeCoSTAR/decostar/adjacencies.txt",
		adj_noscaff_file=config["outputdir"]+"/data/data_DeCoSTAR/decostar/adjacencies-scaff.txt"
	shell:
		"bin/scripts/pipeline_input_decostar/create_ADJfile_for_DeCoSTAR.py {input.GENE} {input.scaff_file} {output.adj_file} {output.adj_noscaff_file} "+config["SEP"]



rule write_1_tree_per_file:
	input:
		gene_trees=config["outputdir"]+"/"+config["gene_trees"],
		GENE=config["outputdir"]+"/data/data_DeCoSTAR/GENE_file"
	output:
		gene_trees_dir=config["outputdir"]+"/data/data_DeCoSTAR/decostar/"+config["decostar_gene_trees"],
		distrib_gene_trees=config["outputdir"]+"/data/data_DeCoSTAR/decostar/distrib_"+config["decostar_gene_trees"]+".txt"
	shell:
		"bin/scripts/pipeline_input_decostar/write_1tree_per_file.py {input.gene_trees} {input.GENE} {output.gene_trees_dir} {output.gene_trees_dir} {output.distrib_gene_trees} "+config["SEP"]
