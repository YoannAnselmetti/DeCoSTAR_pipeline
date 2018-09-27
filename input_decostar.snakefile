configfile: "config_files/snakemake/config_21Anopheles_Xtopo.yaml"
# configfile: "config_files/snakemake/config_18Anopheles_Xtopo.yaml"
# configfile: "config_files/snakemake/config_27avian.yaml"

# Create output directory
snakemake.utils.makedirs(config["outputdir"]+"/data/data_DeCoSTAR/decostar")

rule all:
	input:
		config["outputdir"]+"/data/data_DeCoSTAR/decostar/adjacencies.txt",
		config["outputdir"]+"/data/data_DeCoSTAR/decostar/distrib_"+config["decostar_gene_trees"]+".txt"


rule filter_GENE_with_geneTrees:
	input:
		input_GENE=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF",
		gene_trees=config["outputdir"]+"/"+config["gene_trees"]
	output:
		output_GENE=config["outputdir"]+"/data/data_DeCoSTAR/GENE_file"
	shell:
		"echo -e \"\tscript: python2 bin/scripts/pipeline_input_decostar/filter_GENE_with_geneTrees.py {input.input_GENE} {input.gene_trees} {output.output_GENE} "+config["SEP"]+" no\";\
		python2 bin/scripts/pipeline_input_decostar/filter_GENE_with_geneTrees.py {input.input_GENE} {input.gene_trees} {output.output_GENE} "+config["SEP"]


rule create_CTG_file:
	input:
		GENE=config["outputdir"]+"/data/data_DeCoSTAR/GENE_file",
		SCAFF=config["outputdir"]+"/"+config["ASSEMBLY_dir"]
	output:
		CTG=config["outputdir"]+"/data/data_DeCoSTAR/CTG_file"
	shell:
		"echo -e \"\tscript: python2 bin/scripts/pipeline_input_decostar/create_CTG_file.py {input.GENE} {input.SCAFF} {output.CTG}\";\
		python2 bin/scripts/pipeline_input_decostar/create_CTG_file.py {input.GENE} {input.SCAFF} {output.CTG}"


# CONFIGURED FOR THE SCAFFOLDING TOOL: "BESST"
rule create_scaff_adj_file:
	input:
		BESST_dir=config["outputdir"]+"/"+config["BESST_dir"],
		CTG=config["outputdir"]+"/data/data_DeCoSTAR/CTG_file"
	output:
		scaff_file=config["outputdir"]+"/data/data_DeCoSTAR/scaff_BESST_DeCoSTAR"
	shell:
		"echo -e \"\tscript: python2 bin/scripts/pipeline_input_decostar/create_scaff_adj_file.py {input.BESST_dir} {input.CTG} 1000000000 4 {output.scaff_file}\";\
		python2 bin/scripts/pipeline_input_decostar/create_scaff_adj_file.py {input.BESST_dir} {input.CTG} 1000000000 4 {output.scaff_file}"


rule create_ADJfile_for_DeCoSTAR:
	input:
		GENE=config["outputdir"]+"/data/data_DeCoSTAR/GENE_file",
		scaff_file=config["outputdir"]+"/data/data_DeCoSTAR/scaff_BESST_DeCoSTAR"
	output:
		adj_file=config["outputdir"]+"/data/data_DeCoSTAR/decostar/adjacencies.txt",
		adj_noscaff_file=config["outputdir"]+"/data/data_DeCoSTAR/decostar/adjacencies-scaff.txt"
	shell:
		"echo -e \"\tscript: python2 bin/scripts/pipeline_input_decostar/create_ADJfile_for_DeCoSTAR.py {input.GENE} {input.scaff_file} {output.adj_file} {output.adj_noscaff_file} "+config["SEP"]+"\";\
		python2 bin/scripts/pipeline_input_decostar/create_ADJfile_for_DeCoSTAR.py {input.GENE} {input.scaff_file} {output.adj_file} {output.adj_noscaff_file} "+config["SEP"]


if config["gene_trees"]:
	rule write_1_tree_per_file:
		input:
			gene_trees=config["outputdir"]+"/"+config["gene_trees"],
			GENE=config["outputdir"]+"/data/data_DeCoSTAR/GENE_file"
		output:
			gene_trees_dir=config["outputdir"]+"/data/data_DeCoSTAR/decostar/"+config["decostar_gene_trees"],
			distrib_gene_trees=config["outputdir"]+"/data/data_DeCoSTAR/decostar/distrib_"+config["decostar_gene_trees"]+".txt"
		shell:
			"echo -e \"\tscript: python2 bin/scripts/pipeline_input_decostar/write_1tree_per_file.py {input.gene_trees} {input.GENE} {output.gene_trees_dir} {output.gene_trees_dir} {output.distrib_gene_trees} "+config["SEP"]+"\";\
			python2 bin/scripts/pipeline_input_decostar/write_1tree_per_file.py {input.gene_trees} {input.GENE} {output.gene_trees_dir} {output.gene_trees_dir} {output.distrib_gene_trees} "+config["SEP"]
else:
	rule write_1_tree_per_file:
		input:
			gene_trees=config["outputdir"]+"/data/GENE_TREES/unrooted_trees_filtered.nwk",
			GENE=config["outputdir"]+"/data/data_DeCoSTAR/GENE_file"
		output:
			gene_trees_dir=config["outputdir"]+"/data/data_DeCoSTAR/decostar/"+config["decostar_gene_trees"],
			distrib_gene_trees=config["outputdir"]+"/data/data_DeCoSTAR/decostar/distrib_"+config["decostar_gene_trees"]+".txt"
		shell:
			"echo -e \"\tscript: python2 bin/scripts/pipeline_input_decostar/write_1tree_per_file.py {input.gene_trees} {input.GENE} {output.gene_trees_dir} {output.gene_trees_dir} {output.distrib_gene_trees} "+config["SEP"]+"\";\
			python2 bin/scripts/pipeline_input_decostar/write_1tree_per_file.py {input.gene_trees} {input.GENE} {output.gene_trees_dir} {output.gene_trees_dir} {output.distrib_gene_trees} "+config["SEP"]
