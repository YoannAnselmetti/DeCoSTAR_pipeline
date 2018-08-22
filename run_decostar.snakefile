configfile: "config_files/snakemake/config_21Anopheles_Xtopo.yaml"
# configfile: "config_files/snakemake/config_18Anopheles_Xtopo.yaml"
# configfile: "config_files/snakemake/config_27avian.yaml"


# TO IMPROVE => GET "outputdir_decostar" AND "prefix_decostar" FROM DECOSTAR PARAMETER FILE
# config["outputdir_decostar"]
# config["prefix_decostar"]


# Create output directory
snakemake.utils.makedirs(config["outputdir"]+"/"+config["outputdir_decostar"])


rule all:
	input:
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+".adjacencies.txt",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_kept",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_new_extant_adjacencies",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_new_extant_adjacencies_with_scaff_adj",
		# config["outputdir"]+"/results/FASTA/SCAFF"


rule run_decostar:
	input:
		config_file=config["config_decostar"],
	output:
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+".adjacencies.txt",
	shell:
		"bin/software_libraries/DeCoSTAR/bin/DeCoSTAR parameter.file={input.config_file}"


rule run_linearization_assign_scj:
	input:
		adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+".adjacencies.txt",
	output:
		new_adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_kept",
	shell:
		"python2 bin/scripts/post_decostar/linearize_genomes.py {input.adj_file} ALL "+config["linearization_threshold"]+" "+config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+" M1 1.0 0.001 > "+config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_linearization.log",


rule create_new_ADJ_file:
	input:
		CTG=config["outputdir"]+"/data/data_DeCoSTAR/CTG_file",
		adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_kept",
	output:
		new_adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_new_extant_adjacencies",
	shell:
		"python2 bin/scripts/post_decostar/compute_scaffstats_and_newADJfile.py {input.CTG} {input.adj_file} "+config["SEP"]+" 50 {output.new_adj_file}",

rule add_scaff_adj_infos:
	input:
		new_adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_new_extant_adjacencies",
		scaff_file=config["outputdir"]+"/data/data_DeCoSTAR/scaff_BESST_DeCoSTAR",
	output:
		new_adj_file_with_scaff=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_new_extant_adjacencies_with_scaff_adj",
	shell:
		"python2 bin/scripts/post_decostar/add_scaffADJ_to_new_extant_ADJ.py {input.new_adj_file} {input.scaff_file} {output.new_adj_file_with_scaff}"


# rule produce_scaffolding_FASTA_files:
# 	input:
# 		new_adj=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_new_extant_adjacencies",
# 		input_SCAFF=config["outputdir"]+"/"+config["ASSEMBLY_dir"]
# 	output:
# 		config["outputdir"]+"/results/FASTA/SCAFF",
# 		output_SCAFF=config["outputdir"]+"/results/FASTA/SCAFF/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_"
# 	shell:
# 		"python2 bin/scripts/post_decostar/create_FASTA_from_new_adjacencies.py {input.new_adj} {input.input_SCAFF} {output.output_SCAFF}"

