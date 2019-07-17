configfile: "config_files/snakemake/config_18Anopheles_Xtopo.yaml"
# configfile: "config_files/snakemake/config_18Anopheles_WGtopo.yaml"
# configfile: "config_files/snakemake/config_21Anopheles_Xtopo.yaml"
# configfile: "config_files/snakemake/config_27avian.yaml"
# configfile: "config_files/snakemake/config_12Drosophila.yaml"


# TO IMPROVE => GET "outputdir_decostar" AND "prefix_decostar" FROM DECOSTAR PARAMETER FILE
# config["outputdir_decostar"]
# config["prefix_decostar"]

# Create output directory
snakemake.utils.makedirs(config["outputdir"]+"/"+config["outputdir_decostar"])

rule all:
	input:
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+".adjacencies.txt",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_kept",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_disc",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_scaffolding_stats",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_new_extant_adjacencies",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_new_extant_adjacencies_with_scaff_adj",
		config["outputdir"]+"/results/AGP",
		# config["outputdir"]+"/results/FASTA/SCAFF"


rule run_decostar:
	input:
		config_file=config["config_decostar"],
	output:
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+".adjacencies.txt",
	shell:
		"echo -e \"\tscript: bin/software_libraries/DeCoSTAR/bin/DeCoSTAR parameter.file={input.config_file}\";\
		bin/software_libraries/DeCoSTAR/bin/DeCoSTAR parameter.file={input.config_file}"


rule run_linearization:
	input:
		adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+".adjacencies.txt",
	output:
		kept_adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_kept",
		disc_adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_disc",
	shell:
		"echo -e \"\tscript: python2 bin/scripts/post_decostar/linearize_genomes.py {input.adj_file} ALL "+config["linearization_threshold"]+" "+config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+" "+config["linearization_algo"]+" 1.0 0.001 > "+config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_linearization.log\n\";\
		python2 bin/scripts/post_decostar/linearize_genomes.py {input.adj_file} ALL "+config["linearization_threshold"]+" "+config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+" "+config["linearization_algo"]+" 1.0 0.001 > "+config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_linearization.log"


rule create_new_ADJ_file:
	input:
		FASTA=config["outputdir"]+"/"+config["ASSEMBLY_dir"],
		CTG=config["outputdir"]+"/data/data_DeCoSTAR/CTG_file",
		kept_adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_kept",
		disc_adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_disc",
	output:
		new_adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_new_extant_adjacencies",
		scaff_stats=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_scaffolding_stats",
	shell:
		"echo -e \"\tscript: python2 bin/scripts/post_decostar/compute_scaffstats_and_newADJfile.py {input.FASTA} {input.CTG} {input.kept_adj_file} {input.disc_adj_file} "+config["SEP"]+" {output.new_adj_file} {output.scaff_stats}\";\
		python2 bin/scripts/post_decostar/compute_scaffstats_and_newADJfile.py {input.FASTA} {input.CTG} {input.kept_adj_file} {input.disc_adj_file} "+config["SEP"]+" {output.new_adj_file} {output.scaff_stats}"


rule add_scaff_adj_infos:
	input:
		new_adj_file=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_new_extant_adjacencies",
		scaff_file=config["outputdir"]+"/data/data_DeCoSTAR/scaff_BESST_DeCoSTAR",
	output:
		new_adj_file_with_scaff=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_new_extant_adjacencies_with_scaff_adj",
	shell:
		"echo -e \"\tscript: python2 bin/scripts/post_decostar/add_scaffADJ_to_new_extant_ADJ.py {input.new_adj_file} {input.scaff_file} {output.new_adj_file_with_scaff}\";\
		python2 bin/scripts/post_decostar/add_scaffADJ_to_new_extant_ADJ.py {input.new_adj_file} {input.scaff_file} {output.new_adj_file_with_scaff}"


rule produce_AGP_files:
	input:
		new_adj_file_with_scaff=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_new_extant_adjacencies_with_scaff_adj",
		input_SCAFF=config["outputdir"]+"/"+config["ASSEMBLY_dir"]
	output:
		output_AGP=config["outputdir"]+"/results/AGP"
	shell:
		"echo -e \"\tscript: python2 bin/scripts/post_decostar/create_AGP_from_new_adjacencies.py {input.new_adj_file_with_scaff} {input.input_SCAFF} "+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_ {output.output_AGP}\";\
		python2 bin/scripts/post_decostar/create_AGP_from_new_adjacencies.py {input.new_adj_file_with_scaff} {input.input_SCAFF} "+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_ {output.output_AGP}"


rule produce_newFASTA_files:
	input:
		AGP_dir=config["outputdir"]+"/results/AGP",
		input_SCAFF=config["outputdir"]+"/"+config["ASSEMBLY_dir"]
	output:
		output_SCAFF=config["outputdir"]+"/results/FASTA/SCAFF",
	shell:
		"echo -e \"\tscript: python2 bin/scripts/post_decostar/create_newFASTA_from_AGP.py {input.AGP_dir} "+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_ {input.input_SCAFF} {output.output_SCAFF}\";\
		python2 bin/scripts/post_decostar/create_newFASTA_from_AGP.py {input.AGP_dir} "+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_"+config["linearization_algo"]+"_ {input.input_SCAFF} {output.output_SCAFF}"
