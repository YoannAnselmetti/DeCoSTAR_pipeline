# configfile: "config_files/snakemake/config_18Anopheles_Xtopo.yaml"
configfile: "config_files/snakemake/config_27avian.yaml"
# configfile: "config_files/snakemake/config_9passeriformes.yaml"


# Recover parameter present in DeCoSTAR parameter file
# output.dir.decostar
# output.prefix


# Create output directory
snakemake.utils.makedirs(outputdir+"/"+output.dir.decostar)


rule all:
	input:
		config["outputdir"]+"/"+output.dir.decostar+"/"+output.prefix+"adjacencies.txt",
		# config["outputdir"]+"/"+output.dir.decostar+"/"+output.prefix+"_Linear_"+config["linearization_threshold"]+"_M1_new_extant_adjacencies",
		# config["outputdir"]+"/results/FASTA/SCAFF"


rule run_decostar:
	input:
		config_file=config["outputdir"]+"/"+config["config_decostar"]
	output:
		output_GENE=config["outputdir"]+"/"+output.dir.decostar,
		config["outputdir"]+"/"+output.dir.decostar+"/"+output.prefix+"adjacencies.txt"
	shell:
		"bin/software_libraries/DeCoSTAR/bin/DeCoSTAR parameter.file={input.config_file}"


# rule run_linearization_assign_scj:
# 	input:
# 		outputDeCoSTAR=config["outputdir"]+"/"+output.dir.decostar,
# 		geneWithGF=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF"
# 	output:
# 		output_GENE=config["outputdir"]+"/"+output.dir.decostar,
# 		config["outputdir"]+"/"+output.dir.decostar+"/"+output.prefix+"_Linear_"+config["linearization_threshold"]+"_M1_new_extant_adjacencies",
# 	shell:
# 		"bin/scripts/post_decostar/linearize_generate_stats_decostar.sh {input.outputDeCoSTAR} "+output.prefix+" "+config["linearization_threshold"]+" {input.geneWithGF}"



# rule produce_scaffolding_FASTA_files:
# 	input:
# 		new_adj=config["outputdir"]+"/"+output.dir.decostar+"/"+output.prefix+"_Linear_"+config["linearization_threshold"]+"_M1_new_extant_adjacencies",
# 		SCAFF=config["outputdir"]+"/"+config["ASSEMBLY_dir"]

# 	output:
# 		config["outputdir"]+"/results/FASTA/SCAFF"
# 	shell:
# 		"bin/scripts/post_decostar/generate_FASTA_files.py {input.new_adj} "


# rule library_file_gapFiller:
# 	input:
# 		lib_file=config["orientation_file"]
# 	output:
# 		lib_file=config["outputdir"]+"/results/libraries_file.txt"
# 	shell:
# 		"modif_lib_file.py {input.lib_file} {output.lib_file}"



# rule gap_filling:
	# Moussa




# rule adj_graph


# rule stats_graphics


# ...