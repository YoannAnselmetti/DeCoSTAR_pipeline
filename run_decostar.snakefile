configfile: "config_files/config_18Anopheles_Xtopo+scaff.yaml"


# Create output directory
snakemake.utils.makedirs(config["outputdir"]+"/"+config["outputdir_decostar"])


rule all:
	input:
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"adjacencies.txt",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"genes.txt",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"reconciliations.newick",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"species.txt",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"speciesTree.newick",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_orthogroups",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_disc",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_kept",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_new_extant_adjacencies",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_obs_scaffolds",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scaffolds",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scaffolds_assignment",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scaffolds_log",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scj",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scj_log",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scj_stats"



rule run_decostar:
	input:
		config_file=config["outputdir"]+"/"+config["config_decostar"]
	output:
		output_GENE=config["outputdir"]+"/"+config["outputdir_decostar"],
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"adjacencies.txt",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"genes.txt",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"reconciliations.newick",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"species.txt",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"speciesTree.newick"
	shell:
		"bin/software_libraries/DeCoSTAR/bin/DeCoSTAR parameter.file={input.config_file}"



rule run_linearization_assign_scj:
	input:
		outputDeCoSTAR=config["outputdir"]+"/"+config["outputdir_decostar"],
		geneWithGF=config["outputdir"]+"/data/GFF_to_GENE_files/with_filter/ALL_species_GENE_file_with_GF"
	output:
		output_GENE=config["outputdir"]+"/"+config["outputdir_decostar"],
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_orthogroups",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_disc",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_kept",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_new_extant_adjacencies",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_obs_scaffolds",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scaffolds",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scaffolds_assignment",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scaffolds_log",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scj",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scj_log",
		config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Linear_"+config["linearization_threshold"]+"_M1_scj_stats"
	shell:
		"bin/scripts/post_decostar/linearize_generate_stats_decostar.sh {input.outputDeCoSTAR} "+config["prefix_decostar"]+" "+config["linearization_threshold"]+" {input.geneWithGF}"
		