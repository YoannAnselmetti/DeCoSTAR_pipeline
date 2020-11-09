# configfile: "config_files/snakemake/config_18Anopheles_Xtopo.yaml"
# configfile: "config_files/snakemake/config_18Anopheles_WGtopo.yaml"
# configfile: "config_files/snakemake/config_21Anopheles_Xtopo.yaml"
configfile: "config_files/snakemake/config_27avian.yaml"
# configfile: "config_files/snakemake/config_12Drosophila.yaml"


# TO IMPROVE => GET "outputdir_decostar" 
# AND "prefix_decostar" FROM DECOSTAR PARAMETER FILE
# config["outputdir_decostar"]
# config["prefix_decostar"]

rule all:
	input:
		config["outputdir"]+"/results/ADJ_graph/log_CTGmap.log",
		config["outputdir"]+"/results/ADJ_graph/log_GENEmap.log",

#########
### ADJ graph with CTG node
#########
rule create_adj_CTGgraph_chrMAP:
	input:
		AGPdir=config["outputdir"]+"/"+config["AGP_dir"],
		CTG=config["outputdir"]+"/data/data_DeCoSTAR/CTG_file",
	output:
		dot=config["outputdir"]+"/results/ADJ_graph/DOT/CTG",
		svg=config["outputdir"]+"/results/ADJ_graph/SVG/CTG",
	shell:
		"echo -e \"\tscript:\n\";\
		for species in $(ls {input.AGPdir}); do \
			for AGPfile in $(ls {input.AGPdir}/$species); do \
				echo ./bin/scripts/post_decostar/ADJ_graph/create_CTGgraph_GM.py \
					{input.AGPdir}/$species/$AGPfile \
					{input.CTG} \
					$species \
					{output.dot}/$species/$species\_CTGmap.dot \
					{output.svg}/$species/$species\_CTGmap.svg; \
				./bin/scripts/post_decostar/ADJ_graph/create_CTGgraph_GM.py \
					{input.AGPdir}/$species/$AGPfile \
					{input.CTG} \
					$species \
					{output.dot}/$species/$species\_CTGmap.dot \
					{output.svg}/$species/$species\_CTGmap.svg; \
			done \
		done"

rule add_newADJ_to_adj_CTGgraph:
	input:
		dot=config["outputdir"]+"/results/ADJ_graph/DOT/CTG",
		svg=config["outputdir"]+"/results/ADJ_graph/SVG/CTG",
		new_adj_file_with_scaff=config["outputdir"]+"/"+
			config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+
			config["linearization_threshold"]+"_"+config["linearization_algo"]+
			"_new_extant_adjacencies_with_scaff_adj",
	output:
		log=config["outputdir"]+"/results/ADJ_graph/log_CTGmap.log",
	shell:
		"echo -e \"\tscript:\n\";\
		for species in $(ls {input.dot}); do \
			echo ./bin/scripts/post_decostar/ADJ_graph/create_CTGgraph.py \
				{input.dot}/$species/$species\_CTGmap.dot \
				{input.new_adj_file_with_scaff} \
				$species \
				{input.dot}/$species/$species\_"+config["prefix_decostar"]+
				"_Lin"+config["linearization_threshold"]+"_"+
				config["linearization_algo"]+"_CTGmap.dot \
				{input.svg}/$species/$species\_"+config["prefix_decostar"]+
				"_Lin"+config["linearization_threshold"]+"_"+
				config["linearization_algo"]+"_CTGmap.svg \
				0.0 \
				y; \
			./bin/scripts/post_decostar/ADJ_graph/create_CTGgraph.py \
				{input.dot}/$species/$species\_CTGmap.dot \
				{input.new_adj_file_with_scaff} \
				$species \
				{input.dot}/$species/$species\_"+config["prefix_decostar"]+
				"_Lin"+config["linearization_threshold"]+"_"+
				config["linearization_algo"]+"_CTGmap.dot \
				{input.svg}/$species/$species\_"+config["prefix_decostar"]+
				"_Lin"+config["linearization_threshold"]+"_"+
				config["linearization_algo"]+"_CTGmap.svg \
				0.0 \
				y; \
		done > {output.log}"



#########
### ADJ graph with GENE node
#########
rule create_adj_GENEgraph_chrMAP:
	input:
		AGPdir=config["outputdir"]+"/"+config["AGP_dir"],
		CTG=config["outputdir"]+"/data/data_DeCoSTAR/CTG_file",
	output:
		dot=config["outputdir"]+"/results/ADJ_graph/DOT/GENE",
		svg=config["outputdir"]+"/results/ADJ_graph/SVG/GENE",
	shell:
		"echo -e \"\tscript:\n\";\
		for species in $(ls {input.AGPdir}); do \
			for AGPfile in $(ls {input.AGPdir}/$species); do \
				echo ./bin/scripts/post_decostar/ADJ_graph/create_GENEgraph_GM.py \
					{input.AGPdir}/$species/$AGPfile \
					{input.CTG} \
					$species \
					{output.dot}/$species/$species\_GENEmap.dot \
					{output.svg}/$species/$species\_GENEmap.svg \
					y; \
				./bin/scripts/post_decostar/ADJ_graph/create_GENEgraph_GM.py \
					{input.AGPdir}/$species/$AGPfile \
					{input.CTG} \
					$species \
					{output.dot}/$species/$species\_GENEmap.dot \
					{output.svg}/$species/$species\_GENEmap.svg \
					y; \
			done \
		done"

rule add_newADJ_to_adj_GENEgraph:
	input:
		dot=config["outputdir"]+"/results/ADJ_graph/DOT/GENE",
		svg=config["outputdir"]+"/results/ADJ_graph/SVG/GENE",
		CTG=config["outputdir"]+"/data/data_DeCoSTAR/CTG_file",
		new_adj_file_with_scaff=config["outputdir"]+"/"+
			config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+
			config["linearization_threshold"]+"_"+config["linearization_algo"]+
			"_new_extant_adjacencies_with_scaff_adj",
	output:
		log=config["outputdir"]+"/results/ADJ_graph/log_GENEmap.log",
	shell:
		"echo -e \"\tscript:\n\";\
		for species in $(ls {input.dot}); do \
			echo ./bin/scripts/post_decostar/ADJ_graph/create_GENEgraph.py \
				{input.dot}/$species/$species\_GENEmap.dot \
				{input.CTG} \
				{input.new_adj_file_with_scaff} \
				$species \
				{input.dot}/$species/$species\_"+config["prefix_decostar"]+
				"_Lin"+config["linearization_threshold"]+"_"+
				config["linearization_algo"]+"_GENEmap.dot \
				{input.svg}/$species/$species\_"+config["prefix_decostar"]+
				"_Lin"+config["linearization_threshold"]+"_"+
				config["linearization_algo"]+"_GENEmap.svg \
				0.0 \
				y; \
			./bin/scripts/post_decostar/ADJ_graph/create_GENEgraph.py \
				{input.dot}/$species/$species\_GENEmap.dot \
				{input.CTG} \
				{input.new_adj_file_with_scaff} \
				$species \
				{input.dot}/$species/$species\_"+config["prefix_decostar"]+
				"_Lin"+config["linearization_threshold"]+"_"+
				config["linearization_algo"]+"_GENEmap.dot \
				{input.svg}/$species/$species\_"+config["prefix_decostar"]+
				"_Lin"+config["linearization_threshold"]+"_"+
				config["linearization_algo"]+"_GENEmap.svg \
				0.0 \
				y; \
		done > {output.log}"
