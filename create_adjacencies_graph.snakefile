# configfile: "config_files/snakemake/config_27avian.yaml"
configfile: "config_files/snakemake/config_21Anopheles_Xtopo.yaml"
# configfile: "config_files/snakemake/config_18Anopheles_Xtopo.yaml"


# TO IMPROVE => GET "outputdir_decostar" AND "prefix_decostar" FROM DECOSTAR PARAMETER FILE
# config["outputdir_decostar"]
# config["prefix_decostar"]


rule all:
	input:
		config["outputdir"]+"/results/ADJ_graph/",


#########
### ADJ graph with CTG node
#########
rule create_adj_CTGgraph_chrMAP:
	input:
		AGP_file=config["outputdir"]+"/data/INPUT_DATA/ADJ_graph/Zbor_scaff_anchored_in_Tgut_chr_050618_filt.agp",
		CTG=config["outputdir"]+"/data/data_DeCoSTAR/CTG_file",
		spe="Zosterops_borbonicus",
	output:
		config["outputdir"]+"/results/ADJ_graph/",
		dot=config["outputdir"]+"/results/ADJ_graph/DOT/CTG/Zosterops_borbonicus/Zbor_CTG_filt.dot",
		svg=config["outputdir"]+"/results/ADJ_graph/SVG/CTG/Zosterops_borbonicus/Zbor_CTG_filt.svg",
	shell:
		"python2 bin/scripts/post_decostar/ADJ_graph/create_CTGgraph_GM.py {input.AGP_file} {input.CTG_file} {input.spe} {ouput.dot} {ouput.svg}"


rule add_newADJ_to_adj_CTGgraph:
	input:
		chr_map_dot=config["outputdir"]+"/results/ADJ_graph/DOT/CTG/Zosterops_borbonicus/Zbor_CTG_filt.dot",
		new_adj_file_with_scaff=config["outputdir"]+"/"+config["outputdir_decostar"]+"/"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_new_extant_adjacencies_with_scaff_adj",
		spe="Zosterops_borbonicus",
	output:
		config["outputdir"]+"/results/ADJ_graph/",
		dot=config["outputdir"]+"/results/ADJ_graph/DOT/CTG/Zosterops_borbonicus/Zbor+"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_CTG_filt.dot",
		dot=config["outputdir"]+"/results/ADJ_graph/SVG/CTG/Zosterops_borbonicus/Zbor+"+config["prefix_decostar"]+"_Lin"+config["linearization_threshold"]+"_M1_CTG_filt.svg",
	shell:
		"python2 bin/scripts/post_decostar/ADJ_graph/create_CTGgraph_GM.py {input.chr_map_dot} {input.new_adj_file_with_scaff} {input.spe} {ouput.dot} {ouput.svg}"


#########
### ADJ graph with GENE node
#########


