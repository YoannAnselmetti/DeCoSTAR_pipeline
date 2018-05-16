#!/bin/bash

scriptETE=./bin/scripts/post_decostar/compute_stats_graph/code/draw_tree_ETE.py
Xtree=data/INPUT_DATA/Anopheles_species_tree_X_topology.nwk
WGtree=data/INPUT_DATA/Anopheles_species_tree_WG_topology.nwk
XtopoPREF=results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1
XtopoRAWPREF=results/decostar/Xtopo_RAW/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1
XtopoWITHOUTSCAFFPREF=results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1
WGtopoPREF=results/decostar/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1



$scriptETE $Xtree $XtopoPREF\_0.1_M1_scj_stats $XtopoPREF\_0.1_M1_kept figures/ETE_species_trees/Xtopo+scaff/ETE_species_tree_Xtopo+scaff.pdf

$scriptETE $Xtree $XtopoRAWPREF\_0.1_M1_scj_stats $XtopoRAWPREF\_0.1_M1_kept figures/ETE_species_trees/Xtopo_RAW/ETE_species_tree_Xtopo_RAW.pdf

$scriptETE $Xtree $XtopoWITHOUTSCAFFPREF\_0.1_M1_scj_stats $XtopoWITHOUTSCAFFPREF\_0.1_M1_kept figures/ETE_species_trees/Xtopo-scaff/ETE_species_tree_Xtopo-scaff.pdf

$scriptETE $WGtree $WGtopoPREF\_0.1_M1_scj_stats $WGtopoPREF\_0.1_M1_kept figures/ETE_species_trees/WGtopo+scaff/ETE_species_tree_WGtopo+scaff.pdf
