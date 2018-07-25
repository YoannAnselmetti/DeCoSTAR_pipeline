#!/bin/bash

scriptSCATTERPLOT=./bin/scripts/post_decostar/compute_stats_graph/code/scatterplot_frag_on_scoreAdjDisc.py
XtopoPREF=results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1
XtopoRAWPREF=results/decostar/Xtopo_RAW/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1
WGtopoPREF=results/decostar/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1
CTGfile=data/data_DeCoSTAR/CTG_file
CHRfile=data/INPUT_DATA/18Anopheles_species


$scriptSCATTERPLOT $XtopoPREF\_0.1_M1_disc $XtopoPREF\_0.1_M1_scj_stats $CTGfile $CHRfile $XtopoPREF.species.txt figures/scatterplot/Xtopo/scatterplot_Xtopo

$scriptSCATTERPLOT $XtopoRAWPREF\_0.1_M1_disc $XtopoRAWPREF\_0.1_M1_scj_stats $CTGfile $CHRfile $XtopoRAWPREF.species.txt figures/scatterplot/Xtopo_RAW/scatterplot_Xtopo_RAW

$scriptSCATTERPLOT $WGtopoPREF\_0.1_M1_disc $WGtopoPREF\_0.1_M1_scj_stats $CTGfile $CHRfile $WGtopoPREF.species.txt figures/scatterplot/WGtopo/scatterplot_WGtopo