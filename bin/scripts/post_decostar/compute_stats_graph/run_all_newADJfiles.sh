#!/bin/bash

scriptNewADJ=./bin/scripts/post_decostar/compute_stats_graph/code/compute_scaffstats_and_newADJfile.py
CTGfile=data/data_DeCoSTAR/CTG_file
Nxx=50

Xtopokept=results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_0.1_M1_kept
XtopoRAWkept=results/decostar/Xtopo_RAW/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1_0.1_M1_kept
WGtopokept=results/decostar/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1_0.1_M1_kept

XtoponewADJ=results/decostar/Xtopo+scaff/DeCoSTAR_Anopheles_Xtopo+scaff_Boltz_0.1_0.1_M1_new_extant_adjacencies
XtopoRAWnewADJ=results/decostar/Xtopo_RAW/DeCoSTAR_Anopheles_Xtopo_RAW_Boltz_0.1_0.1_M1_new_extant_adjacencies
WGtoponewADJ=results/decostar/WGtopo+scaff/DeCoSTAR_Anopheles_WGtopo+scaff_Boltz_0.1_0.1_M1_new_extant_adjacencies


$scriptNewADJ $CTGfile $Xtopokept $Nxx $XtoponewADJ
$scriptNewADJ $CTGfile $XtopoRAWkept $Nxx $XtopoRAWnewADJ
$scriptNewADJ $CTGfile $WGtopokept $Nxx $WGtoponewADJ