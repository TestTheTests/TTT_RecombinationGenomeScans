#!/bin/bash
##############################################################################
#
# File    : runAllelStats.sh
# History : 11/26/2018 - Created by Kevin Freeman (KF) 
#
##############################################################################
#
# This is a small wrapper for the python script calcAdditionalStatsFromVCF.py 
# It takes  a single positional argument indicating the simulation number and
# calls the python script on just that file. It always calls the script
# with a 1000bp window and using one core at a time.
#
#############################################################################

usage="\nrunAllelStats.sh [simNumber]\n\nex: runAllelStats.sh 10900\n"

if [ $# -ne 1 ]; then
	echo -e "\nrunAllelStats.sh takes exactly one argument. Usage: "
	echo -e $usage
	exit
fi

simNum=$1
datapath="../results_final/${simNum}_Invers_VCFallFILT.vcf.gz"
cmd="python3 ./calcAdditionalStatsFromVCF.py -d $datapath --ncpus 1"

echo -e "\nRunning command: $cmd\n"
$cmd
