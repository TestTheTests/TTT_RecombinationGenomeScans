#!/usr/bin/bash
###################################################################
# 
# File : run_rda_sims.sh
# 
###################################################################
#
# This bash script is used to run RDA on a set of VCFs representing
# a number of simulations, where for each simulation there is one
# VCF with all alleles included and one VCF where only quasi-indep
# alleles are included. It calls two R scripts: proc_sims_RDA.R, 
# which processes the VCF data and runs the RDA analysis, and 
# combine_RDA_tables.R, which organizes the files after all runs
# have finished so that each sim has one file that includes the
# results for both all alleles and quasi-indep only (PRUNED)
#
###################################################################

myDir="$(pwd)"             ## change this value if you are not running this script from same dir as data is stored
begin=10899                ## number of first sim to analyze
end=10900                  ## number of last sim to analyze
simType="Invers"           

SECONDS=0                  # used to time the analysis              

# make files for organization of all outputs if they don't already exist
mkdir -p ${myDir}/log_files
mkdir -p ${myDir}/plots
mkdir -p ${myDir}/results
mkdir -p ${myDir}/results/raw_data

for i in $(seq $begin $end)
do
	echo $i
	# make sure file exists
	if [ ! -f "${i}_Invers_VCFallFILT.vcf.gz" ]; then
	  echo "VCF for ${i} not found, skipping"
	  continue
	fi
	# process and analyze the pruned and unpruned files simultaneously
	Rscript proc_sims_RDA.R $i $simType $myDir "" >${myDir}"/log_files/${i}.log" &
        Rscript proc_sims_RDA.R $i $simType $myDir "PRUNED_" >${myDir}"/log_files/${i}_pruned.log" &
	wait # make sure processors don't get overloaded, could change this based on computer
done

# combine results
echo -e "\n\nCombining pruned and unpruned results\n"
Rscript combine_RDA_tables.R $begin $end $myDir $simType

# finished, report how long analysis took
echo -e "\n\nDone. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"
