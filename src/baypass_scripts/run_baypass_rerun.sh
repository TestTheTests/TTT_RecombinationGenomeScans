#!/usr/local/sh
######################################################################################
#
# File	  : run_baypass.pl
# History : 5/4/2018  Created by Kevin Freeman(KF)
#           5/11/2018 Added further logic
#	    9/14/2018 changed directory structure
#
#######################################################################################
#
# This script combines multiple scripts to form an automated pipeline for analyzing
# a batch of VCFs with Baypass. The script is meant to be run in a folder called
# run_baypass with a "results_final" directory in its parent directory that contains all the 
# necessary input files. The required scripts should be in the as the run_baypass directory. 
# It takes no arguments.
# 
#######################################################################################
#
# Example file tree with required scripts, input files, and directories:
# (Note --- converted_files, baypass_results, and log_files 
# directories and files are created during execution of the script)
# .
# |------results_final
# |       	|-------xxxx_Invers_indFILT.txt
# |			|------xxxx_Invers_ScanResults.txt
# |			|------xxxx_Invers_VCFallFILT.vcf.gz
# |			|------xxx_baypass_results.txt (created by script)
# |
# |-----src	
#	|---------------run_baypass     
# 			|---converted_files     (created by script)
#   			|	|------xxxx_Invers_VCFallFILT.covar
# 			|	|------xxxx_Invers_VCFallFILT.geno
# 			|	|------xxxx_Invers_VCFallFILT_PRUNED.covar
#			|	|------xxxx_Invers_VCFallFILT_PRUNED.geno
# 			|-----baypasss_results     (created by script)
# 			|	. (many files, see baypass 2.1 documnetation for more info)
# 			|	.
# 			|	. 
# 			|-----log_files            (created by script)
# 			|	|------xxxx_baypass_err.txt
# 			|	|------xxxx_baypass_log.txt					
# 			|		
# 	     		|-----run_baypass.sh
# 			|----vcf2baypass.pl
# 			|----pruneSNPs.pl
# 			|----getTableOfBaypassResults.pl					
#
#######################################################################################
SECONDS=0

# source ~/intel/bin/compilervars.sh intel64    # required to run the intel-compiled version of baypass (i_baypass)
mypath="/shared_lab/20181217_TTT_RecombinationGenomeScans/src/baypass_scripts"
ncpus=60
nsims=206
start=10900
finish=$(($start + $nsims-1))
echo $start $finish
rerun=(10911 10918 10924 10930 10942 10945 10957 10960 10963 10969 10980 10994 11001 11010 11045 11053 11059 11065 11084 11095 11097)
echo ${rerun[@]}
declare -A npopHash    # declare an associative array to store # of populations for each file
function call_baypass_no_mat { 
	cd $mypath
	file_name=$1
	out_suffix=$2
	i=$3
    
    # convert vcf to baypass format, store npop in associative array (hash)
    npopHash=( [${i}]="$(perl ./vcf2baypass.pl -vcf ../../results_final/${i}${file_name} \
    -pop ../../results_final/${i}_Invers_indFILT.txt \
    -colGroup 5 -colEnv 4 -colPheno 3 -colInFinal 6 \
    -outGeno ./converted_files/${i}_Invers_VCFallFILT${out_suffix} \
    -outCovar ./converted_files/${i}_Invers_VCFallFILT${out_suffix})")
    
    echo $npopHash
    # run baypass
    cmd="g_baypass -npop ${npopHash[${i}]} \
    -gfile ./converted_files/${i}_Invers_VCFallFILT${out_suffix}.geno \
    -efile ./converted_files/${i}_Invers_VCFallFILT${out_suffix}.covar -scalecov \
    -outprefix ./baypass_results/${i}${out_suffix} -nthreads $ncpus"

    $cmd >>"./log_files/${i}_baypass.log" 2>>"./log_files/${i}_baypass.err" || echo -e $cmd
}
function call_baypass_no_efile { 
	cd $mypath
	file_name=$1
	out_suffix=$2
	i=$3
    
    # convert vcf to baypass format, store npop in associative array (hash)
    npopHash=( [${i}]="$(perl ./vcf2baypass.pl -vcf ../../results_final/${i}${file_name} \
    -pop ../../results_final/${i}_Invers_indFILT.txt \
    -colGroup 5 -colEnv 4 -colPheno 3 -colInFinal 6 \
    -outGeno ./converted_files/${i}_Invers_VCFallFILT${out_suffix} \
    -outCovar ./converted_files/${i}_Invers_VCFallFILT${out_suffix})")
    
    echo $npopHash
    # run baypass
    cmd="g_baypass -npop ${npopHash[${i}]} \
    -gfile ./converted_files/${i}_Invers_VCFallFILT${out_suffix}.geno \
    -outprefix ./baypass_results/${i}${out_suffix} -nthreads $ncpus"

    $cmd >>"./log_files/${i}_baypass.log" 2>>"./log_files/${i}_baypass.err" || echo -e $cmd
}

cd $mypath
# create necessary directories for organization
mkdir -p ./log_files
mkdir -p ./converted_files
mkdir -p ./baypass_results

# run baypass on all the SNPs
echo -e "\n############## Running scripts on all snps ######################################"
for i in ${rerun[@]}
do
	pwd	
	echo -e "\n$i"
    	gunzip "../../results_final/${i}_Invers_VCFallFILT.vcf.gz"
    	if [ ! -f  "../../results_final/${i}_Invers_VCFallFILT.vcf" ]; then
    		echo "file not found"
    		continue
    	fi
    	call_baypass_no_mat "_Invers_VCFallFILT.vcf" "" $i

# prune the snps and run baypass again on the pruned files
	echo -e "\n############ Pruning SNP files, running scripts on pruned files #################"
	pwd
	echo -e "\n"$i
#	cd ../../results_final                                     # temporarily store pruned files in results_final
    	if [ ! -f  "../../results_final/${i}_Invers_ScanResults.txt" ]; then
    		echo "Scan results for ${i} not found"
    		continue
    	fi 
	perl ${mypath}/pruneSNPs.pl -scan ../../results_final/${i}_Invers_ScanResults.txt \
	 -vcf ../../results_final/${i}_Invers_VCFallFILT.vcf.gz
	
	gunzip ../../results_final/${i}_Invers_VCFallFILT.pruned.vcf.gz	
	call_baypass_no_efile "_Invers_VCFallFILT.pruned.vcf" "_PRUNED" "${i}"

# run baypass again, this time using the pruned covar matrix from the last step
	
    cd $mypath
    gfileName="_Invers_VCFallFILT.geno"
    efileName="_Invers_VCFallFILT.covar"
    omegaFile="_PRUNED_mat_omega.out"

    echo -e "\n########## Running baypass on all SNPs with pruned matrix #########################"
    echo -e "\n"$i
    if [ ! -f  ./converted_files/${i}$gfileName ]; then
    	echo "Geno file for ${i} not found"
    	continue
    fi 
    if [ ! -f  ./converted_files/${i}$efileName ]; then
    	echo "env file for ${i} not found"
    	continue
    fi
    if [ ! -f  ./baypass_results/${i}$omegaFile ]; then
    	echo "Omega matrix file for ${i} not found"
    	continue
    fi

    g_baypass -npop ${npopHash[$i]} -gfile ./converted_files/${i}$gfileName \
        -efile ./converted_files/${i}$efileName -scalecov \
        -omegafile ./baypass_results/${i}$omegaFile \
        -outprefix "./baypass_results/${i}_ALL_PRUNED_MAT" -nthreads $ncpus \
        >>./log_files/${i}_baypass_log.txt 2>>./log_files/${i}_baypass_err.txt

	gzip ../../results_final/${i}_Invers_VCFallFILT.vcf # re-zip for space

done
# pull out analysis results and store in a table
echo -e "\n################ Putting Results into Tables ##############################"
cd ${mypath}/baypass_results
pwd
perl ../getTableOfBaypassResults.pl -start $start -finish $finish -in ${mypath}/../../results_final -out ${mypath}/../../results
echo -e "\n\nDone. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"
exit 0
