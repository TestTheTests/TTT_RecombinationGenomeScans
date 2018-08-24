#!/usr/bin/bash
###################################################################
# 
# File   : run_hapflk.sh
# History: created by Kevin Freeman June 2018
#
###################################################################
#
# run_hapflk.sh is a script that runs hapflk analysis on a set
# of VCFs. It should be run in a folder containing the 1) VCFs 
# 2) "scan results files" indicating independence (for pruning)
# and 3) a "src" subfolder that contains the following scripts:
#   - vcf2plink.pl
#   - splitPedMapbyChrom.pl
#   - pruneSNPs.pl
#   - runHapflkAllChroms.sh 
#   - concat_chr.R
#
###################################################################

myDir="$(pwd)"             ## change this value if you are not running this script from same dir as data is stored
begin=10901                ## number of first sim to analyze
end=10961                  ## number of last sim to analyze
simType="Invers"           
K=15			   # number of hapflk clusters

SECONDS=0                  # used to time the analysis 

cd $myDir
mkdir -p plink_files
mkdir -p results
mkdir -p results/ALL
mkdir -p results/PRUNED
mkdir -p results/ALL_PRUNED_KIN
mkdir -p kinship
mkdir -p logfiles

for i in $(seq $begin $end) 
do
	mkdir -p $myDir/plink_files/$i
	mkdir -p $myDir/results/ALL/$i
	mkdir -p $myDir/results/PRUNED/$i
	mkdir -p $myDir/results/ALL_PRUNED_KIN/$i
	
	echo $i
	curVCF="$myDir/${i}_${simType}_VCFallFILT.vcf"
	curPruned="$myDir/${i}_PRUNED_${simType}_VCFallFILT.vcf"
	curPop="$myDir/${i}_${simType}_indFILT.txt"
	# check if gunzip is needed and if vcf is there
	if [ ! -f $curVCF ]; then
		gunzip ${curVCF}.gz
	fi
	if [  ! -f $curVCF ]; then
		echo "could not find vcf for $i, skipping"
		continue
	fi
	# create pruned file
	${myDir}/src/pruneSNPs.pl -vcf $curVCF -scan $myDir/${i}_${simType}_"ScanResults.txt" -outfile $curPruned
	
	echo -e "\nconverting to PLINK format\n"
	
	allPlink="${myDir}/plink_files/${i}/${i}_${simType}_VCFallFILT"
	prunedPlink="${myDir}/plink_files/${i}/${i}_PRUNED_${simType}_VCFallFILT"
	numGroups=$(perl ${myDir}/src/vcf2plink.pl -vcf $curVCF -pop $curPop -outfile $allPlink 2>>logfiles/$i.log &)
	perl ${myDir}/src/vcf2plink.pl -vcf $curPruned -pop $curPop -outfile $prunedPlink 2>>logfiles/$i.log &
	 wait
	
	# generate a kinship file using all data
	gzip $curVCF
	gzip $curPruned
	hapflk --file $allPlink -p "$myDir/kinship/$i" >>logfiles/$i.log 2>&1 &
	kinshipAll="$myDir/kinship/${i}_fij.txt"
	
	hapflk --file $prunedPlink -p "${myDir}/kinship/${i}_PRUNED" >>logfiles/${i}.log 2>&1 & 
	wait 
	kinshipPruned="$myDir/kinship/${i}_PRUNED_fij.txt"
	
	# split files by chromosome
	perl ${myDir}/src/splitPedMapbyChrom.pl -infile $allPlink -outfile "${myDir}/plink_files/${i}/${i}" >>logfiles/$i.log 2>&1 &
	perl ${myDir}/src/splitPedMapbyChrom.pl -infile $prunedPlink -outfile "${myDir}/plink_files/${i}/${i}_PRUNED" >>logfiles/$i.log 2>&1 &
	wait
	
	################# Main hafplk analysis #############################
	### unpruned
	bash ${myDir}/src/runHapflkAllChroms.sh ${myDir}/plink_files/${i} ${myDir}/results/ALL/$i $i $kinshipAll $K >>logfiles/$i.log 2>&1
	# combine the results from the 10 chromosomes using an Rscript
	Rscript ${myDir}/src/concat_chr.R "${myDir}/results/ALL/$i/${i}_chr" "${myDir}/results/ALL/" $i
	### unpruned with pruned kinship
	bash ${myDir}/src/runHapflkAllChroms.sh ${myDir}/plink_files/${i} ${myDir}/results/ALL_PRUNED_KIN/$i ${i} $kinshipPruned $K >>logfiles/$i.log 2>&1
	Rscript ${myDir}/src/concat_chr.R "${myDir}/results/ALL_PRUNED_KIN/$i/${i}_chr" "${myDir}/results/ALL_PRUNED_KIN/" $i
	#### pruned with pruned kinship
	bash ${myDir}/src/runHapflkAllChroms.sh ${myDir}/plink_files/${i} ${myDir}/results/PRUNED/$i ${i}_PRUNED $kinshipPruned $K >>logfiles/$i.log 2>&1
	Rscript ${myDir}/src/concat_chr.R "${myDir}/results/PRUNED/$i/${i}_PRUNED_chr" "${myDir}/results/PRUNED/" $i
	echo -e "\n\nFinished $i. $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min into analysis"

done

echo -e "\n\nDone. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"
