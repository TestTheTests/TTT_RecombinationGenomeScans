#!/bin/bash
set -e
set -u
set -o pipefail

# This script only works if you are able to run ncpus in the background at once

mypath="/home/k.lotterhos/k_store/TestTheTests/TTT_RecombinationGenomeScans"
#mypath="/Users/katie/Desktop/TestTheTests/TTT_RecombinationGenomeScans"
cd $mypath
ncpus=50
start=10900 #10950
finish=$(($start + $ncpus-1))
echo $start $finish
simType="Invers"


# make a list of numbers to seed the sims 

##############
#### run slim sims
#############
SECONDS = 0 # used to time analyses
for i in $(seq $start $finish)
do	
	echo $i
	outfile=${i}"_Invers_out.txt"
	outpath=${mypath}"/results/"
	outpathfile=${outpath}${outfile}
	echo $outpathfile
	
	#run slim in background
	slim -d seed=$i -d mypath="'${outpath}'" simfiles/Recombination_Inversion_treeseq.slim > $outpathfile &
		# input my seed and path to outputput
done

wait ${!} #wait until the last background process is finished

echo -e "\n\nDone with SLiM sims. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"

##############
#### run python script to process tree sequence results
#############
SECONDS = 0 # used to time analyses
for i in $(seq $start $finish)
do
    echo $i
    python3 a_process_trees.py -s ${i} > ${i}"_Invers_pytree.out" 2> ${i}"_Invers_pytree.error" & echo $!
    sleep 10s
done

wait ${!} #wait until the last background process is finished

echo -e "\n\nDone with processing tree sequences. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"


##############
#### compress vcf files
#############

cd results
gzip -f *.vcf


##############
#### run R script
#############
#Rscript ../../Installs.r > Install_output.txt
#Rscript TO DO DELETE THIS LINE

SECONDS = 0 # used to time analyses
echo "Running R scripts"
for i in $(seq $start $finish)
do
    echo $i
    Rscript --vanilla ../src/b_Proc_Sims.R ${i} $simType > ${i}"_Invers_R.out" 2> ${i}"_Invers_R.error" & echo $!
    sleep 1m
done

wait ${!} #wait until the last background process is finished
echo -e "\n\nDone with processing first R script. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"

## Some of the R simulations failed the first time around, which I think happened because R was trying to read files that weren't completely written yet
## I added some wait times to the code, and all simulations worked the second time around


##############
#### run R script LEA
#############
SECONDS = 0 # used to time analyses
echo "Running R LEA"
for i in $(seq $start $finish)
do
    echo $i
    Rscript --vanilla ../src/c_Proc_Sims_LEA.R ${i} 'Invers' > ${i}"_Invers_R_LEA.out" 2> ${i}"_Invers_R_LEA.error" & echo $!
    sleep 10s
done

wait ${!} #wait until the last background process is finished
echo -e "\n\nDone with processing R script LEA. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"


##############
#### run R script RDA
#############
echo "Running R RDA"
for i in $(seq $start $finish)
do
    echo $i
    # make sure file exists
    if [ ! -f "${i}_Invers_VCFallFILT.vcf.gz" ]; then
        echo "VCF for ${i} not found, skipping"
        continue
    fi

    Rscript --vanilla ../src/d_proc_sims_RDA.R ${i} 'Invers' > ${i}"_Invers_R_RDA.out" 2> ${i}"_Invers_R_RDA.error" & echo $!
    sleep 10s
    # process and analyze the pruned and unpruned files simultaneously
    Rscript proc_sims_RDA.R $i $simType $myDir "" >${myDir}"/log_files/${i}.log" &
    Rscript proc_sims_RDA.R $i $simType $myDir "PRUNED_" >${myDir}"/log_files/${i}_pruned.log" &
    wait # make sure processors don't get overloaded, could change this based on computer
done

wait ${!} #wait until the last background process is finished
echo -e "\n\nDone with processing R script RDA. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"


##############
#### run scikit-allel
#############

##############
#### run baypass
#############
