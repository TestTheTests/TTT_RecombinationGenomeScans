#!/bin/bash
set -e
set -u
set -o pipefail

source activate msprime-env

conda list

#mypath="/home/k.lotterhos/k_store/TestTheTests/TTT_RecombinationGenomeScans"
mypath="/Users/katie/Desktop/Repos/TestTheTests/TTT_RecombinationGenomeScans"
cd $mypath
nreps=210 #number of reps
start=10900  #10900 10950
finish=$(($start + $nreps-1))
echo -e $start $finish
simType="Invers"

# make a list of numbers to seed the sims 


for i in $(seq $start $finish)
do	
    echo -e "\n\n"
    echo $i
	outfile=${i}"_Invers_out.txt"
	outpath=${mypath}"/results/"
	outpathfile=${outpath}${outfile}
    #echo $outpathfile

    ##############
    #### run slim sims (takes 10 min)
    #############
    echo "running SLiM"
    SECONDS=0 # used to time analyses, no spaces around "=" sign
	#run slim in background
	slim -d seed=$i -d mypath="'${outpath}'" simfiles/Recombination_Inversion_treeseq.slim > $outpathfile 2> $outpathfile".error" & echo $!
		# input my seed and path to outputput
    wait
    echo "Done with SLiM sims. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"

    # make sure output from last analysis exists
    if [ ! -f "results/${i}_Invers.trees" ]; then
        echo "trees file for ${i} not found, skipping"
        continue # skip rest of code and go to next $i
    fi

    ##############
    #### run python script to process tree sequence results  (takes 1 min)
    #############
    cd results
    echo "processing tree sequences in python"
    SECONDS=0 # used to time analyses
    python3 ../src/a_process_trees.py -s ${i} > ${i}"_Invers_pytree.out" 2> ${i}"_Invers_pytree.error" & echo $!

    wait #wait until the last background process is finished

    echo "Done with processing tree sequences. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"


    ##############
    #### compress vcf files
    #############
    gzip -f *.vcf

    ##############
    #### run R script  (takes 3 min)
    #############
    SECONDS=0 # used to time analyses
    echo "Running R scripts"
    Rscript --vanilla ../src/b_Proc_Sims.R ${i} $simType > ${i}"_Invers_R.out" 2> ${i}"_Invers_R.error" & echo $!
    wait #wait until the last background process is finished
    echo "Done with processing first R script. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"

    ##############
    #### run R script RDA (takes 1 min)
    #############
    SECONDS=0 # used to time analyses
    echo "Running R RDA"

    # make sure output from last analysis exists
    if [ ! -f "../results_final/${i}_Invers_VCFallFILT.vcf.gz" ]; then
        echo "VCF for ${i} not found, skipping"
        continue
    fi

    Rscript --vanilla ../src/c_proc_sims_RDA.R ${i} 'Invers' > ${i}"_Invers_R_RDA.out" 2> ${i}"_Invers_R_RDA.error" & echo $!

    wait #wait until the last background process is finished
    echo "Done with processing R script RDA. Analysis took $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"

    ##############
    #### run scikit-allel takes about 3 min
    #############
    echo "Running scikit-allel"
    SECONDS=0
    #pwd
	datapath="../results_final/${i}_${simType}_VCFallFILT.vcf.gz"
	cmd="python3 ../src/d_calc_sk-allel_stats.py -d $datapath --ncpus 1"

	echo -e "\nRunning command: $cmd\n" > ${i}"_${simType}_sk-allel.out"
	$cmd > ${i}"_${simType}_sk-allel.out" 2> ${i}"_${simType}_sk-allel.error" & echo $!

    wait # wait until the last background process is finished
    echo "Done with processing scikit-allel. Analysis took $((SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min"

    ##############
    #### remove all temp files in the results folder
    #############
    echo "cleaning up files"
    echo $i
    mv ${i}"_genotypes.lfmm" ../results_final/
    mv ${i}"_gradients.env" ../results_final/
    mv *.error ../results_final/
    cd ../results/ | rm *
    cd ..
done


#cd results_final
#mv *.lfmm ../results
#mv *.env ../results

#transfer files in results_final to cluster
# run LEA and Baypass
# remove .lfmm and .env files
