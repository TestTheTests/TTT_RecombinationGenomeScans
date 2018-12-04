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

# make a list of numbers to seed the sims 

##############
#### run slim sims
#############
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



##############
#### run python script to process tree sequence results
#############
for i in $(seq $start $finish)
do
    echo $i

done

wait ${!} #wait until the last background process is finished


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

echo "Running R scripts"
for i in $(seq $start $finish)
do
    echo $i
    Rscript --vanilla ../src/b_Proc_Sims.R ${i} 'Invers' > ${i}"_Invers_R.out" 2> ${i}"_Invers_R.error" & echo $!
    sleep 1m
done

## Some of the R simulations failed the first time around, which I think happened because R was trying to read files that weren't completely written yet
## I added some wait times to the code, and all simulations worked the second time around


##############
#### run R script LEA
#############

echo "Running R LEA"
for i in $(seq $start $finish)
do
    echo $i
    Rscript --vanilla ../src/c_Proc_Sims_LEA.R ${i} 'Invers' > ${i}"_Invers_R_LEAonlyR.out" 2> ${i}"_Invers_R_LEAonlyR.error" & echo $!
    sleep 10s
done

##############
#### run R script RDA
#############

##############
#### run scikit-allel
#############

##############
#### run baypass
#############
