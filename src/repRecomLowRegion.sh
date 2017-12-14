#!/bin/bash
set -e
set -u
set -o pipefail

# This script only works if you are able to run ncpus in the background at once

#mypath="/shared_lab/TestTheTests/TTT_RecombinationGenomeScans"
mypath="/Users/katie/Desktop/TestTheTests/TTT_RecombinationGenomeScans"
cd $mypath
ncpus=2
start=9900
finish=$(($start + $ncpus-1))
echo $start $finish

# make a list of numbers to seed the sims 

##############
#### run slim sims
#############
for i in $(seq $start $finish)
do	
	echo $i
	outfile=${i}"_RecomLowReg_out.txt"
	outpath=${mypath}"/results/"
	outpathfile=${outpath}${outfile}
	echo $outpathfile
	
	#run slim in background
	slim -d seed=$i -d mypath="'${outpath}'" simfiles/Recombination_LowRegions.slim > $outpathfile &
		# input my seed and path to outputput
done

wait ${!} #wait until the last background process is finished

##############
#### compress vcf files
#############

cd results
gzip *.vcf


##############
#### run H12 in background
#############
for i in $(seq $start $finish)
# read all .MS files
do
    myms=${i}"_RecomLowReg_outputMS.txt"
    echo $myms
	H-scan -m -i $myms > ${i}"_RecomLowReg_H.txt" &
done

wait ${!} #wait until the last background process is finished

# to do: go to shared lab and create TestTheTests folder
# clone github repo to the folder

##############
#### run R script
#############
Rscript Installs.r > Install_output.txt
Rscript 

for i in $(seq $start $finish)
# read all .MS files
do
    Rscript src/Proc_Sims.R ${i} --vanilla > ${i}"_RecomLowReg_Rout.txt" &
done


