#!/bin/bash
set -e
set -u
set -o pipefail

# This script only works if you are able to run ncpus in the background at once

mypath="/home/k.lotterhos/k_store/TestTheTests/TTT_RecombinationGenomeScans"
#mypath="/Users/katie/Desktop/TestTheTests/TTT_RecombinationGenomeScans"
cd $mypath
ncpus=5
start=10900
nreps=25
each=$((nreps / ncpus))
finish=$(($start + $nreps-1))
starteach=$(seq $start $each $finish)
echo $start $each $finish
echo $starteach

# make a list of numbers to seed the sims

cd results

##############
#### run "each" number of reps of R script LEA in serial
#############

run_LEAreps(){
	startthis=$1
	endthis=$2
	for i in $(seq $startthis $endthis)
	do
   	 	echo $i
    	Rscript --vanilla ../src/Proc_Sims_LEA.R ${i} 'Invers' > ${i}"_Invers_R_LEAonlyR.out" 2> ${i}"_Invers_R_LEAonlyR.error" echo $!
    	sleep 10s
	done
}


##############
#### loop through 50 CPUS and call serial loop in background
#############
for i in $starteach
do
	j=$(($i+$each-1))
	echo $i $j
	run_LEAreps $i $j
	sleep 2s
done

#echo "Running R LEA"
#for i in $(seq $start $finish)
#do
#    echo $i
#    Rscript --vanilla ../src/Proc_Sims_LEA.R ${i} 'Invers' > ${i}"_Invers_R_LEAonlyR.out" 2> ${i}"_Invers_R_LEAonlyR.error" & echo $!
#    sleep 10s
#done
