#!/bin/bash
set -e
set -u
set -o pipefail

# This script only works if you are able to run ncpus in the background at once

mypath="/home/k.lotterhos/k_store/TestTheTests/TTT_RecombinationGenomeScans"
#mypath="/Users/katie/Desktop/TestTheTests/TTT_RecombinationGenomeScans"
cd $mypath
ncpus=62
start=10900
finish=$(($start + $ncpus-1))
echo $start $finish

# make a list of numbers to seed the sims

cd results

##############
#### run R script LEA
#############

echo "Running R LEA"
for i in $(seq $start $finish)
do
    echo $i
    Rscript --vanilla ../src/Proc_Sims_LEA.R ${i} 'Invers' > ${i}"_Invers_R_LEAonlyR.out" 2> ${i}"_Invers_R_LEAonlyR.error" & echo $!
    sleep 10s
done
