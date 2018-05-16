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
#### run R script
#############
#Rscript ../../Installs.r > Install_output.txt
#Rscript TO DO DELETE THIS LINE

echo "Running R scripts"
for i in $(seq $start $finish)
do
    echo $i
    Rscript --vanilla ../src/Proc_Sims.R ${i} 'Invers' > ${i}"_Invers_R_onlyR.out" 2> ${i}"_Invers_R_onlyR.error" & echo $!
    sleep 1m
done

## Some of the R simulations failed the first time around, which I think happened because R was trying to read files that weren't completely written yet
## I added some wait times to the code, but all simulations worked the second time around
##
