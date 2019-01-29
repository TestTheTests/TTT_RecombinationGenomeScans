#!/bin/bash

#  Script.sh
#  
#
#  Created by katielotterhos on 1/4/19.
#  

awk -F' ' '{print NF; exit}' results_final_old2/10900_Invers_ScanResults.txt


#### FIX 1 ###################
results_final_old2/10900_Invers_ScanResults.txt

for i in {10900..11109}
do
    echo $i

awk '{$44=$45=$46=$47=$48=$49=$50=""; print $0}' results_final_old2/${i}_Invers_ScanResults.txt > results_final/${i}_Invers_ScanResults.txt

awk '{print $44, $45, $46, $47, $48, $49, $50}' results_final_old2/${i}_Invers_ScanResults.txt > results_final/${i}_Invers_scikitallel.txt

done


#### FIX 2 ###################
mkdir results_final_old3
mv results_final/*Invers_ScanResults.txt results_final_old3
# remove white space
for i in {10900..11109}
do
echo $i
sed 's/ *$//' results_final_old3/${i}_Invers_ScanResults.txt > results_final/${i}_Invers_ScanResults.txt
done



sed -e 's/^[ \t]*//' results_final_old3/10900_Invers_scikitallel.txt > results_final/10900_Invers_scikitallel.txt

for i in {10900..11109}
do
echo $i
sed -e 's/^[ \t]*//' results_final_old3/${i}_Invers_scikitallel.txt > results_final/${i}_Invers_scikitallel.txt
done


#### FIX 3 ###################
### Had to remove phenotype columns from Baypass file because these
### were actually run with the environment and not the phenotype

cd ~/Desktop/Repos/TestTheTests/TTT_RecombinationGenomeScans/results_baypass_old
for i in {10900..11109}
do
echo $i
awk '{print $1,$3,$4,$6,$7,$8, $10}' ${i}_baypass_results.txt > ../results_final/${i}_baypass_results.txt
done
