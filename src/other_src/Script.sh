#!/bin/bash

#  Script.sh
#  
#
#  Created by katielotterhos on 1/4/19.
#  

awk -F' ' '{print NF; exit}' results_final_old2/10900_Invers_ScanResults.txt

results_final_old2/10900_Invers_ScanResults.txt

for i in {10900..11109}
do
    echo $i

awk '{$44=$45=$46=$47=$48=$49=$50=""; print $0}' results_final_old2/${i}_Invers_ScanResults.txt > results_final/${i}_Invers_ScanResults.txt

awk '{print $44, $45, $46, $47, $48, $49, $50}' results_final_old2/${i}_Invers_ScanResults.txt > results_final/${i}_Invers_scikitallel.txt

done
