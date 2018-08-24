#!/usr/bin/bash
############################################################################
#
# File   : get_p_values.sh
# History: Created by Kevin Freeman (KF) 7/10/2018
#
# Description : bash script to run scaling_chi2_hapflk.py on all hapflk
#               files in the current directory -- allows p-values to be 
#               calculated easily
#
# 
# 
# Variables that may need to be changed in different environments:
#   -- K  : number of clusters used in the hapflk analysis
#   -- N  : number of groups (populations) in the data set
#   -- src: absolute path to the scaling_chi2_hapflk.py script 
#
###########################################################################

K=15
N=39
src="/home/kevin/LOTTERHOS_LAB/UsefulScripts/FreemanScripts/run_hapflk/src"
mapfile -t hapflk_files_arr < <( ls *.hapflk )
arr_length=${#hapflk_files_arr[@]}

for i in $(seq 0 $arr_length)
do
	hapflk_file=${hapflk_files_arr[$i]}
	python2.7 ${src}/scaling_chi2_hapflk.py $hapflk_file $K $N
done
