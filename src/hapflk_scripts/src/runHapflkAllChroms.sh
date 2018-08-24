#!/usr/bin/bash
chrStart=1
chrFinish=10
ncpu=8


### command line args
inputDir=$1
outputDir=$2
prefix=$3
kinship=$4
K=$5

echo "chr $chrStart to  $chrFinish"
echo "running hapflk with K = $K"
echo

for i in $(seq $chrStart $chrFinish)
do
	hapflk --file "${inputDir}/${prefix}_chr$i" \
		-p ${outputDir}/${prefix}_chr$i --kinship $kinship \
		-K $K --nfit=10 --ncpu=$ncpu
	wait
done

