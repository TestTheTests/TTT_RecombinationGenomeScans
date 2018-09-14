#!/usr/bin/python3
# coding: utf-8

import subprocess
import numpy as np
import pandas as pd
import numpy.ma as ma
import os.path
import math
import glob
import sys
import warnings
from decimal import *

######## readFeatures #################

# readFeatures is a subroutine that takes the name of a file and a list of the indices (columns) for the metrics that
# you want to use as features. Indices are 0 based. Files must be tab or space delimited. It returns a list that contains 
# a modified header followed by all the indicated features at each snp



def readFeatures(file, metricIndices):
    allFeatures = []
    with open (file, 'r') as f:
        first = True
        for line in f:
            splitLine = line.split()
            #features = '\t'.join(splitLine[i] for i in metricIndices)
            features  = [ splitLine[i] for i in metricIndices ]
            if first:
                features = [f.strip('"') for f in features]
                first = False
            allFeatures.append(features)
    return pd.DataFrame(np.array(allFeatures[1:]), columns = allFeatures[0])

##### writeToFile

# write to file takes the summaries calculated in summarizeByWindow and writes
# them to the appropriate files

def writeToFile(summaries, outDir, fileNameAndClass, header, nwin):
    filename        = outDir + "/" + fileNameAndClass[0]
    # classLabel      = fileNameAndClass[1]
    
    fileHeaderList  = []
    windows         = list(range(0,nwin))
        
    #### create new header, but only if file doesn't already exist
    if not os.path.exists(filename):                  
        for stat in header:
        # create new header entries, example of format: 'Hscan_v1.3_H12_win0'
            for win in windows:
                fileHeaderList.append(stat + '_win' + str(win))
        fileHeader = "\t".join(fileHeaderList)
        # fileHeader = "classLabel\t" + fileHeader  # add another column for classLabel
        with open(filename, "w") as myfile:
            myfile.write(fileHeader)
            
    #### write data to file
    numStats        = summaries.shape[1]
    line            = []# [classLabel]
    with open(filename, "a") as myfile:
        myfile.write("\n")
        for i in range(0, numStats):
            currStat   = summaries[:,i].tolist()
            line.extend(currStat)
            
        stringLine = "\t".join(map(str, line))
        myfile.write(stringLine)


####### summarizeByWindow ##########################

# This function takes a pandas datatframe, where row in the dataframe corresponds to the statistics for one SNP. It divides this list into n windows, each representing an equal genetic distance along the chromosome. The windows do not necessarily have an equal number of SNPs. The average value for each statistic of interest is calculated over the window. The function returns an n (# of windows) x 9 (# of stats) numpy array

def summarizeByWindow(snps, nwin, scaled):
    
    ## divide the input chromosome into n equally sized windows
    snps  = snps.apply(pd.to_numeric, errors='coerce')
    snps  = snps.replace([np.inf, -np.inf], np.nan)
    
    firstSnpPos = snps['pos'].iloc[0]
    winBegin    = round(firstSnpPos, -3)
    lastSnpPos  = snps['pos'].iloc[-1]
    winEnd      = round(lastSnpPos, -3)
    windowSize  = (winEnd - winBegin) / nwin
    
    ## initialize window
    currWinMin   = winBegin
    currWinMax   = winBegin + windowSize
    winSummaries = []
    
    for i in range(0, nwin):
        maxMask = snps['pos'] <= currWinMax
        minMask = snps['pos'] >= currWinMin
        
        currWinDf = snps[maxMask & minMask]
        
        currWinDf = currWinDf.drop(columns = 'pos')
        winSummaries.append(currWinDf.mean().values)
        
        ### move to new window
        currWinMin = currWinMax + 1
        currWinMax = currWinMin + windowSize

    stack = np.vstack(winSummaries)
    if scaled:
    # calculate the % of total statistic contributed by each window
        scaledStats = []
        for column in stack.T:
            minStat = np.nanmin(column)
            if minStat < 0:
                # if any columns are negative, scale by addition so that min column = 0
                column = column - minStat
            scaledStats.append(np.divide(column, np.nansum(column)))
         # put columns back together and transpose
        scaledStats = np.vstack(scaledStats).T
        return scaledStats
    else:
         return stack

######## string2float #########################

# This is a helper function for use in summarizing the statistics. Some of the statistics contain missing values,
# represented by 'NA' in the data. This function takes a list of strings (statistics for one SNP), converts all 
# the numbers in the list to 'float' and converts any non-numbers ('NA') to numpy.nan


def string2float(dataList):
    for i, x in enumerate(dataList):
        try:
            dataList[i] = float(x)
        except ValueError:
            dataList[i] = np.nan

