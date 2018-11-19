#!/usr/bin/env python
# coding: utf-8

# In[79]:


#!/usr/bin/python
import re, subprocess, os.path, time, argparse, pkg_resources, sys
import numpy as np, pandas as pd
import allel
from allel.model.ndarray import SortedIndex
from allel.util import asarray_ndim
from multiprocessing import Pool, cpu_count

#############################################################################################
# File    : calcAdditionalStatsFromVCF.py
# History : 10-17-2018 - Created by Kevin Freeman (KF)
#           11-19-2018 - updated to use multiprocessing and break windows at chroms
#--------------------------------------------------------------------------------------------
# This is a script that uses scikit-allel to calculate several summary statistics on a set
# of simulations, stored in VCFs. The only required argument is --dataDir [path to vcf files]
# but there are 2 optional args: --ncpus and --winodow. Run with the -help option to
# see full usage information.
#
# Stats the scripts calculates (updated 11-19):
#  - Tajima's D
#  - pi (sequence diversity)
#  - watterson's theta
#  - Garud's H12
#  - Garud's H2/H1
#  - IHS (integrated haplotype score)
#
###############################################################################################

### small function to print to stderr instead of stdout
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# In[65]:


#-----------------------------------------------------------------------------------------
# calcAndAppendStatValForScan:
#
# The calculations for each statistic. Takes the name of the statistic and the relevant 
# data structures and performs the right calculations. Appends the result to statVals
# ----------------------------------------------------------------------------------------
#
# snpLocs : a numpy array with the positions in bp for each snp in the window
# statName: which stat to calculate 
# statVals: dictionary with stat names as the keys and calculated stats for every snp in
#           order as the values
# subWinStart : start of the subwindow (bp)
# subWinEnd   : end of the subwindow (bp)
# alleleCounts: numpy array representing the # of calls of each allele per variant 
#
#----------------------------------------------------------------------------------------

def calcAndAppendStatValForScan(snpLocs, statName, hapsInSubWin, statVals,
                                subWinStart = None, subWinEnd = None, alleleCounts = None):
# modified code from https://github.com/kern-lab/diploSHIC/blob/master/fvTools.py
    if statName == "tajD":
        statVals[statName].append(allel.stats.diversity.tajima_d(
            alleleCounts, pos=snpLocs, start=subWinStart, stop=subWinEnd))
    elif statName == "pi":
        statVals[statName].append(allel.stats.diversity.sequence_diversity(
            snpLocs, alleleCounts, start=subWinStart, stop=subWinEnd))
    elif statName == "thetaW":
        statVals[statName].append(allel.stats.diversity.watterson_theta(
            snpLocs, alleleCounts, start=subWinStart, stop=subWinEnd)) 
    elif statName == 'H2/H1' or statName == 'H12':
        h1,h12,h123,h21 = allel.stats.selection.garud_h(hapsInSubWin)
        if statName == 'H2/H1':
            statVals[statName].append(h21)
        else:
            statVals[statName].append(h12)
    elif statName == 'ihs':
        try:
            statVals[statName].extend(allel.ihs(hapsInSubWin, snpLocs, map_pos = None, use_threads = False))
        except: 
            print("hapsInSubWin: " + str(len(hapsInSubWin)) + "  snpLocs: " + str(len(snpLocs)))
    else:
        print(statName + " not found")


# In[66]:


#------------------------------------------------------------------------------------------------------
# processFile:
# 
# Where the data processing gets done. Takes a single vcf file, divides it into chromosomes and
# then into windows, and sends it to calcAndAppendStatValForScan, which does the actual calculations
# It then takes the results of these calculations and prints them to a new output file
# -----------------------------------------------------------------------------------------------------
# Inputs  : the name of a vcf file (string)
# Outputs : a new file called "XXXX_Invers_ScanResults_with_sk-allel.txt", where XXXX is the sim number
#           The file contains any stats that were previously calculated in XXXX_Invers_ScanResults.txt
#           along with the new stats calculated in sci-kit allel
#------------------------------------------------------------------------------------------------------

def processFile(f):
    eprint("Processing " + f)
    vcf     = allel.read_vcf(f, fields = ["CHROM", "POS", "GT"])
    m       = re.search('[0-9]{5}', f)      # m is a regex match object
    simNum  = m.group(0)                    # sim number is the first set of 5 numbers in the vcf name
    
    g           = allel.GenotypeArray(vcf["calldata/GT"])
    ac          = g.count_alleles()
    haps        = g.to_haplotypes()
    posSeries   = pd.Series(vcf['variants/CHROM'], index = vcf['variants/POS'])
    # cast chroms to int so we get numeric sort instead of alphabetic sort in groupby
    groupedPos = posSeries.groupby(by = [int(c) for c in vcf['variants/CHROM']])
    
    statvals = { "tajD"   : [],
                 "pi"     : [],
                 "thetaW" : [],
                 "H12"    : [],
                 "H2/H1"  : [],
                 "ihs"    : [] }
    
    for chrom, group in groupedPos:             # only consider 1 chrom at a time to preserve breaks
        pos = list(group.index)
        chromStart = np.searchsorted(vcf['variants/POS'], pos[0])
        chromEnd   = np.searchsorted(vcf['variants/POS'], pos[-1], side = 'right')
        chromHap   = haps.subset(list(range(chromStart, chromEnd)))
        
        # now calc ihs, it can just take the whole chromosome without need for a sliding window
        calcAndAppendStatValForScan(snpLocs = list(pos), statName = 'ihs', statVals = statvals,
                                    hapsInSubWin = chromHap)
        # calc the rest of the stats over a 1000bp sliding window
        for SNP in pos:
            winStart = SNP - int(winsize/2)
            winEnd   = SNP + int(winsize/2)
        
            if winStart < pos[0]:
                winStart = pos[0]
            if winEnd > pos[-1]:
                winEnd = pos[-1]
        
            ####### subset the haplotype array ##############
            # search sorted finds the place in a sorted array where a number can be inserted
            startInd = np.searchsorted(pos, winStart)                   
            endInd   = np.searchsorted(pos, winEnd)
        
            hapsInSubWin = chromHap.subset(list(range(startInd, endInd + 1)))
    
            
            for statName in [key for key in statvals.keys() if key !='ihs']:
                calcAndAppendStatValForScan(list(pos), statName, hapsInSubWin, statvals, winStart, winEnd, ac)
        
    
    scanResultsFile = datadir + simNum + "_Invers_ScanResults.txt"
    outfile         = datadir + simNum + "_Invers_ScanResults_with_sk-allel.txt"

    try:
        newStatDf         = pd.DataFrame.from_dict(statvals)
    except:
        eprint(["length of stat: " + stat + " = " + len(statvals[stat]) for stat in list(statvals.keys())])
        sys.exit()
    # add version to stat column names
    newStatDf.columns = [stat + "_sk-allel_v" + skVersion for stat in list(statvals.keys())] 
    
    allFeatures = []
    try:
        with open (scanResultsFile, 'r') as f:
            for line in f:
                features = line.split()
                features = [f.strip('"') for f in features]
                allFeatures.append(features)
            stats  = allFeatures[1:]
            header = allFeatures[0]
            featDf =  pd.DataFrame(stats, columns = header)
            featDf = featDf[featDf.keep_loci == "TRUE"].reset_index(drop = True)
        featDf = featDf.join(newStatDf)
    except: 
        featDf = newStatDf
        
    featDf.to_csv(outfile, sep = " ", index = False)


# In[67]:


############## Main ############################################################
if __name__ == '__main__':
    
    start_time = time.time()
    
    ######## Parse arguments ######

    parser = argparse.ArgumentParser(
        description = 'Calculate summary statistics using scikit allel on the input VCFs',
        usage       = "calcAdditionalStatsFromVCF.py -d path/directory -n 4 -w 1000")
    parser.add_argument("-d", "--datadir", required = True,
                    help = "full path to the directory with the VCFs")
    parser.add_argument("-n", "--ncpus", default = cpu_count(), type = int,
                    help = "number of cpus to use (default: all of them)")  
    parser.add_argument("-w", "--window", default = 1000, type = int,
                    help = "size (in bp) of sliding window around each SNP for calculations (default: 1000)")

    args     = vars(parser.parse_args())
    datadir  = args['datadir'] 
    nproc    = args['ncpus']
    winsize  = args['window']

    #### generate file list, find sk-a version ################
    skVersion = pkg_resources.get_distribution("scikit-allel").version # for labeling output
    cmd      = 'ls ' + datadir + '*.vcf.gz'
    fileList =  subprocess.run(cmd, shell = True, stdout=subprocess.PIPE).stdout.decode('utf-8')
    fileList = fileList.split("\n")
    fileList = list(filter(None, fileList))            # remove blank entries
    
    ####  analyze vcfs in fileList in parallel ###########
    with Pool(processes = nproc) as pool:             # launch parallel processes
         pool.map(processFile, fileList)
    secs = time.time() - start_time
    
    eprint("\nIt took {} hrs {} minutes {} seconds to process {} simulation files".format(secs // 360, 
                                                                                   secs % 360 // 60,
                                                                                   round(secs % 360 % 60, 2), 
                                                                                   len(fileList)))

