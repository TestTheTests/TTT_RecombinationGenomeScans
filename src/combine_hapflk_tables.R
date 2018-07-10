#################################################################
#
# File: combine_hapflk_tables.R
#
#################################################################
#
# combine_hapflk_tables.R takes a number of arguments that point it
# towards a set of results from an RDA analysis on pruned and
# unpruned data sets. It combines the pruned, unpruned, and pruned
# matrix results into one table based on the SNPs so 
# that SNPs can be easily compared between pruned and unpruned.
# It is meant to be run on the output of run_hapflk.sh
# 
#################################################################


### input arguments ####
# args    <- commandArgs(TRUE)
start   <-  10900 #args[1]             # first sim in set
end     <-  10961 #args[2]             # last sim in set
myDir   <-  "/home/kevin/LOTTERHOS_LAB/UsefulScripts/FreemanScripts/hapflk_analysis" #args[3]             # directory /results/ is in, which contains the tables to combine
simType <-  "Invers"                           # args[4]             # ex: 'Invers'

for (sim in start:end){
  print(sim)
  # skip sim if the corresponding files do not exist
  vcf <- paste(myDir,"/", sim, "_", simType, "_VCFallFILT.vcf.gz", sep="")
  if(!file.exists(vcf)){
    print(paste(sep = " ", "Skipping sim number", sim, "file not found:", vcf))
    next()
  }
  # unpruned sim
  ALL <- read.table(paste(myDir, "/results/ALL/", sim, "_allChr.hapflk", sep = ""), header = TRUE)
  ALL <- as.data.frame(ALL)
  colnames(ALL) <- c("rs", "chr", "pos", "hapflk_v1.4_ALL")
  
  # pruned sim
  PRUNED <- read.table(paste(myDir, "/results/PRUNED/", sim, "_allChr.hapflk", sep = ""), header = TRUE)
  PRUNED <- as.data.frame(PRUNED)
  colnames(PRUNED) <- c("rs", "chr", "pos","hapflk_v1.4_PRUNED")
  PRUNED <- PRUNED[, !(names(PRUNED) %in% c("chr", "pos"))]            # drop pos and col, would be dups with ALL df
  
  # all data, pruned kinship sim
  ALL_PRUNED_KIN <- read.table(paste(myDir, "/results/ALL_PRUNED_KIN/", sim, "_allChr.hapflk", sep = ""), header = TRUE)
  ALL_PRUNED_KIN <- as.data.frame(ALL_PRUNED_KIN)
  colnames(ALL_PRUNED_KIN) <- c("rs", "chr", "pos","hapflk_v1.4_ALL_PRUNED_KIN")
  ALL_PRUNED_KIN <- ALL_PRUNED_KIN[, !(names(ALL_PRUNED_KIN) %in% c("chr", "pos"))]  # drop pos and col, would be dups with ALL df
  
  #combine data frames
  combinedDF <- merge(ALL, PRUNED, by = 'rs', all = T)
  combinedDF <- merge(combinedDF, ALL_PRUNED_KIN, by = 'rs', all = T)
  write.table(combinedDF, paste(sep="",
                                myDir, 
                                "/results/", 
                                sim, 
                                "_", 
                                simType, 
                                "_hapflk_scores.txt"), 
              row.names = FALSE)
}




