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
library(vcfR)

########## collect and check command args ##########################
args    <- commandArgs(TRUE)
if (length(args) == 0) {
  print('Using default values')
  args[1] <- 10900
  args[2] <- 10961
  args[3] <- "/home/kevin/LOTTERHOS_LAB/UsefulScripts/FreemanScripts/run_hapflk"
  args[4] <- "Invers"
  args[5] <- "/home/kevin/LOTTERHOS_LAB/UsefulScripts/data/TTT_RecombinationGenomeScans/results_final"
} else if (length(args) != 4){
  stop("Wrong number of arguments. This script requires 4 arguments:
       1) start (first sim)                   ---> default: 10900
       2) end   (last sim)                    ---> default: 10961
       3) directory containing results tables ---> default: ~/LOTTERHOS_LAB/UsefulScripts/FreemanScripts/run_hapflk
       4) type of sim                         ---> default: Invers
       5) directory with indexes remaining    ---> default: /home/kevin/LOTTERHOS_LAB/UsefulScripts/data/TTT_RecombinationGenomeScans/results_final

      To use defaults, simply provide 0 arguments")
}
start   <-  args[1]             # first sim in set
end     <-  args[2]             # last sim in set
myDir   <-  args[3]             # directory /results/ is in, which contains the tables to combine
simType <-  args[4]             # ex: 'Invers'
indDir  <-  args[5]             # directory where "indexes remaining" files are, telling script how to combine pruned and unpruned tables

# indDir is directory that contains files indicating which indexes remain after pruning (tab-delimited text files, one for each sim)


add1 <- function(x){
  y <- x+1
  return(y)
}

############## Loop through specified files ####################################################
vers <- "hapflk_v1.4"
for (sim in start:end){
  print(sim)
  # skip sim if the corresponding files do not exist
  vcf <- paste(myDir,"/", sim, "_", simType, "_VCFallFILT.vcf.gz", sep="")
  if(!file.exists(vcf)){
    print(paste(sep = " ", "Skipping sim number", sim, "file not found:", vcf))
    next()
  }
  ################### Add indexes pruned table, ids SNPs for combining ################################
  allChr <- paste(myDir,"/results/PRUNED/", sim, "_", "allChr.hapflk", sep="")
  position_index <- scan(paste(sep = "", indDir, "/", sim, "_Invers_VCFallFILT_indexes_remaining.txt"))
  position_index <- sapply(position_index, add1)
  
  pruned_table <- read.table(allChr, header = TRUE)
  pruned_table$rs <- position_index
  
  write.table(pruned_table, allChr, row.names = FALSE)
  print(paste(sep =" ", "Created", allChr))
  
  ################### Read all 3 tables into memory ##################################################
  # unpruned sim
  ALL <- read.table(paste(myDir, "/results/ALL/", sim, "_allChr.hapflk_sc", sep = ""), header = TRUE)
  ALL <- as.data.frame(ALL)
  colnames(ALL) <- c("rs", "chr", "pos", 
                     paste(sep ="", vers, "_ALL"), 
                     paste(sep="", vers, "_ALL_scaled"), 
                     paste(sep = "", vers, "_ALL_log_pvalue"))
  # -log10(pvalues)
  ALL$hapflk_v1.4_ALL_log_pvalue <- -log10(ALL$hapflk_v1.4_ALL_log_pvalue)
  
  # pruned sim
  PRUNED <- read.table(paste(myDir, "/results/PRUNED/", sim, "_allChr.hapflk_sc", sep = ""), header = TRUE)
  PRUNED <- as.data.frame(PRUNED)
  colnames(PRUNED) <- c("rs", "chr", "pos", 
                        paste(sep ="", vers, "_PRUNED"), 
                        paste(sep ="", vers, "_PRUNED_scaled"), 
                        paste(sep ="", vers, "_PRUNED_log_pvalue"))
  PRUNED <- PRUNED[, !(names(PRUNED) %in% c("chr", "pos"))]            # drop pos and col, would be dups with ALL df
  # -log10(pvalues)
  PRUNED$hapflk_v1.4_PRUNED_log_pvalue <- -log10(PRUNED$hapflk_v1.4_PRUNED_log_pvalue)
  
  # all data, pruned kinship sim
  ALL_PRUNED_KIN <- read.table(paste(myDir, "/results/ALL_PRUNED_KIN/", sim, "_allChr.hapflk_sc", sep = ""), header = TRUE)
  ALL_PRUNED_KIN <- as.data.frame(ALL_PRUNED_KIN)
  colnames(ALL_PRUNED_KIN) <- c("rs", "chr", "pos", 
                                paste(sep = "", vers, "_PRUNED_KIN"),
                                paste(sep = "", vers, "_PRUNED_KIN_scaled"),
                                paste(sep = "", vers, "_PRUNED_KIN_log_pvalue"))
  ALL_PRUNED_KIN <- ALL_PRUNED_KIN[, !(names(ALL_PRUNED_KIN) %in% c("chr", "pos"))]  # drop pos and col, would be dups with ALL df
  # -log10(pval)
  ALL_PRUNED_KIN$hapflk_v1.4_PRUNED_KIN_log_pvalue <- -log10(ALL_PRUNED_KIN$hapflk_v1.4_PRUNED_KIN_log_pvalue)
  
  ##############  Combine data frames and print file ################################################
  combinedDF <- merge(ALL, PRUNED, by = 'rs', all = T)
  combinedDF <- merge(combinedDF, ALL_PRUNED_KIN, by = 'rs', all = T)
  write.table(combinedDF, paste(sep="",
                                myDir, 
                                "/results/", 
                                sim, 
                                "_", 
                                simType, 
                                "_hapflk_scores_sc.txt"), 
              row.names = FALSE)
  print(paste(sep = "", "Created ",myDir, "/results/",sim, "_", simType, "_hapflk_scores_sc.txt"))
}




