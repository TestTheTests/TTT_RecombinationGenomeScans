#################################################################
#
# File: combine_RDA_tables.R
#
#################################################################
#
# combine_RDA_tables.R takes a number of arguments that point it
# towards a set of results from an RDA analysis on pruned and
# unpruned data sets. It combines the pruned and unpruned results
# into one table based on the "position" column, so the SNPs at
# each position can be easily compared between pruned and unpruned.
# It is meant to be run on the output of proc_sims_RDA.R
# 
#################################################################


### input arguments ####
args    <- commandArgs(TRUE)
start   <- args[1]             # first sim in set
end     <- args[2]             # last sim in set
myDir   <- args[3]             # directory /results/ is in, which contains the tables to combine
simType <- args[4]             # ex: 'Invers'

for (sim in start:end){
  print(sim)
  # skip sim if the corresponding files do not exist
  unpruned.file <- paste(myDir,"/results/raw_data/", sim, "_", simType, "_loadings_raw.txt", sep="")
  if(!file.exists(unpruned.file)){
    print(paste(sep = " ", "Skipping sim number", sim, "file not found:", unpruned.file))
    next()
  }
  pruned.file <- paste(myDir,"/results/raw_data/", sim, "_PRUNED_", simType, "_loadings_raw.txt", sep="")
  if(!file.exists(pruned.file)){
    print(paste(sep = " ", "Skipping sim number", sim, "file not found:", pruned.file))
    next()
  }
  
  # unpruned sim
  load.rda <- read.table(paste(sep="", myDir, "/results/raw_data/", sim, "_", simType, "_loadings_raw.txt"), header = TRUE)
  load.rda <- as.data.frame(load.rda)
  colnames(load.rda) <- c("position", "RDAvegan_v2.5.2_ALL_loading_RDA1", "RDAvegan_v2.5.2_ALL_loading_RDA2",
                                 "RDAvegan_v2.5.2_ALL_loading_PCA1", "RDAvegan_v2.5.2_ALL_loading_PCA2", 
                                 "RDAvegan_v2.5.2_ALL_loading_PCA3")

  # pruned sim
  load.rda.pruned <- read.table(paste(sep="", myDir, "/results/raw_data/", sim, "_PRUNED_", simType, "_loadings_raw.txt"), header = TRUE)
  load.rda.pruned <- as.data.frame(load.rda.pruned)
  colnames(load.rda.pruned) <- c("position", "RDAvegan_v2.5.2_PRUNED_loading_RDA1", "RDAvegan_v2.5.2_PRUNED_loading_RDA2",
                                 "RDAvegan_v2.5.2_PRUNED_loading_PCA1", "RDAvegan_v2.5.2_PRUNED_loading_PCA2", 
                                 "RDAvegan_v2.5.2_PRUNED_loading_PCA3")

  #combine data frames
  combinedDF <- merge(load.rda.pruned, load.rda, by.x = 'position', all = T)
  write.table(combinedDF, paste(sep="", myDir, "/results/", sim, "_", simType, "_RDA_loadings.txt"))
}




