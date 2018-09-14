#################################################################
#
# File: combine_RDA_tables.R
#
#################################################################
#
# combine_RDA_tables.R takes a number of arguments that point it
# towards a set of results from an RDA analysis on pruned and
# unpruned data sets. It combines the pruned and unpruned results
# into one table based on SNPs, so the SNPs at
# each position can be easily compared between pruned and unpruned.
# It is meant to be run on the output of proc_sims_RDA.R
# 
#################################################################


### input arguments ####
# args    <- commandArgs(TRUE)
start   <-  args[1]             # first sim in set
end     <-  args[2]             # last sim in set
myDir   <-  args[3]             # directory /results/ is in, which contains the tables to combine
simType <-  args[4]             # ex: 'Invers'

for (sim in start:end){
  print(sim)
  # skip sim if the corresponding files do not exist
  unpruned.file <- paste(myDir,"/RDA_envi_results/raw_data/", sim, "_", simType, "_loadings_raw.txt", sep="")
  if(!file.exists(unpruned.file)){
    print(paste(sep = " ", "Skipping sim number", sim, "file not found:", unpruned.file))
    next()
  }
  pruned.file <- paste(myDir,"/RDA_envi_results/raw_data/", sim, "_PRUNED_", simType, "_loadings_raw.txt", sep="")
  if(!file.exists(pruned.file)){
    print(paste(sep = " ", "Skipping sim number", sim, "file not found:", pruned.file))
    next()
  }
  
  # unpruned sim
  load.rda <- read.table(unpruned.file, header = TRUE)
  load.rda <- as.data.frame(load.rda)
  colnames(load.rda) <- c("position_index", "position", 
                          "RDAvegan_v2.5.2_ALL_loading_RDA1_envi", 
                          "RDAvegan_v2.5.2_ALL_loading_PCA1_envi", 
                          "RDAvegan_v2.5.2_ALL_loading_PCA2_envi", 
                          "RDAvegan_v2.5.2_ALL_loading_PCA3_envi")

  # pruned sim
  load.rda.pruned <- read.table(pruned.file, header = TRUE)
  load.rda.pruned <- as.data.frame(load.rda.pruned)
  load.rda.pruned <- subset(load.rda.pruned, select=-c(position))
  colnames(load.rda.pruned) <- c("position_index", 
                                 "RDAvegan_v2.5.2_PRUNED_loading_RDA1_envi",
                                 "RDAvegan_v2.5.2_PRUNED_loading_PCA1_envi", 
                                 "RDAvegan_v2.5.2_PRUNED_loading_PCA2_envi", 
                                 "RDAvegan_v2.5.2_PRUNED_loading_PCA3_envi")

  #combine data frames
  combinedDF <- merge(load.rda, load.rda.pruned, by= 'position_index', all = T)
  write.table(combinedDF, paste(sep="",
                                myDir, 
                                "/RDA_envi_results/", 
                                sim, 
                                "_", 
                                simType, 
                                "_RDA_loadings_envi.txt"), 
              row.names = FALSE)
}




