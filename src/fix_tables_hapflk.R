library(vcfR)

first <- 10921
last  <- 10961

rawDataDir <- "/home/kevin/LOTTERHOS_LAB/UsefulScripts/FreemanScripts/hapflk_analysis/results/PRUNED"
indDir <- "/home/kevin/LOTTERHOS_LAB/UsefulScripts/data/TTT_RecombinationGenomeScans/results_final"

add1 <- function(x){
  y <- x+1
  return(y)
}

for (i in first:last){
  allChr <- paste(rawDataDir,"/", i, "_", "allChr.hapflk", sep="")
  print(allChr)
  if(!file.exists(allChr)){
    print(paste(sep = " ", "Skipping sim number", i, "file not found:", allChr))
    next()
  }
  position_index <- scan(paste(sep = "", indDir, "/", i, "_Invers_VCFallFILT_indexes_remaining.txt"))
  position_index <- sapply(position_index, add1)
  
  pruned_table <- read.table(allChr, header = TRUE)
  pruned_table$rs <- position_index
  
  write.table(pruned_table, allChr, row.names = FALSE)
  print(paste(sep =" ", "Created", allChr))
}

