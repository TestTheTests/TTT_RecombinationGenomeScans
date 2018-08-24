args   <- commandArgs(TRUE)
inFile <- args[1]
outDir <- args[2]
simNum <- args[3]
allResults <- data.frame(rs = numeric(), chr = numeric(), pos = numeric(),
                         hapflk = numeric())
for (i in 1:10){
  currentDf <- read.table(paste(inFile, i, ".hapflk", sep = ""), header = TRUE)
  allResults <- rbind(allResults, currentDf)
}

write.table(allResults, paste(sep="", outDir, simNum, "_allChr.hapflk"), sep = " ", quote=FALSE, row.names = FALSE)


