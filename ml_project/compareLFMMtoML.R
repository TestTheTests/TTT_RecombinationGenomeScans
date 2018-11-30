#################################################################
#
# File    : compareLFMMtoML.R
# History : 11-8-2018 Created by Kevin Freeman (KF)
#
################################################################
#
# This R script reads in the sim results in order to 
# evaluate how well the LFMM ridge method is identifying
# QTLs. It corrects the p values using the benjamini-hotchberg
# method and then calculates precision, recall, and F1
# for the corrected values, considering a positive to be 
# q < 0.05 and a negative to be q >= 0.05
#
###############################################################

# ---------------------------------------------------------------
#  calcPrecRecallF1: a function to calculate precision/recall/F1 score
# ---------------------------------------------------------------
#  Inputs : a dataframe containing lfmm q values
#  Outputs: none 
#  Prints a string describing the sensitivity/specifity/F1 score
#  Assuming true positives are labeled as MT=2 and a positive q 
#  value is < 0.05
#---------------------------------------------------------------

calcPrecRecallF1 <- function(qValDf, mutColName = "muttype"){
  truePositives <- qValDf[as.integer(qValDf[[mutColName]]) == 2 & !is.na(qValDf$ridgeQvals) & qValDf$ridgeQvals < 0.05, "pos"]
  falsePositives <- qValDf[as.integer(qValDf[[mutColName]]) != 2 & !is.na(qValDf$ridgeQvals) & qValDf$ridgeQvals < 0.05, "pos"]

  ### calc true negative
  trueNegatives <- qValDf[as.integer(qValDf[[mutColName]]) != 2 & !is.na(qValDf$ridgeQvals) & qValDf$ridgeQvals >= 0.05, "pos"]
  falseNegatives <- qValDf[as.integer(qValDf[[mutColName]]) == 2 & !is.na(qValDf$ridgeQvals) & qValDf$ridgeQvals >= 0.05, "pos"]
  
  ### true positive rate/sensitivity/recall
  
  truePosRate <- length(truePositives)/(length(truePositives) + length(falseNegatives))
  
  ## true negative rate/specificity
  
  trueNegRate <- length(trueNegatives)/(length(trueNegatives) + length(falsePositives))
  
  ## precision/positive predictive value
  
  precision <- length(truePositives)/(length(truePositives) + length(falsePositives))
  
  ## F1 (harmonic mean of precision and sensitivity)
  
  f1 <- 2 * precision * truePosRate / (precision + truePosRate)
  print(paste("precision:", precision,  "   recall:", truePosRate, "    F1 score: ", f1))
  print(paste("true positives:", length(truePositives), "    false positives:", length(falsePositives),
              "    false negatives:", length(falseNegatives), "    true negatives:", length(trueNegatives)))
}

# -----------------------------------------------------------------
# relabel: small function to relabel the SNPs close to a given pos
# -----------------------------------------------------------------
# Inputs: row of matrix, df
# Output: what the linked label should be (MT=2 if there is a 
# QTN w/in 200bp, the same label as before if not)
# -----------------------------------------------------------------
relabel <- function(row, fullDf){
  first <- as.integer(row["pos"]) - 200
  last   <- as.integer(row["pos"]) + 200
  linkedMuts <- fullDf[fullDf$pos <= last & fullDf$pos >= first, "muttype"]
  if ('MT=2' %in% linkedMuts){
    return("MT=2")
  }
  else {
    return(row["linkedLabel"])
  }
}

######## Read in all the sim results

first = TRUE
for (sim in 10900:10961){
  tryCatch({
              scanResults <- read.csv(paste0("/media/kevin/TOSHIBA_EXT/TTT_RecombinationGenomeScans/results_final/", sim, "_Invers_ScanResults.txt"), 
                        header = TRUE, sep = " ")
              ridgeLog10p <- scanResults$LFMM_ridge_0.0_ALL_log10p
              ridgePvals <- 10^(-ridgeLog10p)
    
              ridgeQvals <- p.adjust(ridgePvals, method = "hochberg")
      
              scanResults <- cbind(scanResults, ridgeQvals)
              
              ### Add extra column that labels linked QTNs
              scanResults$linkedLabel <- scanResults$muttype
              scanResults$linkedLabel <- apply(scanResults, 1, relabel, fullDf = scanResults)
              scanResults$linkedLabel <- as.factor(scanResults$linkedLabel)

              if (first){
                allResults <- scanResults
                first <- FALSE
              }
              else{
                allResults <- rbind(allResults, scanResults)
              }
              }, error = function(e){
                message(paste("file does not exitst:", sim, "skipping"))
                return(NULL)
              },
              warning = function(w) {
                message(paste("Sim", sim, "caused a warning: ", w))
                return(NULL)
              },
              finally={
                message(paste("processed sim: ", sim))
              })
  
}

######### Calc sensitivity/Specificity for anything labeled as a QTL
print("For all QTLs:")
calcPrecRecallF1(qValDf = allResults, mutCol = "muttype")

##### calc sensitiviy/specificity for just small QTLs
smallQTLs <- allResults
smallQTLs$prop[is.na(smallQTLs$prop)] <- 0
smallQTLs <- smallQTLs[smallQTLs$prop < .20, ]

print("For just small QTLs:")
calcPrecRecallF1(smallQTLs, "muttype")

#### calc for just large QTLs

lgQTLs <- allResults
#lgQTLs$prop[(is.na(lgQTLs$prop))] <- 0
lgQTLs <- lgQTLs[lgQTLs$prop >= .20 | is.na(lgQTLs$prop), ]

print("for just large QTLs:")
calcPrecRecallF1(lgQTLs, "muttype")

print("for all QTLs and linked mutations: ")
calcPrecRecallF1(allResults, "linkedLabel")

calcPrecRecallF1(smallQTLs, "linkedLabel")
calcPrecRecallF1(lgQTLs, "linkedLabel")
