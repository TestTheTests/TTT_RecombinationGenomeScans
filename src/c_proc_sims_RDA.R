#############################################################################
#
# File: proc_sims_RDA.R
#
#############################################################################
#
# This script uses the 'vegan' package to perform RDA on a set of SNPs to 
# identify those that are likely to be under selection. It takes 4 arguments 
# from the command and outputs a file with the position of each SNP analyzed
# and the first 4 loadings (RDA1,PCA1, PCA2, PCA3) for each SNP.
# It also produces a graph to visualize the effects of RDA1 and RDA2.
#
#############################################################################

### load libraries ####
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA
library(vcfR)

sessionInfo()

### input arguments ####
args   <- commandArgs(TRUE)
seed   <- args[1]             ### Number corresponding to the name of the sim to be analyzed
type   <- args[2]             ### Type of sim, ex: "Invers". Used to identify the file to be read
myDir  <- "../results_final/"            ### Directory that the files to be analyzed are contained in 
date()

### function to convert 0 based to 1 based numbering ###
add1 <- function(x){
  y <- x+1
  return(y)
}

### read files ####

vcf <- read.vcfR(paste(myDir, seed, "_", type, "_VCFallFILT.vcf.gz", sep=""))
env <- read.table(paste(myDir,"/", seed, "_", type, "_indFILT.txt", sep = ""), header = TRUE)
final_df <- read.table(paste("../results_final/", seed, "_", 
                             type, "_ScanResults.txt", sep=""),
                       header=TRUE)

### Convert VCF to 012 format ####
geno <- vcf@gt[,-1] 
# Character matrix containing the genotypes
# individuals in columns
# Remove 1st column, which is
position <- getPOS(vcf)     # Positions in bp
chromosome <- getCHROM(vcf) # Chromosome information
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

rdaG <- t(G)            # must transpose table to use vegan for RDA

### select envi column
env <- env[env$infinal=="TRUE",]
pred <- subset(env, select=envi)

##### Calculate RDA #####################
if (sum(is.na(rdaG)) != 0){
  stop("Data contains missing values")
}


#### OLD CODE FROM KEVIN #####
# sim.rda <- rda(rdaG~., data = pred, scale = T)
# get loadings from rda
# load.rda <- as.data.frame(scores(sim.rda, choices = c(1:4), display="species"), row.names = FALSE)

###### Process and output loading data ################
# Get indexes of snps, this will be used to combine data later
# if (pruned == "PRUNED"){                                       # use output of pruneSNPs.pl
#  position_index <- scan(paste(sep = "", myDir, "/", seed, "_Invers_VCFallFILT_indexes_remaining.txt"))
#  position_index <- sapply(position_index, add1)
#}else {                                                       # table with all snps, can just indexes of dataframe
#  position_index <- as.numeric(rownames(load.rda)) 
#}

### Run RDA on all data
sim.rda.all <- rda(rdaG~., data = pred, scale = T)

final_df_inG <- final_df[which(final_df$keep_loci),]

if(length(sim.rda.all$colsum)!= nrow(final_df_inG)){
  print("Error: number of loci doesn't match")
}


load.rda.all <- as.data.frame(scores(sim.rda.all, choices = c(1:4), display="species"), row.names = FALSE)
colnames(load.rda.all) <- paste0("RDAvegan_2.52_ALL_", colnames(load.rda.all))
load.rda.all$unique <- final_df_inG$unique

### Run RDA on pruned data
sim.rda.pruned <- rda(rdaG[,which(final_df_inG$quasi_indep)]~., data = pred, scale = T)

load.rda.pruned <- as.data.frame(scores(sim.rda.pruned, choices = c(1:4), display="species"), row.names = FALSE)
colnames(load.rda.pruned) <- paste0("RDAvegan_2.52_PRUNED_", colnames(load.rda.pruned))
load.rda.pruned$unique <- final_df_inG$unique[which(final_df_inG$quasi_indep)]

# merge tables and add to final_df
load.rda <- merge(load.rda.all, load.rda.pruned, by="unique", all.x=TRUE)
dim(load.rda)
dim(final_df)

final_df2 <- merge(final_df, load.rda, all.x=TRUE, by="unique")
head(final_df2)

#plot(final_df2$RDAvegan_2.52_ALL_RDA1, final_df2$RDAvegan_2.52_PRUNED_RDA1)
#  abline(0,1)
#plot(final_df2$RDAvegan_2.52_ALL_PC1)
#plot(final_df2$RDAvegan_2.52_PRUNED_PC1)

## export data frame of loadings
write.table(final_df2, paste("../results_final/", seed, "_", 
                            type, "_ScanResults.txt", sep=""),row.names=FALSE)
print(paste(sep =" ", "Created outfile"))


# ## plot
# pdf(paste(sep = "", myDir,"/plots/", seed, "_", pruned, type, "_plot.pdf"))
# plot(sim.rda, type="n", scaling=3)
# points(sim.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
# text(sim.rda, scaling=3, display="bp", col="#0868ac", cex=1) 
# title(paste(sep="", seed, "_", pruned, type))
# dev.off()

