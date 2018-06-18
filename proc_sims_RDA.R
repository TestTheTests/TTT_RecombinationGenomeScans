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
myDir  <- args[3]             ### Directory that the files to be analyzed are contained in 
pruned <- args[4]             ### Have non-quasi-indep alleles been removed? if so, pruned <- "_PRUNED". If not, 
                                # pruned must be set to "" (not NA, not UNDEF, not 0....)
date()

### function to convert 0 based to 1 based numbering ###
add1 <- function(x){
  y <- x+1
  return(y)
}

### read files ####

vcf <- read.vcfR(paste(myDir,"/", seed, "_",
                        pruned, type, "_VCFallFILT.vcf.gz", sep=""))
env <- read.table(paste(myDir,"/", seed, "_", type, "_indFILT.txt", sep = ""), header = TRUE)

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

sim.rda <- rda(rdaG~., data = pred, scale = T)
# get loadings from rda
load.rda <- as.data.frame(scores(sim.rda, choices = c(1:4), display="species"), row.names = FALSE)

###### Process and output loading data ################
# Get indexes of snps, this will be used to combine data later
if (pruned == "PRUNED"){                                       # use output of pruneSNPs.pl
  position_index <- scan(paste(sep = "", myDir, "/", seed, "_Invers_VCFallFILT_indexes_remaining.txt"))
  position_index <- sapply(position_index, add1)
}else {                                                       # table with all snps, can just indexes of dataframe
  position_index <- as.numeric(rownames(load.rda)) 
}

# add position and position_index to loadings table
load.rda <- cbind(position, load.rda)
load.rda <- cbind(position_index, load.rda)

## export data frame of loadings
outfile <- paste(sep="", myDir, "/results/raw_data/", seed, "_", pruned, type, "_loadings_raw.txt")
write.table(load.rda, outfile, row.names = FALSE)
print(paste(sep =" ", "Created", outfile))

# ## plot
# pdf(paste(sep = "", myDir,"/plots/", seed, "_", pruned, type, "_plot.pdf"))
# plot(sim.rda, type="n", scaling=3)
# points(sim.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
# text(sim.rda, scaling=3, display="bp", col="#0868ac", cex=1) 
# title(paste(sep="", seed, "_", pruned, type))
# dev.off()

