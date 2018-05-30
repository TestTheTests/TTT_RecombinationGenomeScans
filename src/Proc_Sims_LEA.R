### Bayesian (LEA) 1.2.0 model genotype ~ environment ####
# Creation of the genotype file: "genotypes.lfmm"
# 400 SNPs for 50 individuals.

#plot(pc$sdev[1:10]^2, pch = 19)
library(LEA)

### input arguments ####
args <- commandArgs(TRUE)
seed <- args[1]
type <- args[2] #"RecomLowReg" or "Invers"
date()

final_df <- read.table(paste("../results_final/", seed, "_", 
                             type, "_ScanResults.txt", sep=""),
                       header=TRUE
                      )

gename <- paste0(seed, "_genotypes.lfmm")
envname <- paste0(seed, "_gradients.env")

################
# running lfmm #
################
project1 = snmf(gename,
               K = 1:10, 
               entropy = TRUE, 
               repetitions = 3,
               project = "new")

# plot cross-entropy criterion of all runs of the project
plot(project1, cex = 1.2, col = "lightblue", pch = 19)

project = lfmm(gename, 
               envname, 
               K = 3, 
               repetitions = 3, 
               project = "new")
#project = load.lfmmProject(paste0(sub(".lfmm","",gename),"_",seed,"_gradients.lfmmProject"))
# summary of the project
summary(project)
# get adjusted p-values using all runs
#pv1 <-  adjusted.pvalues(project, K = 1)
#pv2 <-  adjusted.pvalues(project, K = 2)
#pv3 <-  adjusted.pvalues(project, K = 3)
# above 3 lines don't work on version on cluster
# bioconductor loads different version on my laptop than on cluster
# maybe has to do with version of R?
#pv1 <-  rowMeans(p.values(project, K = 1))
#pv2 <-  rowMeans(p.values(project, K = 2))
pv3 <-  rowMeans(p.values(project, K = 3))
# get the z-scores
#z_ave1 <- rowMeans(z.scores(project, K = 1))
#z_ave2 <- rowMeans(z.scores(project, K = 2))
z_ave3 <- rowMeans(z.scores(project, K = 3))

#plot(training$position, -log10(pv1$p.values), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA",  xaxs="i")
#plot(training$position, -log10(pv2$p.values), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA", xaxs="i")
#plot(training$position, -log10(pv3$p.values), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA", xaxs="i")

#plot(training$position, abs(z_ave1), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA", xaxs="i")
#plot(training$position, abs(z_ave2), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA", xaxs="i")
#plot(training$position, abs(z_ave3), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="LFMM GEA", xaxs="i")

final_df_LEA <- data.frame(#LEA_1.2.0_ALL_K1_log10p =  -log10(pv1),
                       #LEA_1.2.0_ALL_K2_log10p =  -log10(pv2),
                       LEA_1.2.0_ALL_K3_log10p =  -log10(pv3),
                       #LEA_1.2.0_ALL_K1_z = z_ave1,
                       #LEA_1.2.0_ALL_K2_z = z_ave2,
                       LEA_1.2.0_ALL_K3_z = z_ave3)


write.table(final_df_LEA, file = paste("../results_final/", seed, "_", type, "_Results_LEA.txt", sep=""), row.names = FALSE)
