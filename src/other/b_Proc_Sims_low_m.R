### KE Lotterhos
### 20180602
### This script is use to process the raw output 
###### from the LOW RECOMBINATION simulations. It produces
###### a table (for each replication simulation)
###### with test statistics and P-values (in columns)
###### for each locus (in rows). 
### The column name format is <method>_<version>_<condition>_<stat>
### The script does one simulation at a time, 
###### so it can be pseudo-parallelized on a cluster.
### Filtering: loci with MAF < 0.01 are removed
### Script assumes libraries are already installed
### This file is in the "/src" folder
### Output written to "/results_final" folder
### Then concatenated into one large file

### load libraries ####
library(OutFLANK)
library(vcfR)
library(adegenet)
library(pcadapt)
library(lfmm)
library(LEA) 
library(bigsnpr)
library(bigstatsr)
#library(Relatedness)
#library(Demerelate)
library(related)
#library(RColorBrewer)
#library(ggplot2)
#library(fields)
library(rehh)

sessionInfo()

### input arguments ####
args <- commandArgs(TRUE)
seed <- args[1]
type <- args[2] #"Invers"
date()

### read files ####
folder <- "results/"
# this script will be called from working directory in this folder
vcf <- read.vcfR(paste( seed, "_", type, ".vcf.gz", sep=""))
ind0 <- read.table(paste( seed,  "_", type,"_outputIndAll.txt", sep=""), header=TRUE)
muts <- read.table(paste(seed, "_", type,"_outputMuts.txt", sep=""), header=TRUE)
phen_env <-  read.table(paste( seed, "_", type,"_outputPhenEnv.txt", sep=""), header=TRUE)
sim <- system(paste("grep recomb ",  seed, "_", type,"_outputSim.txt", sep=""), TRUE)

### Mutation table ####
muts$simID <- seed
#muts$pos <- muts$position+1

### Fix "CHROM" in VCF ####
# Output "CHROM" in VCF from SliM does not show where the 
# recombination breakpoints are.
# Also output is unordered
### Replace chromsome 1 with actual chromosome positions ####
ends=c(0,  seq(50000, 450000, by=50000))
# linkage groups recombination breakpoints 0.5
# hard-coded from simulations
dim(vcf@gt)
vcf@fix[,"CHROM"] <- NA
POS <- as.numeric(vcf@fix[,"POS"])
for (i in 1:(length(ends)-1)){
  cond <- POS >= ends[i] &  POS < ends[i+1]
  print(c(ends[i], ends[i+1], sum(cond)))
  vcf@fix[cond,"CHROM"] = i
}
#table(vcf@fix[,"CHROM"])
my_ord <- order(as.numeric(vcf@fix[,"POS"]))
vcf2 <- vcf
vcf2 <- vcf[my_ord,]


# The inversion tracking mutation had these weird genotypes
invloc <- which(vcf2@fix[,"POS"]==320000)
vcf2@gt[invloc ,]
levels(as.factor(vcf2@gt[invloc ,]))
vcf2@gt[invloc ,] <- sub("2|2", "1|1",vcf2@gt[invloc ,], fixed=TRUE)
vcf2@gt[invloc ,] <- sub("0|2", "0|1",vcf2@gt[invloc ,], fixed=TRUE)
vcf2@gt[invloc ,] <- sub("2|0", "1|0",vcf2@gt[invloc ,], fixed=TRUE)
#vcf2@gt[invloc ,] <- sub("3|3", "1|1", vcf2@gt[invloc ,], fixed=TRUE)
#vcf2@gt[invloc ,] <- sub("2|3", "0|1", vcf2@gt[invloc ,], fixed=TRUE)
#vcf2@gt[invloc ,] <- sub("3|2", "1|0", vcf2@gt[invloc ,], fixed=TRUE)
vcf <- vcf2
vcf@gt[invloc ,]


### Assign individuals to populations ####
toclust <- ind0[,c("x","y")]
d <- dist(toclust)
hc <- hclust(d, method="ward.D")
k <- 17
group <- cutree(hc, k=k)
table(group)
# label populations by their environment
ind0$group <- group
plot(ind0$x, ind0$y, col=group)
text(ind0$x, ind0$y, group)
# the id should put them back into order
ind <- ind0
head(ind)



### Convert VCF to 012 format ####
# Remove fixed loci (or heterozygote in every individual)
# Remove loci with He <0.01
geno <- vcf@gt[,-1] 
#grep("3", geno)

# Character matrix containing the genotypes
# individuals in columns
# Remove 1st column, which is
position <- getPOS(vcf) # Positions in bp
chromosome <- getCHROM(vcf) # Chromosome information
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2


a_freq <- rowSums(G)/(2*ncol(G))
length(a_freq)

G[invloc,]
a_freq[invloc]
tail(muts)
tapply(group, group, length)
print("allele freq of inversion across populations")
tapply(G[invloc,], group, FUN = function (x){sum(x)/(2*length(x))})


final_df0 <- data.frame(position=position, chrom=chromosome, 
                        a_freq_old = a_freq, 
                        vcf_ord = 1:length(position),
                        unique = paste(position, a_freq, sep="_")
                        #muttype = unlist(lapply(strsplit(getINFO(vcf),split = ";"),function(x){x[6]}))
)


muts$unique <- paste(muts$position, muts$freq, sep="_")
# some mutations may be at the same site, so we match here based on frequency
head(final_df0)
head(muts)

#merge(final_df0, muts, by="unique", all.y=TRUE )

final_df02 <- merge(final_df0, muts, by="unique", all.x=TRUE )
final_df02[!is.na(final_df02$descrip),]

### Mutation information ###

if (length( which(duplicated(final_df0$unique))) > 0){
  print("Error: duplicated unique locus names")
  break;
} # sanity check, should be 0


### Calculate Relatedness ####

# total # derived alleles / total # individuals
#rem = c(which(rowSums(G)==0), which(rowSums(G-2)==0)) ## fixed loci
keep1 <- which(final_df0$a_freq_old > 0.05 & final_df0$a_freq_old< 0.95)

set.seed(99999)
# seed set here to make sure we subsample the same 800 loci for relatedness
# Can be variation in relatedness values and filtering based on 
keep2 <- sort(sample(keep1, 800, replace = FALSE))

final_df0$for_relatedness <- FALSE
final_df0$for_relatedness[keep2] <- TRUE

### organize data for relatedness ####
genosub <- geno[keep2,]
genosub <- genosub
G_relate <- matrix(NA, nrow=ncol(genosub), ncol=nrow(genosub)*2)
abob <- t(substr(genosub,start = 1, stop=1))
bbob <- t(substr(genosub,start = 3, stop=3))
odd <- seq(1, ncol(G_relate), by=2)
G_relate[,odd] <- abob
G_relate[,odd+1] <- bbob
rownames(G_relate) <- rownames(abob)
colnames(G_relate) <- rep(position[keep2], each=2)

class(G_relate) <- "numeric"
G_relate[1:5,1:10]
t(genosub[1:5,1:5])

G_relate <- data.frame(Individual = rownames(G_relate), population=ind$group, G_relate)
head(G_relate[,1:20])
#Loci.test(G_relate, bt=1000, ref.pop=NA,object=TRUE, value="rxy", file.output=TRUE)

#Demerelate package
#pop.relate <- Emp.calc(G_relate, value="rxy",ref.pop="NA")
#pop.relate

#related package
G_relate2 <- G_relate[,-2]
#str(GenotypeData)
G_relate2$Individual <- as.character(G_relate$Individual)
G_relate2[,2:ncol(G_relate2)] <- G_relate2[,2:ncol(G_relate2)]+1 # 0s read as missing data
head(G_relate2[,1:11])
#whichlocs <- sample(2:ncol(G_relate2), 100, replace = FALSE)
poprelate2 <- coancestry(G_relate2 ,  lynchrd =1 )
#read.table("output.txt", header=TRUE)
str(poprelate2)
head(poprelate2$relatedness, 50)
#boxplot(poprelate2$relatedness[,6:10], las=2)

#names of highly related individuals
(rem_ind <- unique(poprelate2$relatedness$ind2.id[which(poprelate2$relatedness$lynchrd>0.5)]))
# indexes of unrelated individuals
ind_keep <- which(!(G_relate2$Individual %in% rem_ind))
length(ind_keep)

# specify which individuals are kept for analysis
head(ind)
ind$infinal <- FALSE
ind$infinal[ind_keep] <- TRUE
cbind(table(ind$group), table(ind$group[ind$infinal]))
hist(ind$phenotype1)
hist(ind$phenotype1[ind$infinal])


### Phenotype-environment correlation
print("phenotype-environment correlation")
(phen_env_corr <- cor.test(ind$phenotype1[ind$infinal], ind$envi[ind$infinal]))
#qplot(ind$x, ind$y, colour=ind$phenotype1)
#qplot(ind$x, ind$y, colour=ind$envi)

# First, remove the related individuals and re-calculate allele frequency
if (length(ind_keep)>0){
  G_sub <- G[,ind_keep]
}else{
  G_sub <- G
}
a_freq2 <- rowSums(G_sub)/(2*ncol(G_sub))
final_df0$a_freq_final <- a_freq2

plot(final_df0$a_freq_old, final_df0$a_freq_final, xlim=c(0,1))

# Next, decide which loci to filter based on the re-calculated allele frequency 
keep_loci = which(a_freq2 > 0.01 & a_freq2 < 0.99)
final_df0$keep_loci <- FALSE
final_df0$keep_loci[keep_loci] <- TRUE
sum(is.na(final_df0$keep_loci))

G_sub2 <- G_sub[keep_loci,]

vcf_filt <- vcf[keep_loci,c(1, ind_keep+1)]
# first column has FORMAT

dim(G_sub2)
dim(vcf_filt)
muts$freq_old <- muts$freq

muts$freq_final <- final_df0$a_freq_final[match(muts$unique, final_df0$unique)]
identical(muts$freq_final, final_df0$a_freq_final[match(muts$unique, final_df0$unique)])

muts <- merge(muts, final_df0, all.x=TRUE, all.y=FALSE)

# Update mutation table to have updated allele frequencies after related individuals removed
muts$pa2[muts$type=="m2"] <- muts$selCoef[muts$type=="m2"]^2*
  (muts$freq_final[muts$type=="m2"])*
  (1-muts$freq_final[muts$type=="m2"])
muts$prop <- round(muts$pa2/sum(muts$pa2, na.rm=TRUE),5)
muts$He <- 2*muts$freq_final*(1-muts$freq_final)
muts[order(muts$He),]
#plot(muts$He[muts$He>0.005], abs(muts$selCoef[muts$He>0.005]))
#text(x = muts$He[muts$He>0.005], 
#     y = abs(muts$selCoef[muts$He>0.005]), 
#     labels = muts$position[muts$He>0.005], cex=0.5)
#sum(muts$prop[muts$count==TRUE])

head(final_df)

### Obtain a quasi-independent set of loci ####
Sys.sleep(20) # this helped offset multiple processes on the cluster when
# calling many replicates in the background
training <- list(G = G_sub2, position = final_df0$pos[which(final_df0$keep_loci)], 
                 chromosome = final_df0$chrom[which(final_df0$keep_loci)])
options(bigstatsr.typecast.warning = FALSE)
G_coded <- add_code256(big_copy(t(training$G),
                                type="raw"),
                       code=bigsnpr:::CODE_012)
# puts it in the raw format and stores likelihood genotype probability
newpc<-snp_autoSVD(G=G_coded,
                   infos.chr = as.integer(training$chromosome),
                   infos.pos = training$position)
# this is doing SNP pruning - removing correlated SNPs
# take snps with highest MAF and correlate snps around it
# Snps with R^2 > 0.2 are removed
# the subset is the indexes of the remaining SNPs
# common error: TridiagEigen: failed to compute all the eigenvalues
# missing values?
# sum(is.na(training$G))
# columns with no variation (zero or very low std)
# hist(apply(training$G, 1, sd), xlim=c(0,1))
# I have noticed that this error arises when sd < 0.1
print("Yoyo")    
#str(newpc)
#plot(newpc)
which_pruned =  attr(newpc, which="subset")
# which_pruned are the indexes of the quasi-independent 
# set of loci that are kept after pruning for LD
training$G_coded <- G_coded
training$G_pruned <- training$G[which_pruned,]#
training$which_pruned <- which_pruned
rm(which_pruned)
#str(training)

final_df <- final_df0[which(final_df0$keep_loci),]
head(final_df)
sum(is.na(final_df0$keep_loci))
sum(is.na(final_df0$a_freq_final))

final_df$quasi_indep <- FALSE
final_df$quasi_indep[training$which_pruned] <- TRUE

### PCA loadings and scores if all data is used ####
gename <- paste0(seed, "_genotypes.lfmm")
envname <- paste0( seed, "_gradients.env")
write.lfmm(t(training$G), gename)
write.env(ind$envi[ind_keep],envname)

### PCADAPT VERSION 3.0.4
Sys.sleep(20) # it appears R on the cluster continues to run before files are written
pcafile <- read.pcadapt(gename, type="lfmm")
#Sys.sleep(60) # it appears R on the cluster continues to run before files are written
pca_all <- pcadapt(pcafile,K=3)#, data.type="genotype")
#Sys.sleep(60)
head(pca_all$loadings)
final_df$pca_ALL_PC1_loadings <- pca_all$loadings[,1]
final_df$pca_ALL_PC2_loadings <- pca_all$loadings[,2]
final_df$pca_ALL_PC3_loadings <- pca_all$loadings[,3]
ind$pca_ALL_PC1_scores <- NA
ind$pca_ALL_PC2_scores <- NA
ind$pca_ALL_PC1_scores[ind_keep] <- pca_all$scores[,1]
ind$pca_ALL_PC2_scores[ind_keep] <- pca_all$scores[,2]
head(final_df)
head(ind)

#plot(final_df$pos, final_df$pca_ALL_PC1_loadings)
#plot(final_df$pos, final_df$pca_ALL_PC2_loadings)    
#plot(final_df$pos, final_df$pca_ALL_PC3_loadings) 

#qplot(ind$pca_ALL_PC1_scores, ind$pca_ALL_PC2_scores, colour=ind$envi, 
#      main="Individual scores without LD pruning") #+ 
#qplot(ind$pca_ALL_PC1_scores, ind$pca_ALL_PC2_scores, colour=ind$phenotype1, 
#      main="Individual scores without LD pruning")


### PCA loadings if pruned data is used ####
gename2 <- paste0(seed, "_genotypes_pruned.lfmm")
write.lfmm(t(training$G_pruned), gename2)
pcafile2 <- read.pcadapt(gename2, type="lfmm")
pca_pruned <- pcadapt(pcafile2,K=3)

final_df$pca_PRUNED_PC1_loadings <- NA # prevents error in next 3 lines of code
final_df$pca_PRUNED_PC2_loadings <- NA
final_df$pca_PRUNED_PC1_loadings[training$which_pruned] <- pca_pruned$loadings[,1]
final_df$pca_PRUNED_PC2_loadings[training$which_pruned] <- pca_pruned$loadings[,2]  

ind$pca_PRUNED_PC1_scores <-NA
ind$pca_PRUNED_PC2_scores <- NA
ind$pca_PRUNED_PC1_scores[ind_keep] <-  pca_pruned$scores[,1] 
ind$pca_PRUNED_PC2_scores[ind_keep] <-  pca_pruned$scores[,2] 

cor.test(final_df$pca_ALL_PC1_loadings, final_df$pca_PRUNED_PC1_loadings)
cor.test(final_df$pca_ALL_PC2_loadings, final_df$pca_PRUNED_PC2_loadings)
#plot(final_df$pos, final_df$pca_PRUNED_PC1_loadings)
#plot(final_df$pos, final_df$pca_PRUNED_PC2_loadings)

#qplot(ind$pca_PRUNED_PC1_scores, ind$pca_PRUNED_PC2_scores, colour=ind$envi, 
#      main="Individual scores WITH LD pruning") #+ 
#qplot(ind$pca_PRUNED_PC1_scores, ind$pca_PRUNED_PC2_scores, colour=ind$phenotype1, 
#      main="Individual scores WITH LD pruning")

#plot(ind$pca_ALL_PC1_scores, ind$envi)
#  abline(lm(ind$envi~ind$pca_ALL_PC1_scores))
#plot(ind$pca_PRUNED_PC1_scores, ind$envi)
#  abline(lm(ind$envi~ind$pca_PRUNED_PC1_scores))

### PCADAPT 3.0.4 all data ####
final_df$pcadapt_3.0.4_ALL_chisq <- as.numeric(pca_all$chi2.stat)
final_df$pcadapt_3.0.4_ALL_log10p <- -log10(pca_all$pvalues)
#plot(final_df$pos, final_df$pcadapt_3.0.4_ALL_chisq)
#plot(final_df$pos, final_df$pcadapt_3.0.4_ALL_log10p)
### PCADAPT pruned data ####    
test <- snp_gc(snp_pcadapt(training$G_coded, U.row = newpc$u[,1]))
final_df$pcadapt_3.0.4_PRUNED_log10p <- -predict(test,log10=T)
#plot(final_df$pos, final_df$pcadapt_3.0.4_PRUNED_log10p )
cor.test(final_df$pcadapt_3.0.4_ALL_log10p, final_df$pcadapt_3.0.4_PRUNED_log10p, method = "spearman")
#plot(final_df$pos,final_df$pcadapt_3.1.0_PRUNED_log10p)

print("StartOutFLANK")
### OutFLANK v0.2 all data ####
FstDataFrame <- MakeDiploidFSTMat(t(training$G),final_df$vcf_ord,
                                  ind$group[ind_keep])
out_ini <- OutFLANK(FstDataFrame, NumberOfSamples=k) 
str(out_ini)
#OutFLANKResultsPlotter(out_ini, withOutliers = TRUE,
#                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
#                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
#OutFLANKResultsPlotter(out_ini, withOutliers = TRUE,
#                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
#                         TRUE, RightZoomFraction = 0.15, titletext = NULL)
#hist(out_ini$results$pvaluesRightTail)
#plot(final_df$position, final_df$OutFLANK_0.2_ALL_log10p)

### OutFLANK pruned data ####    
out_pruned <- OutFLANK(FstDataFrame[training$which_pruned,], NumberOfSamples=k)     
str(out_pruned)

P1 <- pOutlierFinderChiSqNoCorr(FstDataFrame, 
                                Fstbar = out_pruned$FSTNoCorrbar, 
                                dfInferred = out_pruned$dfInferred, Hmin=0.1)
P1 <- P1[order(P1$LocusName),]
identical(P1$LocusName, final_df$vcf_ord) # should be TRUE
identical(out_ini$results$LocusName, final_df$vcf_ord) # should be TRUE
forfinal <- data.frame(vcf_ord = out_ini$results$LocusName,
                       OutFLANK_0.2_FST = out_ini$results$FST,
                       OutFLANK_0.2_He = out_ini$results$He,
                       OutFLANK_0.2_ALL_log10p = -log10(out_ini$results$pvaluesRightTail),
                       OutFLANK_0.2_PRUNED_log10p = -log10(P1$pvaluesRightTail))
final_df_temp <- merge(final_df, forfinal, by="vcf_ord", all.x=TRUE)
dim(final_df)
dim(final_df_temp)
head(final_df_temp)
final_df <- final_df_temp
#plot(final_df$pos, final_df$OutFLANK_0.2_He)
# plot(final_df$pos, final_df$OutFLANK_0.2_ALL_log10p) 
#plot(final_df$pos, final_df$OutFLANK_0.2_PRUNED_log10p) 

#plot(final_df$OutFLANK_0.2_ALL_log10p, final_df$OutFLANK_0.2_PRUNED_log10p)

head(final_df)

### LFMM 0.0 Ridge regression model genotype ~ phenotype ####
# extract scaled genotypes
scaled.genotype <- scale(as.matrix(t(training$G)))
#scaled.genotype <- as.matrix(t(sim1$G))
# extract scaled phenotypes
phen <- scale(as.matrix(ind$phenotype1[ind_keep]))
# ridge regression
lfmm.ridge <- lfmm::lfmm_ridge(Y = scaled.genotype, X = phen, K = 3, lambda = 1e-4)
#The lfmm.ridge object contains estimates for the latent variables and for the effect sizes. Here, the estimates are used for computing calibrated significance values and for testing associations between the response matrix Y and the explanatory variable x. It can be done as follows:
lfmm.test.ridge <- lfmm::lfmm_test(Y = scaled.genotype, X = phen, lfmm = lfmm.ridge, calibrate = "gif")
p.values.ridge <- lfmm.test.ridge$calibrated.pvalue
print(c("lfmm.test.ridge$gif", lfmm.test.ridge$gif))
#hist(p.values.ridge, col = "lightgreen", main="LFMM ridge")
#qval <- qvalue::qvalue(p.values)
#plot(qval)
#The plot suggests that setting fdr.level = 0.025 warrant few false positives.
#qval <- qvalue::qvalue(p.values, fdr.level = 0.005)
#candidates <- which(qval$significant)
#plot(training$position, -log10(p.values.ridge), cex = .5, pch = 19, col = "black", main="LFMM ridge", ylim=c(0, 60))
#plot_layers(y_head=55, y_arrows=c(10, 0))

final_df$LFMM_ridge_0.0_ALL_log10p <- -log10(as.numeric(p.values.ridge))
#plot(final_df$pos, final_df$LFMM_ridge_0.0_ALL_log10p) 

### LFMM LASSO model model genotype ~ phenotype ####    
#LFMM parameters can alternatively be estimated by solving regularized least-squares mimimization, with lasso penalty as follows.
lfmm.lasso <- lfmm::lfmm_lasso(Y = scaled.genotype, X = phen, K = 3, nozero.prop = 0.02)
# The lfmm.lasso object contains new estimates for the latent variables and for the effect sizes. 
# Note that for lasso, we didn't set the value of a regularization parameter. 
# Instead, we set the proportion of non-null effects (here 2 percent).
lfmm.test.lasso <- lfmm::lfmm_test(Y = scaled.genotype, X = phen, lfmm = lfmm.lasso, calibrate = "gif")
p.values.lasso <- lfmm.test.lasso$calibrated.pvalue
print(c("lfmm.test.lasso$gif", lfmm.test.lasso$gif))
#hist(p.values.lasso, col = "lightblue")
#plot(training$position, -log10(p.values.lasso), cex = .5, pch = 19, col = "black", main="LFMM lasso", xaxs="i")
#plot_layers(y_head=45, y_arrows=c(5,0))

final_df$LFMM_lasso_0.0_ALL_log10p <- -log10(as.numeric(p.values.lasso))
#plot(final_df$pos, final_df$LFMM_lasso_0.0_ALL_log10p)

### iHS 2.0.2 convert ####
# (nlociperchrom <- table(vcf_filt@fix[,"CHROM"]))
# 
# substring("0|1", 1, last=1)
# substring("0|1", 3, last=3)
# 
# head(vcf_filt)
# keepihs <- which(!(duplicated(vcf_filt@fix[,"POS"])))
# vcf_filt2 = vcf_filt[keepihs,]
# 
# ### Get into right format
# for (i in 1:length(nlociperchrom)){
#   keep <- which((vcf_filt2@fix[,"CHROM"]==i))
#   #head(vcf_filt2@gt[keep,1:10], 10)
#   hap1 <- apply(vcf_filt2@gt[keep,-1], 1, FUN=function(x) substring(x,1,1))
#   #dim(hap1)
#   #head(hap1[,1:10])
#   hap2 <-  apply(vcf_filt2@gt[keep,-1], 1, FUN=function(x) substring(x,3,3))
#   #head(hap2[,1:10])
#   
#   nind <- nrow(ind[ind_keep,])
#   hapt_out <- matrix(NA, nrow=2*nind, ncol=length(keep))
#   odd <- seq(1,(2*nind), by=2)
#   even <- odd +1
#   hapt_out[odd,] <- as.numeric(hap1)
#   rownames(hapt_out) <- rep("", nrow(hapt_out))
#   rownames(hapt_out)[odd] <- rownames(hap1)
#   rownames(hapt_out)[even] <- rownames(hap2)
#   hapt_out[even,] <- as.numeric(hap2)
#   #head(hapt_out[,1:10])
#   write.table(cbind(rownames(hapt_out), hapt_out+1), paste0(seed,"_chrom",i,".thap"), row.names=F, col.names=F, quote = FALSE)
#   Sys.sleep(20) # it appears R on the cluster continues to run before files are written
# }
# #a<- read.table("chrom1.thap")
# 
# ### Also need to convert map.inp
# #Each line contains five columns corresponding to:
# #1. the SNP name
# #2. the SNP chromosome (or scaffold) of origin
# #3. the SNP position on the chromosome (or scaffold). Note that it is up to the user to choose either
# #physical or genetic map positions (if available).
# #4. the SNP ancestral allele (as coded in the haplotype input file)
# #5. the SNP derived alleles (as coded in the haplotype input file)
# map <- data.frame(name=1:nrow(vcf_filt2), chrom=as.numeric(vcf_filt2@fix[,"CHROM"]), 
#                   pos=as.numeric(vcf_filt2@fix[,"POS"]), anc=1, derived=2)
# # setting anc=0 and derived = 1 thinks missing data
# head(map)
# which(duplicated(map$pos))
# write.table(map, paste0(seed,"_map.inp"), row.names=F, col.names=F)
# 
# ### iHS run ####
# cnt=0
# wgscan <- NULL
# for(i in 1:length(nlociperchrom)){
#   cnt=cnt+1
#   tmp.hapfile=paste0(seed,"_chrom",i,".thap")
#   
#   tmp.hap=data2haplohh(hap_file=tmp.hapfile, map_file=paste0(seed,"_map.inp"), chr.name=i,haplotype.in.columns=FALSE)
#   
#   tmp.scan=scan_hh(tmp.hap,threads=4)
#   
#   if(cnt==1){wgscan=tmp.scan}else{wgscan=rbind(wgscan,tmp.scan)}
# }
# 
# head(wgscan)
# tail(wgscan)
# dim(wgscan)
# dim(final_df)
# 
# ihs=ihh2ihs(wgscan,minmaf=0.05,freqbin=0.05)
# head(ihs$iHS,25)
# colnames(ihs$iHS)[4] <- "rehh_2.0.2_ALL_log10p"
# colnames(ihs$iHS)[3] <- "rehh_2.0.2_ALL_iHS"
# #tail(ihs$iHS)
# ihs$frequency.class
# 
# #plot(ihs$iHS$POSITION,ihs$iHS$`-log10(p-value)`, col=rgb(0,0,0,0.2), pch=20, main="REHH iHS")
# #plot(ihs$iHS$POSITION,ihs$iHS$iHS, col=rgb(0,0,0,0.2), pch=20, main="REHH iHS")
# 
# final_dfc <- final_df[keepihs,] # remove the duplicated loci for merger
# final_dfb <- merge(final_dfc, ihs$iHS, by.x = "pos", by.y="POSITION", all.x = TRUE)
# dim(final_df)
# dim(final_dfb)
# head(final_dfb)
# 
# final_dfd <- merge(final_df, final_dfb, all.x=TRUE)
# dim(final_dfd)
# dim(final_df)
# head(final_dfd)
# final_df <- final_dfd
# head(final_df)  
# all(final_df$vcf_ord == cummax(final_df$vcf_ord))
# check if vcf_ord is monotonically increasing

### Raw associations with environment ####
spearmans <- apply(training$G, 1, FUN=function(x){
  a <- cor.test(x, ind$envi[ind_keep], method = "spearman", )
  return(a$estimate)
})
length(spearmans)
dim(final_df)
sum(is.na(spearmans))
#hist(spearmans)
#plot(training$position, abs(spearmans), cex = .5, pch = 19, col = rgb(0,0,0,0.2), main="raw association", xaxs="i")

final_df$Spearmans_ALL_rho <- spearmans
#plot(final_df$pos, abs(final_df$Spearmans_ALL_rho))
head(final_df)

### Tables to write to file ####
ind$simID <- seed

head(final_df)
head(muts)

print(c("dim final_df", dim(final_df)))
print(c("dim muts", dim(muts)))

final_df$simID <- seed

muts <- muts[,-which(colnames(muts) %in% c("freq", "freq_old", "freq_final"))]

write.vcf(vcf_filt, file = paste("../results_final/", seed, "_", type, "_VCFallFILT.vcf.gz", sep=""))
write.table(ind, file = paste("../results_final/", seed, "_", type, "_indFILT.txt", sep=""), row.names=FALSE)
write.table(final_df, file = paste("../results_final/", seed, "_", type, "_ScanResults.txt", sep=""), row.names = FALSE)
write.table(muts, file = paste("../results_final/", seed, "_", type, "_muts.txt", sep=""), row.names = FALSE)



